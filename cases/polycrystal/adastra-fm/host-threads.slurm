#!/bin/bash
# Setup of the host shared between the threads of OpenMP (option
# enable-openmp of mfem-mgis, JOURNAL.md, 2026-11-16): the first assembly of
# the tangent operator of the final computation of the polycrystal
# (FM_GRAINS grains on all the GPUs of the allocation, mesh partitioned
# beforehand, FM_NSTEPS steps, 1 by default) with the executable
# FM_EXECUTABLE, for each number of threads per process of FM_THREADS
# (default "1 2 4 8 12 24"; OMP_NUM_THREADS, the processes keeping all
# their cores), FM_REPS times (default 1), and with the executable
# FM_REFERENCE, if given, before each repetition.
#
# The job first writes where the launcher places the processes of a GPU
# computation and of a CPU computation (cores allowed to each process), as
# cases/measure-polycrystal-grains.sh launches them, and the settings of
# OpenMP that the processes see.
#
# A line per computation in <results>/series.txt (see
# cases/measure-polycrystal-grains.sh), and in <results>/threads.txt:
# computation, threads, time of the first assembly (first NS::getJacobian,
# maximum over the processes), of the mesh (Mesh::Constructor) and of the
# whole computation (s), Newton iterations and axial stress (Pa).
# The job fails if a computation fails or does not end its steps, or if an
# axial stress differs from the first one by more than FM_TOLERANCE
# (relative, default 1e-5, see fm-chain/perf-ab.slurm); otherwise the file
# <results>/passed is written.
#SBATCH --job-name=fm-host-threads
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
n=${FM_GRAINS:?}
e=${FM_EXECUTABLE:?}
np=$((SLURM_JOB_NUM_NODES * 4))
out=$R/host-threads-g$n-$SLURM_JOB_ID
mkdir -p "$out"
export mesh=${FM_MESH:-$W/meshes/part/v$n-np$np/mesh.}
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
export refinement=0 nsteps=${FM_NSTEPS:-1} duration=${FM_DURATION:-${FM_NSTEPS:-1}}
export modes=gpu solver_label=best run_timeout_gpu=${FM_TIMEOUT_GPU:-2400}
export solver_gpu=${FM_SOLVER:---non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01}
export extra_options=${FM_EXTRA_OPTIONS:---single-material --line-search Backtracking --prediction LinearExtrapolation --mesh_mode Restart}
log "start, $n grains, $SLURM_JOB_NUM_NODES nodes, $nsteps steps," \
  "threads ${FM_THREADS:-1 2 4 8 12 24} x${FM_REPS:-1}, executable $e," \
  "reference ${FM_REFERENCE:-none}, scripts $P"
# placement: the environment of the computations (scripts/env.sh of the
# deployed prototype, as cases/measure-common.sh)
(
  set +u
  source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null
  set -u
  export PATH=${PROTO_MPI_LAUNCHERS:+$PROTO_MPI_LAUNCHERS:}$PATH
  unset OMP_PROC_BIND OMP_PLACES
  cd "$out" || exit 1
  # shellcheck disable=SC2016
  show='echo "rank $SLURM_PROCID on $(hostname -s): cores $(grep Cpus_allowed_list /proc/self/status | cut -f2), OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}, OMP_PROC_BIND=${OMP_PROC_BIND:-unset}, OMP_PLACES=${OMP_PLACES:-unset}"'
  {
    echo "GPU computation: 4 processes per node, 24 cores each"
    env OMP_NUM_THREADS=24 mpirun --map-by slot:PE=24 --bind-to core \
      -n "$np" bash -c "$show" 2> /dev/null | sort -t' ' -k2 -n
    echo "CPU computation: 96 processes per node, first node"
    env OMP_NUM_THREADS=1 mpirun -n $((SLURM_JOB_NUM_NODES * 96)) \
      bash -c "$show" 2> /dev/null | sort -t' ' -k2 -n | head -96
  } > placement.txt
)
log "placement: $(sed -n 2,5p "$out/placement.txt" | tr '\n' ';')" \
  "CPU: $(grep -c 'OMP_NUM_THREADS=1,' "$out/placement.txt") processes," \
  "$(awk '/^CPU/ {c = 1; next} c {print $6}' "$out/placement.txt" | sort -u | wc -l) different sets of cores"
failed=0
# usage: run_one name executable [threads]
run_one() {
  local name=$1 threads=${3:-} line f
  executable=$2 label=$name \
    gpu_environment=${threads:+OMP_NUM_THREADS=$threads} \
    "$P/cases/measure-polycrystal-grains.sh" "$out" > "$out/$name.out" 2>&1
  line=$(tail -1 "$out/series.txt")
  f=$out/$name-r0-gpu-n$SLURM_JOB_NUM_NODES-best.txt
  # sections: name | calls | minimum | mean | maximum
  section() {
    awk -F'|' -v s="--> $1 " 'index($0, s) {gsub(/ /, "", $(NF - 3)); print $(NF - 3); exit}' "$f" 2> /dev/null
  }
  echo "$name ${threads:--} $(section NS::getJacobian) $(section Mesh::Constructor)" \
    "$(awk '{print $7, $11, $8}' <<< "$line")" >> "$out/threads.txt"
  log "$name: $(tail -1 "$out/threads.txt")"
  [ "$(awk '{print $5, $9}' <<< "$line")" = "0 $nsteps" ] || failed=1
}
echo "# computation threads first_assembly(s) mesh(s) total(s) newton_iterations axial_stress(Pa)" \
  > "$out/threads.txt"
for i in $(seq "${FM_REPS:-1}"); do
  [ -z "${FM_REFERENCE:-}" ] || run_one "reference-rep$i" "$FM_REFERENCE"
  for t in ${FM_THREADS:-1 2 4 8 12 24}; do
    run_one "threads$t-rep$i" "$e" "$t"
  done
done
energy "$out/energy.txt"
if ! awk -v tol="${FM_TOLERANCE:-1e-5}" '$2 == "gpu" {
    if (ref == "") ref = $8
    d = ($8 - ref) / ref; if (d < 0) d = -d
    if (d > tol) { print $1, $8, "differs from", ref; bad = 1 } }
  END { exit bad }' "$out/series.txt" > "$out/stresses.txt" 2>&1; then
  log "axial stresses: $(tr '\n' ';' < "$out/stresses.txt")"
  failed=1
fi
if [ $failed = 0 ]; then
  touch "$out/passed"
  log "end, all the computations ended, same axial stresses within ${FM_TOLERANCE:-1e-5}"
else
  log "end, FAILURES (the next jobs of the chain do not start)"
  exit 1
fi
