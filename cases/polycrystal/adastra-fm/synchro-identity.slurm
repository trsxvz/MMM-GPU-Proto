#!/bin/bash
# Copies of the arrays of the materials made by the codes which read them on
# the host (mfem-mgis patches 0302 and 0303, JOURNAL.md, 2026-11-18): the
# executable FM_EXECUTABLE_B gives the outputs of FM_EXECUTABLE_A bit for
# bit, on FM_GRAINS grains (default 600) partitioned for the processes, with
# the options of the measurements of the FM (one material for all the
# grains, line search, prediction by linear extrapolation, NewtonK0 with the
# first choice of Eisenstat and Walker from 0.01, mean stresses computed by
# kernels):
# - on the GPUs (FM_GPU_STEPS steps of 1 s, default 20), the element
#   matrices being added in a fixed order
#   (MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered), which makes a
#   computation on the GPUs reproducible: A is run twice to check it;
# - on the cores (FM_CPU_STEPS steps, default 3), with the assembly of the
#   measurements (FM_CPU_STEPS=0: none).
# Compared: <computation>.res, the averages per grain
# (<computation>/avgStressPolycrystal) and the lines of the iterations of
# the logs. The copies of the arrays of the materials between the host and
# the GPUs, written at the end of the logs, are reported in the journal.
# The GPUs are not reproducible bit for bit with the options of the FM,
# even with the ordered assembly (hypre on the GPUs, see JOURNAL.md,
# 2026-11-18): FM_GPU_OPTIONS adds options to the GPU computations, for
# instance --hypre-on-host, with FM_MESH_GPU (mesh read by all the
# processes and partitioned by the program, --mesh_mode FromScratch, for
# instance a smaller RVE), FM_GPU_PROCESSES (processes per GPU, default 1;
# with 24, the GPU computation has as many processes as the CPU one, which
# makes hypre on the host fast enough) and FM_GPU_RUNS (default "A1 B A2").
# The job fails if a computation fails or if outputs differ; otherwise the
# file <results>/passed is written.
#SBATCH --job-name=fm-synchro-identity
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --time=01:30:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
A=${FM_EXECUTABLE_A:?}
B=${FM_EXECUTABLE_B:?}
n=${FM_GRAINS:-600}
out=$R/synchro-identity-g$n-$SLURM_JOB_ID
mkdir -p "$out"
log "start, $n grains, A $A, B $B, GPU ${FM_GPU_STEPS:-20} steps (ordered" \
  "assembly${FM_GPU_OPTIONS:+, $FM_GPU_OPTIONS}${FM_MESH_GPU:+, $FM_MESH_GPU}," \
  "${FM_GPU_PROCESSES:-1} processes per GPU)," \
  "CPU ${FM_CPU_STEPS:-3} steps, scripts $P"
export vectors=$MESHES/v$n/vectors_${n}crystals.txt refinement=0
export solver_label=best run_timeout_gpu=${FM_TIMEOUT_GPU:-2400} run_timeout_cpu=5400
solver="--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01"
export solver_gpu=$solver solver_cpu=$solver
options="--single-material --line-search Backtracking --prediction LinearExtrapolation"
export MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels
# lines of the iterations of a log
lines() {
  grep -E 'Fixed point iteration|Newton iteration|Solving time step|Solution at time|resolution failed|dividing|integration fail|iterations of the linear solver' "$1.txt"
}
# usage: run label executable mode steps [variables...]
run() {
  local label=$1 e=$2 mode=$3 steps=$4 np mesh extra
  shift 4
  np=$((SLURM_JOB_NUM_NODES * 4 * ${FM_GPU_PROCESSES:-1}))
  [ "$mode" = cpu ] && np=$((SLURM_JOB_NUM_NODES * 96))
  mesh=$W/meshes/part/v$n-np$np/mesh.
  extra="$options --mesh_mode Restart"
  if [ "$mode" = gpu ]; then
    if [ -n "${FM_MESH_GPU:-}" ]; then
      mesh=$FM_MESH_GPU
      extra="$options --mesh_mode FromScratch"
    fi
    extra="$extra ${FM_GPU_OPTIONS:-}"
  fi
  env "$@" executable="$e" label="$label" modes=$mode nsteps="$steps" \
    gpu_processes="${FM_GPU_PROCESSES:-1}" \
    duration="$steps" mesh="$mesh" extra_options="$extra" \
    "$P/cases/measure-polycrystal-grains.sh" "$out" > "$out/$label.out" 2>&1
  log "$label: $(tail -1 "$out/series.txt")"
  if [ "$mode" = gpu ]; then
    log "$label: copies $(grep -A2 'copies of the arrays of the materials' \
      "$out/$label-r0-gpu-n$SLURM_JOB_NUM_NODES-best.txt" | tail -2 | tr -s ' ' | tr '\n' ';')"
  fi
}
differ=0
# usage: compare a b: outputs of two computations
compare() {
  local a=$out/$1 b=$out/$2
  if [ -s "$a.res" ] && cmp "$a.res" "$b.res" > "$out/diff-$2.txt" 2>&1 &&
    cmp "$a/avgStressPolycrystal" "$b/avgStressPolycrystal" >> "$out/diff-$2.txt" 2>&1 &&
    diff <(lines "$a") <(lines "$b") >> "$out/diff-$2.txt" 2>&1; then
    log "$2: same outputs as $1 ($(wc -c < "$a.res") bytes of .res," \
      "$(wc -c < "$a/avgStressPolycrystal") bytes of averages per grain," \
      "$(lines "$a" | wc -l) lines of iterations)"
  else
    log "$2: outputs differ from $1 (see $out/diff-$2.txt)"
    differ=1
  fi
}
steps=${FM_GPU_STEPS:-20}
ordered=MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered
for r in ${FM_GPU_RUNS:-A1 B A2}; do
  e=$A
  [ "${r:0:1}" = B ] && e=$B
  run "gpu-$r" "$e" gpu "$steps" "$ordered"
done
for r in ${FM_GPU_RUNS:-A1 B A2}; do
  [ "$r" = A1 ] || compare gpu-A1-r0-gpu-n1-best "gpu-$r-r0-gpu-n1-best"
done
steps=${FM_CPU_STEPS:-3}
if [ "$steps" -gt 0 ]; then
  run cpu-A "$A" cpu "$steps"
  run cpu-B "$B" cpu "$steps"
  compare cpu-A-r0-cpu-n1-best cpu-B-r0-cpu-n1-best
fi
energy "$out/energy.txt"
if [ $differ = 1 ]; then
  log "end, outputs differ (the next jobs of the chain do not start)"
  exit 1
fi
touch "$out/passed"
log end
