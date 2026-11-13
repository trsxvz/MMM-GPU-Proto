#!/bin/bash
# Polycrystal of mm-opera-hpc on a given RVE (mesh and orientations), on all
# the GPUs and then on all the cores of the nodes of the Slurm allocation:
# comparison at equal resources, node for node (convergence study in mesh
# and series in number of grains, see results/ and JOURNAL.md).
#
# usage: cases/measure-polycrystal-grains.sh <results directory>
#
# Variables read (the first two are required):
# - mesh, vectors: mesh (gmsh, periodic) and orientations of the grains
#   (two vectors per grain), absolute paths;
# - label: name of the computations (default: name of the mesh);
# - refinement: uniform refinements of the mesh (default 0);
# - order: order of the finite elements (default 2);
# - nsteps, duration: steps and simulated time in s (default 5 and 5);
# - modes: `gpu` and/or `cpu` (default `gpu cpu`);
# - nodes: nodes used (default SLURM_JOB_NUM_NODES, or 1);
# - gpus_per_node, cores_per_node: default 4 and 96 (MI300A of Adastra);
# - run_timeout_gpu, run_timeout_cpu: time limits of the computations in s
#   (default 2400 and 10800; run_timeout, if given, for both);
# - cpu_environment: variables given to the CPU computation only (for
#   instance MPICH_GPU_SUPPORT_ENABLED=0);
# - solver_gpu, solver_cpu: options of the solvers of the polycrystal on the
#   GPU and on the CPU (default: NewtonK0, Eisenstat-Walker 1 from 0.1, the
#   best on both sides on Jean Zay, results/L5-polycrystal-2026-10-10-0502);
# - profiler: command placed before the executable (for instance
#   `rocprofv3 --kernel-trace ... --`);
# - solver_label: suffix of the names of the computations (comparison of
#   solvers in the same results directory);
# - executable: executable of the polycrystal (default: the one of the build
#   of mm-opera-hpc), for instance built against another mfem-mgis;
# - extra_options: options added to the ones of the polycrystal (for
#   instance --single-material);
# - behaviour_library_gpu, behaviour_library_cpu: libraries of the behaviour
#   (generic-parallel interface on the GPU, generic interface on the CPU)
#   instead of the ones of the builds, for instance variants of the law
#   (polycrystal/law-variants/).
# The GPU computation runs one process per GPU (cores_per_node /
# gpus_per_node cores each) with GPU-aware MPI, the CPU computation one
# process per core. The fixed point on the lateral stresses uses the
# secant acceleration, the averages of the stresses per grain are written
# at each step (post-processings without ParaView export) in
# <results>/<computation>/avgStressPolycrystal. A line per computation is
# appended to <results>/series.txt: label, mode, nodes, processes, status,
# unknowns, total time (s, section root of the profiling: setup included),
# axial stress at the end (Pa), time steps (with the divided ones),
# resolutions of the fixed point on the lateral stresses, Newton iterations,
# failed resolutions, solver label, time of the mesh (s, Mesh::Constructor:
# serial reading, partition, refinements), time of the loop over the steps
# (s, NLEP::solve), energy (J) and duration (s) of the Slurm step of the
# computation (the energy is read right after it, `-` outside Slurm), time
# of the macroscopic stresses after each resolution (s,
# computeMacroscopicCauchyStress, `-` if the executable does not time it),
# time of the post-processings (s, common::post_processing_step), and wall
# clock time from the end of the first step to the end of the last one (s,
# timestamps of <computation>.steps).
set -u
work_package=grains
mesh=${mesh:?mesh (absolute path) required}
vectors=${vectors:?vectors (absolute path) required}
label=${label:-$(basename "$mesh" .msh)}
refinement=${refinement:-0}
order=${order:-2}
nsteps=${nsteps:-5}
duration=${duration:-$nsteps}
nodes=${nodes:-${SLURM_JOB_NUM_NODES:-1}}
gpus_per_node=${gpus_per_node:-4}
cores_per_node=${cores_per_node:-96}
run_timeout_gpu=${run_timeout_gpu:-${run_timeout:-2400}}
run_timeout_cpu=${run_timeout_cpu:-${run_timeout:-10800}}
read -r -a cpu_environment <<< "${cpu_environment:-}"
post_processing=""
source "$(dirname "$0")/measure-common.sh" polycrystal "$@"
# command placed before the executable inside mpirun (`profiler`, for
# instance rocprofv3 and its options followed by --), empty by default
read -r -a launcher <<< "${profiler:-}"
summary="$summary|number of finite element unknowns"
executable=${executable:-$B/uniaxial-polycrystal}
case_options="-m $mesh -f $vectors --order $order -r $refinement"
case_options="$case_options --macroscopic-secant-acceleration"
case_options="$case_options --disable-paraview-export"
case_options="$case_options${extra_options:+ $extra_options}"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
gpu=${PROTO_GPU_DEVICE:-cuda}
gpu_library=$gp/generic-parallel-$gpu/src/libGenericParallelBehaviour-$gpu.so
gpu_library=${behaviour_library_gpu:-$gpu_library}
generic=${behaviour_library_cpu:-$generic}
default_solver="--non-linear-solver NewtonK0"
default_solver="$default_solver --adaptive-linear-solver-tolerance EisenstatWalker1"
default_solver="$default_solver --adaptive-linear-solver-initial-tolerance 0.1"
solver_gpu=${solver_gpu:-$default_solver}
solver_cpu=${solver_cpu:-$default_solver}
suffix=${solver_label:+-$solver_label}
{
  echo "mesh: $mesh, refinement $refinement, order $order"
  echo "steps: $nsteps over $duration s, nodes: $nodes"
  echo "solvers${suffix}: GPU $solver_gpu; CPU $solver_cpu"
} | tee -a "$out/summary.txt"
# usage: watch_steps name
# starts a watcher of the macroscopic stresses of the computation (one line
# per accepted step, written by the process 0 at the end of the step), which
# copies each new line to <results>/<name>.steps preceded by the time (s
# since the epoch, read every second), to give the time of each step (steady
# regime, computations stopped by their time limit); stop_watch name stops
# it after a last reading
watch_steps() {
  local res=$out/$1.res steps=$out/$1.steps stop=$out/$1.stop
  : > "$res"
  : > "$steps"
  rm -f "$stop"
  (
    n=0
    while :; do
      m=$(wc -l < "$res")
      if [ "$m" -gt "$n" ]; then
        t=$(date +%s.%N)
        tail -n +$((n + 1)) "$res" | head -n $((m - n)) | sed "s/^/$t /" \
          >> "$steps"
        n=$m
      fi
      [ -f "$stop" ] && break
      sleep 1
    done
  ) &
  watcher=$!
}
stop_watch() {
  touch "$out/$1.stop"
  wait "$watcher" 2> /dev/null
  rm -f "$out/$1.stop"
}
# usage: record name mode processes
record() {
  local name=$1 mode=$2 np=$3 status unknowns time szz steps resolutions
  local newton failures mesh_time loop_time energy elapsed cauchy_time
  local pp_time wall
  status=missing
  if [ -f "$out/$name.txt" ]; then
    status=$(sed -n 's/^   status //p' "$out/summary.txt" | tail -1)
  fi
  unknowns=$(sed -n 's/.*number of finite element unknowns: *//p' \
    "$out/$name.txt" | head -1)
  time=$(python3 -I "$P/cases/summary-profiling.py" "$out/$name.txt" \
    2> /dev/null | awk -F'|' 'NR > 2 {gsub(/ /, "", $2); sub(/\(.*/, "", $2);
      print $2; exit}')
  szz=$(awk 'NF >= 7 {s = $7} END {print s}' "$out/$name.res" 2> /dev/null)
  steps=$(awk 'NF >= 7' "$out/$name.res" 2> /dev/null | wc -l)
  resolutions=$(grep -c '^Fixed Point iteration' "$out/$name.txt" 2> /dev/null)
  newton=$(grep -c '||r||/||r_0||' "$out/$name.txt" 2> /dev/null)
  failures=$(grep -c 'the resolution failed' "$out/$name.txt" 2> /dev/null)
  # maximum over the processes of a first level section of the profiling
  # (line ` | |--> name | calls | min | average | max | ...`)
  section() {
    awk -F'|' -v s="$1" '$2 == " " && $3 ~ ("^--> " s " *$") {
      gsub(/ /, "", $7); print $7; exit}' "$out/$name.txt" 2> /dev/null
  }
  mesh_time=$(section 'Mesh::Constructor')
  loop_time=$(section 'NLEP::solve')
  cauchy_time=$(section 'computeMacroscopicCauchyStress')
  pp_time=$(section 'common::post_processing_step')
  wall=$(awk 'NR == 1 {t0 = $1} {t = $1} END {if (NR > 1) printf "%.1f", t - t0}' \
    "$out/$name.steps" 2> /dev/null)
  # energy and duration of the last step of the job (this computation),
  # waiting for the accounting (a few seconds)
  energy=- elapsed=-
  if [ -n "${SLURM_JOB_ID:-}" ]; then
    local i line
    for i in $(seq 12); do
      line=$(sacct -j "$SLURM_JOB_ID" -P -n \
        --format=JobID,ConsumedEnergyRaw,ElapsedRaw,State 2> /dev/null |
        awk -F'|' '$1 ~ /\.[0-9]+$/' | sort -t. -k2 -n | tail -1)
      case "$line" in
        *'|COMPLETED'|*'|FAILED'|*'|CANCELLED'*|*'|TIMEOUT'|*'|OUT_OF_ME'*)
          energy=$(cut -d'|' -f2 <<< "$line")
          elapsed=$(cut -d'|' -f3 <<< "$line")
          break
          ;;
      esac
      sleep 5
    done
  fi
  echo "$label $mode $nodes $np ${status:-?} ${unknowns:--} ${time:--}" \
    "${szz:--} ${steps:-0} ${resolutions:-0} ${newton:-0} ${failures:-0}" \
    "${solver_label:-default} ${mesh_time:--} ${loop_time:--}" \
    "${energy:--} ${elapsed:--} ${cauchy_time:--} ${pp_time:--} ${wall:--}" |
    tee -a "$out/series.txt"
}
for mode in ${modes:-gpu cpu}; do
  if [ "$mode" = gpu ]; then
    np=$((nodes * gpus_per_node))
    threads=$((cores_per_node / gpus_per_node))
    # order of the GPUs of a node (the launcher takes the local rank modulo
    # the length of the list)
    gpu_order=$(python3 -I "$P/scripts/gpus-close-to-ranks.py" \
      "$gpus_per_node" "$threads")
    if [ -n "$gpu_order" ]; then
      if [ -n "${PROTO_GPU_PER_RANK:-}" ]; then
        export PROTO_GPU_RANK_ORDER=$gpu_order
      else
        export PROTO_CUDA_VISIBLE_DEVICES=$gpu_order
      fi
    fi
    device=$gpu
    name=$label-r$refinement-gpu-n$nodes$suffix
    read -r -a gpu_aware <<< "${PROTO_GPU_AWARE_MPI_ENVIRONMENT:-}"
    options="$solver_gpu --gpu-aware-mpi"
    run_directory=$out/$name
    run_timeout=$run_timeout_gpu
    watch_steps "$name"
    run "$name" "$threads" "$np" "$gpu_library" "$batch" "${kernels[@]}" \
      "${gpu_aware[@]}"
    stop_watch "$name"
    unset PROTO_CUDA_VISIBLE_DEVICES PROTO_GPU_RANK_ORDER PROTO_GPU_PER_RANK
    device=""
  else
    np=$((nodes * cores_per_node))
    name=$label-r$refinement-cpu-n$nodes$suffix
    options=$solver_cpu
    run_directory=$out/$name
    run_timeout=$run_timeout_cpu
    watch_steps "$name"
    run "$name" 1 "$np" "$generic" "$batch" "${kernels[@]}" \
      "${cpu_environment[@]}"
    stop_watch "$name"
  fi
  record "$name" "$mode" "$np"
done
echo "results in $out"
