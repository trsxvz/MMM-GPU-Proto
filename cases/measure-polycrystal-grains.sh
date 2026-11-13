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
# - run_timeout: time limit of each computation in s (default 7200);
# - solver_gpu, solver_cpu: options of the solvers of the polycrystal on the
#   GPU and on the CPU (default: NewtonK0, Eisenstat-Walker 1 from 0.1, the
#   best on both sides on Jean Zay, results/L5-polycrystal-2026-10-10-0502);
# - solver_label: suffix of the names of the computations (comparison of
#   solvers in the same results directory).
# The GPU computation runs one process per GPU (cores_per_node /
# gpus_per_node cores each) with GPU-aware MPI, the CPU computation one
# process per core. The fixed point on the lateral stresses uses the
# secant acceleration, the averages of the stresses per grain are written
# at each step (post-processings without ParaView export) in
# <results>/<computation>/avgStressPolycrystal. A line per computation is
# appended to <results>/series.txt: label, mode, nodes, processes, status,
# unknowns, time of the steps (s), axial stress at the end (Pa), time steps
# (with the divided ones), resolutions of the fixed point on the lateral
# stresses, Newton iterations, failed resolutions, solver label.
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
run_timeout=${run_timeout:-7200}
post_processing=""
source "$(dirname "$0")/measure-common.sh" polycrystal "$@"
summary="$summary|number of finite element unknowns"
executable=$B/uniaxial-polycrystal
case_options="-m $mesh -f $vectors --order $order -r $refinement"
case_options="$case_options --macroscopic-secant-acceleration"
case_options="$case_options --disable-paraview-export"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
gpu=${PROTO_GPU_DEVICE:-cuda}
gpu_library=$gp/generic-parallel-$gpu/src/libGenericParallelBehaviour-$gpu.so
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
# usage: record name mode processes
record() {
  local name=$1 mode=$2 np=$3 status unknowns time szz steps resolutions
  local newton failures
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
  echo "$label $mode $nodes $np ${status:-?} ${unknowns:--} ${time:--}" \
    "${szz:--} ${steps:-0} ${resolutions:-0} ${newton:-0} ${failures:-0}" \
    "${solver_label:-default}" | tee -a "$out/series.txt"
}
for mode in ${modes:-gpu cpu}; do
  if [ "$mode" = gpu ]; then
    np=$((nodes * gpus_per_node))
    threads=$((cores_per_node / gpus_per_node))
    # order of the GPUs of a node (the launcher takes the local rank modulo
    # the length of the list)
    order=$(python3 -I "$P/scripts/gpus-close-to-ranks.py" "$gpus_per_node" \
      "$threads")
    if [ -n "$order" ]; then
      if [ -n "${PROTO_GPU_PER_RANK:-}" ]; then
        export PROTO_GPU_RANK_ORDER=$order
      else
        export PROTO_CUDA_VISIBLE_DEVICES=$order
      fi
    fi
    device=$gpu
    name=$label-r$refinement-gpu-n$nodes$suffix
    read -r -a gpu_aware <<< "${PROTO_GPU_AWARE_MPI_ENVIRONMENT:-}"
    options="$solver_gpu --gpu-aware-mpi"
    run_directory=$out/$name
    run "$name" "$threads" "$np" "$gpu_library" "$batch" "${kernels[@]}" \
      "${gpu_aware[@]}"
    unset PROTO_CUDA_VISIBLE_DEVICES PROTO_GPU_RANK_ORDER PROTO_GPU_PER_RANK
    device=""
  else
    np=$((nodes * cores_per_node))
    name=$label-r$refinement-cpu-n$nodes$suffix
    options=$solver_cpu
    run_directory=$out/$name
    run "$name" 1 "$np" "$generic" "$batch" "${kernels[@]}"
  fi
  record "$name" "$mode" "$np"
done
echo "results in $out"
