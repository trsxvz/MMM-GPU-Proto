#!/bin/bash
# Scaling (measurement campaign, case 5): beam refined three, four then
# five times (16,000 to 840,000 unknowns, 10 steps), best method of the
# work package L5 (NewtonK0, conjugate gradient preconditioned by the
# multigrid of K0, prediction, tolerance of Eisenstat-Walker 1), all the
# kernels:
# - on the GPU, aggressive BoomerAMG multigrid or AmgX tuned in single
#   precision, assembled tangent operator, and aggressive BoomerAMG with the
#   tangent operator applied without assembly (only K0 is assembled), with
#   the peak memory used on the GPU (sampled by nvidia-smi every 200 ms,
#   other processes included);
# - on 8 CPU processes, aggressive BoomerAMG.
# On the GPU, the memory pool of hypre (Umpire, 4 GiB by default, allocated
# at once) is reduced to 0.5 GiB (option -hdps of mfem-mgis): times within
# a few %, 3.5 GiB less. The next refinement (6.5 million unknowns) does not
# fit in the 16 GB of the GPU, even so (allocation refused around 15.3 GiB,
# with the tangent operator assembled or not).
# Before each computation on the GPU, the script waits at most 10 minutes
# for another process to free the GPU, then reports it in the summary.
# To be run on a quiet machine, after build-mfem-mgis.sh.
#
# usage: cases/measure-beam-scaling.sh [results directory]
work_package=scaling
source "$(dirname "$0")/measure-common.sh" beam "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
matrix_free=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
linear_tolerance=1e-4
base="-nls fixed-point -fpm newton-k0 -k0 -pk0 -alst EisenstatWalker1 -alsti 0.1"
pool="-hdps 0.5"
amg="-els amg -amgo AggressiveCoarseningLevels=1"
amgx="-els amgx -amgxs -amgxo strength_threshold=0.8,interpolator=MULTIPASS"

# waits at most 10 minutes for the GPU to be free of other processes
wait_for_gpu() {
  local i
  for ((i = 0; i != 120; ++i)); do
    if [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; then
      return 0
    fi
    sleep 5
  done
  echo "   GPU shared with another process" | tee -a "$out/summary.txt"
}

# run on the GPU with the peak of the memory used on the GPU
# usage: run_gpu name settings...
run_gpu() {
  local name=$1
  shift
  wait_for_gpu
  nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -lms 200 \
    > "$out/$name.memory.txt" &
  local sampler=$!
  run "$name" 1 1 "$cuda" "$@"
  kill "$sampler"
  wait "$sampler" 2> /dev/null
  echo "   memory of the GPU at most $(sort -n "$out/$name.memory.txt" |
    tail -1) MiB" | tee -a "$out/summary.txt"
}

for refinement in 3 4 5; do
  device=cuda
  options="$base $pool $amg"
  run_gpu "beam-r$refinement-gpu-amg-aggressive" "$batch" "${kernels[@]}"
  options="$base $pool $amgx"
  run_gpu "beam-r$refinement-gpu-amgx-simple" "$batch" "${kernels[@]}"
  options="$base $pool $amg"
  run_gpu "beam-r$refinement-gpu-amg-aggressive-matrix-free" "$batch" \
    "${matrix_free[@]}"
  device=""
  options="$base $amg"
  run "beam-r$refinement-np8-amg-aggressive" 1 8 "$generic" "$batch" \
    "${kernels[@]}"
done
echo "results in $out"
