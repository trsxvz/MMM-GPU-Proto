#!/bin/bash
# Measurements of the notched bar ssna303 in 3D (measurement campaign, case
# 2): example ex4 of mfem-mgis-examples, variant solved by hypre (Newton,
# FGMRES of hypre preconditioned by BoomerAMG), finite strain plasticity
# (logarithmic strain of Hencky), complete loading in 10 steps, mesh refined
# zero, one and two times (28,000 to 1.7 million unknowns):
# - point by point then by batches and kernels on 8 CPU processes;
# - by kernels on the GPU (behaviour `cuda`), with the memory pool of hypre
#   reduced to 0.5 GiB;
# then the same computations with the FBar formulation (`--use-fbar`), by
# kernels since the work package L7.9.
# Before each computation on the GPU, the script waits at most 10 minutes
# for another process to free the GPU, then reports it in the summary.
# At refinement 2 (1.46 million unknowns), the integration rule with 27
# points of the hexahedra does not fit in the memory of the GPU (16 GB): the
# computations on the GPU, and one computation by kernels on 8 processes for
# comparison, use the rule of order 3 (8 points, option -iro 3); FBar on
# 8 processes exhausts the memory of the machine (30 GB) and is not run.
# The variable SSNA303_REFINEMENTS gives the refinements (0 1 2 by default).
# To be run on a quiet machine, after build-mfem-mgis.sh (examples and cases).
#
# usage: cases/measure-ssna303.sh [results directory]
work_package=ssna303
source "$(dirname "$0")/measure-common.sh" ssna303 "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so

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

for refinement in ${SSNA303_REFINEMENTS:-0 1 2}; do
  for fbar in "" "--use-fbar"; do
    suffix=${fbar:+-fbar}
    if [ "$refinement" != 2 ]; then
      device=""
      options="$fbar"
      run "ssna303-r$refinement$suffix-np8-point-by-point" 1 8 "$generic"
      run "ssna303-r$refinement$suffix-np8-kernels" 1 8 "$generic" "$batch" \
        "${kernels[@]}"
      device=cuda
      options="$fbar -hdps 0.5"
      wait_for_gpu
      run "ssna303-r$refinement$suffix-gpu-kernels" 1 1 "$cuda" "$batch" \
        "${kernels[@]}"
      continue
    fi
    # refinement 2: rule of order 3 on the GPU and for the comparison, no
    # FBar on 8 processes
    device=""
    if [ -z "$fbar" ]; then
      options=""
      run "ssna303-r2-np8-point-by-point" 1 8 "$generic"
      run "ssna303-r2-np8-kernels" 1 8 "$generic" "$batch" "${kernels[@]}"
      options="-iro 3"
      run "ssna303-r2-iro3-np8-kernels" 1 8 "$generic" "$batch" \
        "${kernels[@]}"
    fi
    device=cuda
    options="$fbar -iro 3 -hdps 0.5"
    wait_for_gpu
    run "ssna303-r2$suffix-iro3-gpu-kernels" 1 1 "$cuda" "$batch" \
      "${kernels[@]}"
  done
done
echo "results in $out"
