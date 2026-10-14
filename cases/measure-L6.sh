#!/bin/bash
# Measurements of the GPU slice (work packages L3 to L6): reference
# computation on the CPU (8 processes, point by point and everything by
# kernels), then on the GPU with all the kernels, the behaviour `generic`
# (integrated on the host) or the behaviour `cuda` (integrated on the GPU),
# on 1 and 2 processes. To be run on a quiet machine, after
# build-mfem-mgis.sh.
#
# usage: cases/measure-L6.sh [beam|polycrystal] [results directory]
# (see cases/measure-common.sh)
work_package=L6
source "$(dirname "$0")/measure-common.sh" "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
device=""
# reference of the comparisons of the polycrystal: computation point by
# point on the same number of processes (see measure-common.sh)
run "$case_name-np8-point-by-point" 1 8 "$generic"
run "$case_name-np1-point-by-point" 1 1 "$generic"
run "$case_name-np2-point-by-point" 1 2 "$generic"
run "$case_name-np8-kernels" 1 8 "$generic" "$batch" "${kernels[@]}"
device=cuda
run "$case_name-np1-gpu-kernels-behaviour-generic" 1 1 "$generic" "$batch" \
  "${kernels[@]}"
run "$case_name-np1-gpu-kernels-behaviour-cuda" 1 1 "$cuda" "$batch" "${kernels[@]}"
run "$case_name-np2-gpu-kernels-behaviour-cuda" 1 2 "$cuda" "$batch" "${kernels[@]}"
echo "results in $out"
