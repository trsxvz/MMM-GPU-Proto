#!/bin/bash
# Measurements of the work package L3 (computations by kernels): integration
# by batches with the gradients, the inner forces and the element matrices
# of the tangent operator computed element by element or by kernels, on 1
# and 8 processes, on the CPU. To be run on a quiet machine, after
# build-mfem-mgis.sh.
#
# usage: cases/measure-L3.sh [beam|polycrystal] [results directory]
# (see cases/measure-common.sh)
work_package=L3
source "$(dirname "$0")/measure-common.sh" "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
gradients=MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
forces=MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
tangent=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
for np in 8 1; do
  run "$case_name-np$np-point-by-point" 1 "$np" "$generic"
  run "$case_name-np$np-batches" 1 "$np" "$generic" "$batch"
  run "$case_name-np$np-batches-kernels-gradients" 1 "$np" "$generic" \
    "$batch" "$gradients"
  run "$case_name-np$np-batches-kernels-forces" 1 "$np" "$generic" \
    "$batch" "$forces"
  run "$case_name-np$np-batches-kernels-tangent" 1 "$np" "$generic" \
    "$batch" "$tangent"
  run "$case_name-np$np-batches-kernels" 1 "$np" "$generic" \
    "$batch" "$gradients" "$forces" "$tangent"
done
echo "results in $out"
