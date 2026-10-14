#!/bin/bash
# End to end measurements on the GPU of the inverses of the jacobians
# computed once (work package 7j): beam r4 and thermal case r5 on one
# process and the GPU, everything by kernels (behaviour `cuda`), with the
# build of the cases linked to the mfem-mgis installed in
# install/mfem-mgis-registers (before, build/cases-registers) and in
# install/mfem-mgis-jacobians (after, build/cases-jacobians). The
# libraries of the behaviours are the ones of build/cases. To be run on a
# quiet machine, with the GPU free.
#
# usage: cases/measure-jacobians-gpu.sh [results directory]
work_package=jacobians
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/jacobians-gpu-$(date +%Y-%m-%d-%H%M)}
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
for case_name in beam thermal; do
  source "$P/cases/measure-common.sh" "$case_name" "$out"
  if [ "$case_name" = thermal ]; then
    linear_solver=gmres-amg
    linear_tolerance=1e-8
  fi
  cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
  device=cuda
  # measured versions, each one with its build directory
  # build/cases-<version>
  for version in registers jacobians; do
    B=$P/build/cases-$version/$case_name
    run "$case_name-np1-gpu-kernels-$version" 1 1 "$cuda" "${kernels[@]}"
  done
done
echo "results in $out"
