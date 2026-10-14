#!/bin/bash
# Measurements of the matrix free tangent operator (work package L4) on the
# beam: tangent operator assembled by kernels and PCG + BoomerAMG of hypre,
# against the tangent operator applied without assembly by kernels
# (MatrixFree) and the conjugate gradient of MFEM preconditioned by
# BoomerAMG built on the assembled tangent operator, rebuilt every 1, 2 or
# 4 Newton iterations (solver cg-amg, option -pu). On 8 CPU processes then
# on the GPU (behaviour `cuda`). To be run on a quiet machine, after
# build-mfem-mgis.sh.
#
# usage: cases/measure-L4.sh [results directory]
work_package=L4
source "$(dirname "$0")/measure-common.sh" beam "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels)
assembled=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
matrix_free=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
for device in "" cuda; do
  if [ -z "$device" ]; then
    np=8
    library=$generic
    where=np8
  else
    np=1
    library=$cuda
    where=np1-gpu
  fi
  linear_solver=pcg-amg
  options=""
  run "beam-$where-assembled-pcg-amg" 1 $np "$library" "$batch" \
    "${kernels[@]}" "$assembled"
  linear_solver=cg-amg
  for period in 1 2 4; do
    options="-pu $period"
    run "beam-$where-matrix-free-cg-amg-pu$period" 1 $np "$library" \
      "$batch" "${kernels[@]}" "$matrix_free"
  done
done
echo "results in $out"
