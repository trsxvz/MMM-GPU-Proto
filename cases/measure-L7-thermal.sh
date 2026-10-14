#!/bin/bash
# Measurements of the stationary non linear heat transfer by kernels (work
# package L7) in the beam (cases/thermal, r5, 280,000 unknowns): on 8 CPU
# processes, point by point then by kernels; on the GPU (behaviour `cuda`),
# by kernels, tangent operator assembled (GMRES of hypre and BoomerAMG,
# fixed tolerance or Eisenstat-Walker 1) or applied without assembly (GMRES
# of MFEM and Jacobi). The thermal tangent operator is negative definite and
# the imposed degrees of freedom have a unit diagonal: GMRES instead of the
# conjugate gradient. To be run on a quiet machine, after
# build-mfem-mgis.sh.
#
# usage: cases/measure-L7-thermal.sh [results directory]
work_package=L7
source "$(dirname "$0")/measure-common.sh" thermal "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
common=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels)
assembled=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
matrix_free=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_tolerance=1e-8
linear_solver=gmres-amg
run thermal-np8-point-by-point 1 8 "$generic"
run thermal-np8-kernels 1 8 "$generic" "$batch" "${common[@]}" "$assembled"
device=cuda
run thermal-np1-gpu-kernels 1 1 "$cuda" "$batch" "${common[@]}" "$assembled"
options="-alst EisenstatWalker1"
run thermal-np1-gpu-kernels-ew1 1 1 "$cuda" "$batch" "${common[@]}" \
  "$assembled"
options=""
linear_solver=gmres-jacobi
run thermal-np1-gpu-matrix-free 1 1 "$cuda" "$batch" "${common[@]}" \
  "$matrix_free"
echo "results in $out"
