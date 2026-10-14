#!/bin/bash
# Measurements of the tangent operator applied without assembly (work
# package L5) on the beam, on the GPU (behaviour `cuda`): conjugate gradient
# of MFEM preconditioned by the Jacobi or Chebyshev smoother (orders 2 to 5),
# and NewtonK0 (conjugate gradient preconditioned by BoomerAMG built once on
# K0, with or without Eisenstat-Walker 1), compared to NewtonK0 with the
# assembled tangent operator. To be run on a quiet machine, after
# build-mfem-mgis.sh.
#
# usage: cases/measure-L5-matrix-free.sh [results directory]
work_package=L5
source "$(dirname "$0")/measure-common.sh" beam "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
common=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels)
assembled=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
matrix_free=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
device=cuda
linear_tolerance=1e-4
k0="-nls fixed-point -fpm newton-k0 -els amg -k0 -pk0"
ew="-alst EisenstatWalker1 -alsti 0.1"
# name:linear solver:tangent operator:options
variants=(
  "matrix-free-cg-jacobi:cg-jacobi:$matrix_free:"
  "matrix-free-cg-chebyshev-2:cg-chebyshev:$matrix_free:-cho 2"
  "matrix-free-cg-chebyshev-3:cg-chebyshev:$matrix_free:-cho 3"
  "matrix-free-cg-chebyshev-4:cg-chebyshev:$matrix_free:-cho 4"
  "matrix-free-cg-chebyshev-5:cg-chebyshev:$matrix_free:-cho 5"
  "matrix-free-newton-k0:pcg-amg:$matrix_free:$k0"
  "matrix-free-newton-k0-ew1:pcg-amg:$matrix_free:$k0 $ew"
  "assembled-newton-k0-ew1:pcg-amg:$assembled:$k0 $ew")
for variant in "${variants[@]}"; do
  IFS=: read -r name linear_solver tangent options <<< "$variant"
  run "beam-np1-gpu-$name" 1 1 "$cuda" "$batch" "${common[@]}" "$tangent"
done
echo "results in $out"
