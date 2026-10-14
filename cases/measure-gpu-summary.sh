#!/bin/bash
# Summary on the GPU of the work packages of 11/10/2026 (positions, matrices
# and action without forming B, inverses of the jacobians, global matrix):
# beam r4 on one process and the GPU, everything by kernels (behaviour
# `cuda`), tangent operator assembled or applied without assembly, Newton
# or NewtonK0 with Eisenstat-Walker, with the build of the cases linked to
# the mfem-mgis installed in install/mfem-mgis (before, build/cases) and in
# install/mfem-mgis-global (after, build/cases-global). To be run on a quiet
# machine, with the GPU free.
#
# usage: cases/measure-gpu-summary.sh [results directory]
work_package=gpu-summary
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/gpu-summary-$(date +%Y-%m-%d-%H%M)}
source "$P/cases/measure-common.sh" beam "$out"
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels)
assembled=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
matrix_free=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
device=cuda
linear_tolerance=1e-4
k0="-nls fixed-point -fpm newton-k0 -els amg -k0 -pk0"
ew="-alst EisenstatWalker1 -alsti 0.1"
variants=(
  "assembled-newton:pcg-amg:$assembled:"
  "assembled-newton-k0-ew1:pcg-amg:$assembled:$k0 $ew"
  "matrix-free-newton-k0-ew1:pcg-amg:$matrix_free:$k0 $ew")
for version in before after; do
  if [ "$version" = before ]; then
    B=$P/build/cases/beam
  else
    B=$P/build/cases-global/beam
  fi
  for variant in "${variants[@]}"; do
    IFS=: read -r name linear_solver tangent options <<< "$variant"
    run "beam-np1-gpu-$name-$version" 1 1 "$cuda" "${kernels[@]}" \
      "$tangent"
  done
done
echo "results in $out"
