#!/bin/bash
# Ordering of the unknowns on the GPU (parameter FiniteElementSpaceOrdering
# of mfem-mgis): beam r4 (111,843 unknowns), 10 steps, Newton, on the GPU
# (`--execution gpu`, one process, behaviour integrated by batches,
# generic-parallel interface, cuda backend, finite elements computed by
# kernels), the unknowns being ordered by nodes (byNODES, default of
# mfem-mgis) or by vector dimension (byVDIM):
#
# - pcg-amg: PCG of hypre and BoomerAMG (on the GPU, hypre has no
#   elasticity options: mfem-mgis uses the systems options in both
#   orderings);
# - cg-amgx: conjugate gradient of MFEM and one V-cycle of the classical
#   multigrid of AmgX;
# - cg-amgx-aggregation: the same with the aggregation multigrid (-amgxa),
#   without and, in byVDIM only, with nodal blocks of 3 x 3 values
#   (-amgxnb).
#
# The time of each phase (cases/summary-profiling.py <directory>/*.txt),
# in particular of the setup of the preconditioner and of the resolutions,
# decided whether the unknowns are ordered by vector dimension by default
# when the finite elements are computed on the GPU: they are not (JOURNAL.md,
# 10/11/2026).
#
# usage: cases/measure-ordering-gpu.sh [results directory]
work_package=ordering-gpu
source "$(dirname "$0")/measure-common.sh" beam "${1:-}"
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_tolerance=1e-4
summary="^dofs|^time|iterations"
settings=(MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered)
for ordering in byNODES byVDIM; do
  for solver in pcg-amg cg-amgx; do
    linear_solver=$solver
    options="--execution gpu -ord $ordering"
    run "ordering-$solver-$ordering" 1 1 "$cuda" "${settings[@]}"
  done
  linear_solver=cg-amgx
  options="--execution gpu -ord $ordering -amgxa"
  run "ordering-cg-amgx-aggregation-$ordering" 1 1 "$cuda" "${settings[@]}"
done
linear_solver=cg-amgx
options="--execution gpu -ord byVDIM -amgxa -amgxnb"
run "ordering-cg-amgx-aggregation-nodal-blocks-byVDIM" 1 1 "$cuda" \
  "${settings[@]}"
echo "results in $out"
