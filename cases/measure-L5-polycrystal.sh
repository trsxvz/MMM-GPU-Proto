#!/bin/bash
# Measurements of NewtonK0 and of the inexact Newton (work package L5) on the
# polycrystal of mm-opera-hpc: Newton (PCG + BoomerAMG of hypre, new setup
# of BoomerAMG at each iteration) against NewtonK0 (method `NewtonK0` of
# the solver `FixedPoint` of mfem-mgis: conjugate gradient preconditioned
# by BoomerAMG built once on the elastic operator, obtained by an
# integration of the behaviour), with the tolerance of the linear solver of
# the polycrystal (1e-12) or a tolerance adapted to the decrease of the
# residual (Eisenstat-Walker 1, initial tolerance 0.1). All the kernels, on
# 8 CPU processes (behaviour `generic`) then on the GPU (behaviour `cuda`).
# The macroscopic stresses are compared to the ones of Newton in the same
# configuration. To be run on a quiet machine, after build-mfem-mgis.sh.
#
# usage: cases/measure-L5-polycrystal.sh [results directory]
work_package=L5
source "$(dirname "$0")/measure-common.sh" polycrystal "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
ew="--adaptive-linear-solver-tolerance EisenstatWalker1"
ew="$ew --adaptive-linear-solver-initial-tolerance 0.1"
# name:options
variants=(
  "newton:"
  "newton-ew1:$ew"
  "newton-k0:--non-linear-solver NewtonK0"
  "newton-k0-ew1:--non-linear-solver NewtonK0 $ew")
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
  for variant in "${variants[@]}"; do
    IFS=: read -r name options <<< "$variant"
    reference=""
    if [ "$name" != newton ]; then
      reference="polycrystal-$where-newton"
    fi
    run "polycrystal-$where-$name" 1 $np "$library" "$batch" "${kernels[@]}"
  done
done
echo "results in $out"
