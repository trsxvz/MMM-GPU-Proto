#!/bin/bash
# Measurements of the inexact Newton (work package L5) on the beam: fixed
# tolerance of the linear solver (1e-4) or tolerance adapted to the decrease
# of the residual (choices 1 and 2 of Eisenstat and Walker, option -alst),
# for Newton (PCG + BoomerAMG rebuilt at each iteration) and for NewtonK0
# (conjugate gradient preconditioned by BoomerAMG built once on K0,
# prediction by this conjugate gradient). On 8 CPU processes then on the
# GPU (behaviour `cuda`). The fixed tolerance is also used for the
# predictions. To be run on a quiet machine, after build-mfem-mgis.sh.
#
# usage: cases/measure-L5-ew.sh [results directory]
work_package=L5
source "$(dirname "$0")/measure-common.sh" beam "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
k0="-nls fixed-point -fpm newton-k0 -els amg -k0 -pk0"
# name:tolerance of the linear solver:options
variants=(
  "newton:1e-4:"
  "newton-ew1:1e-4:-alst EisenstatWalker1"
  "newton-ew2:1e-4:-alst EisenstatWalker2"
  "newton-k0:1e-4:$k0"
  "newton-k0-ew1:1e-4:$k0 -alst EisenstatWalker1"
  "newton-k0-ew2:1e-4:$k0 -alst EisenstatWalker2")
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
    IFS=: read -r name linear_tolerance options <<< "$variant"
    run "beam-$where-$name" 1 $np "$library" "$batch" "${kernels[@]}"
  done
done
echo "results in $out"
