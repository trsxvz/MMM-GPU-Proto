#!/bin/bash
# Measurements of the options of BoomerAMG (work package L5) on the beam,
# for the best method: NewtonK0 (conjugate gradient preconditioned by
# BoomerAMG built once on K0, prediction by this conjugate gradient) with
# the tolerance of Eisenstat-Walker 1 (initial tolerance 0.1). Options of
# hypre given by -amgo (strength threshold, aggressive coarsening,
# smoother, interpolation, number of sweeps). On the GPU (behaviour `cuda`)
# then on 8 CPU processes. To be run on a quiet machine, after
# build-mfem-mgis.sh.
#
# usage: cases/measure-L5-amg.sh [results directory]
work_package=L5
source "$(dirname "$0")/measure-common.sh" beam "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
linear_tolerance=1e-4
base="-nls fixed-point -fpm newton-k0 -els amg -k0 -pk0 -alst EisenstatWalker1 -alsti 0.1"
# name:options of BoomerAMG (empty: default options of mfem-mgis)
variants=(
  "default:"
  "threshold-05:StrengthThreshold=0.5"
  "threshold-07:StrengthThreshold=0.7"
  "aggressive-1:AggressiveCoarseningLevels=1"
  "hmis:CoarseningType=10"
  "jacobi:RelaxationType=7"
  "chebyshev:RelaxationType=16"
  "sweeps-2:RelaxationSweeps=2"
  "interp-ext:InterpolationType=14"
  "threshold-05-aggressive-1:StrengthThreshold=0.5,AggressiveCoarseningLevels=1")
for device in cuda ""; do
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
    IFS=: read -r name amg <<< "$variant"
    options="$base${amg:+ -amgo $amg}"
    run "beam-$where-amg-$name" 1 $np "$library" "$batch" "${kernels[@]}"
  done
done
echo "results in $out"
