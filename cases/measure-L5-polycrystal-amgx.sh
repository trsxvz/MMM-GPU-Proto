#!/bin/bash
# Measurements of AmgX (work package L5) against BoomerAMG on the
# polycrystal, on the GPU: NewtonK0 (conjugate gradient preconditioned by
# one cycle of the multigrid built once on the elastic operator) with the
# tolerance of Eisenstat-Walker 1, the multigrid being BoomerAMG or AmgX
# (default settings of mfem-mgis, threshold 0.8 and interpolation
# MULTIPASS, in double or in single precision). The Newton computation is
# the reference of the macroscopic stresses. To be run on a quiet machine,
# after build-mfem-mgis.sh.
#
# usage: cases/measure-L5-polycrystal-amgx.sh [results directory]
work_package=L5
source "$(dirname "$0")/measure-common.sh" polycrystal "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
device=cuda
k0="--non-linear-solver NewtonK0"
k0="$k0 --adaptive-linear-solver-tolerance EisenstatWalker1"
k0="$k0 --adaptive-linear-solver-initial-tolerance 0.1"
amgx="--elastic-operator-linear-solver AmgX"
tuned="--amgx-parameters strength_threshold=0.8,interpolator=MULTIPASS"
variants=(
  "newton:"
  "k0-amg:$k0"
  "k0-amgx-default:$k0 $amgx"
  "k0-amgx-tuned:$k0 $amgx $tuned"
  "k0-amgx-tuned-simple:$k0 $amgx $tuned --amgx-single-precision")
for variant in "${variants[@]}"; do
  IFS=: read -r name options <<< "$variant"
  reference=""
  if [ "$name" != newton ]; then
    reference="polycrystal-np1-gpu-newton"
  fi
  run "polycrystal-np1-gpu-$name" 1 1 "$cuda" "$batch" "${kernels[@]}"
done
echo "results in $out"
