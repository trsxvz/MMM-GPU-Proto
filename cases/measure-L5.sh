#!/bin/bash
# Measurements of the fixed point (work package L5) on the beam: best Newton
# (tangent operator assembled by kernels, PCG + BoomerAMG of hypre,
# tolerance 1e-4) against the fixed point on the elastic operator K0,
# accelerated (Anderson of depth 8, Cast3M), K0 solved by one cycle of
# BoomerAMG kept for the whole computation (default smoother or Chebyshev
# smoother), or by PCG + BoomerAMG kept, and against the accelerated
# modified Newton. For the fixed point, the prediction is solved up to 1e-8
# (otherwise the elastic steps require hundreds of iterations). On 8 CPU
# processes then on the GPU (behaviour `cuda`). To be run on a quiet
# machine, after build-mfem-mgis.sh.
#
# usage: cases/measure-L5.sh [results directory]
work_package=L5
source "$(dirname "$0")/measure-common.sh" beam "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
# name:tolerance of the linear solver:options
fp="-nls fixed-point"
variants=(
  "newton:1e-4:"
  "fixed-point-anderson8-amg-k0:1e-8:$fp -acc anderson -ad 8 -els amg -k0"
  "fixed-point-anderson8-amg-chebyshev-k0:1e-8:$fp -acc anderson -ad 8 -els amg -k0 -amgo RelaxationType=16"
  "fixed-point-cast3m-amg-k0:1e-8:$fp -acc cast3m -els amg -k0"
  "fixed-point-anderson8-pcg-amg-k0:1e-8:$fp -acc anderson -ad 8 -els pcg-amg -k0"
  "newton-modified-anderson8:1e-8:$fp -fpm modified-newton -acc anderson -ad 8")
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
