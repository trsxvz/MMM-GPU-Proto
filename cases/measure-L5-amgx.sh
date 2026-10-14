#!/bin/bash
# Measurements of AmgX (work package L5) against BoomerAMG on the beam:
# 1. elastic operator K0 alone (beam-k0): setup and resolutions of the
#    conjugate gradient preconditioned by each multigrid, on the GPU;
# 2. complete computation by the best method, NewtonK0 (conjugate gradient
#    preconditioned by one cycle of the multigrid built once on K0,
#    prediction by this conjugate gradient) with the tolerance of
#    Eisenstat-Walker 1, the multigrid of K0 being BoomerAMG or AmgX, on the
#    GPU then on 8 CPU processes (AmgX then computes on the GPU, shared by
#    the 8 processes).
# Parameters of AmgX given by -amgxo (names of AmgX): default settings of
# mfem-mgis (the ones of MFEM: strength threshold 0.25, interpolation D2),
# threshold 0.8 and interpolation MULTIPASS (best setting of the sweep on
# K0), with two sweeps of the smoother, and in single precision (-amgxs).
# On 8 processes, each process has its AmgX resources on the same GPU.
# To be run on a quiet machine, after build-mfem-mgis.sh.
#
# usage: cases/measure-L5-amgx.sh [results directory]
work_package=L5
source "$(dirname "$0")/measure-common.sh" beam "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
amgx_tuned=strength_threshold=0.8,interpolator=MULTIPASS
# 1. operator K0 alone, on the GPU (name|options of beam-k0)
k0_variants=(
  "amg|-s cg-amg"
  "amg-aggressive|-s cg-amg -amgo AggressiveCoarseningLevels=1"
  "pcg-amg|-s pcg-amg"
  "amgx-default|-s cg-amgx"
  "amgx-tuned|-s cg-amgx -amgxo $amgx_tuned"
  "amgx-tuned-sweeps2|-s cg-amgx -amgxo $amgx_tuned,presweeps=2,postsweeps=2"
  "amgx-tuned-simple|-s cg-amgx -amgxs -amgxo $amgx_tuned"
  "amgx-default-simple|-s cg-amgx -amgxs"
  "amgx-pcg-tuned|-s amgx-pcg -amgxo $amgx_tuned"
  "jacobi|-s cg-jacobi")
for variant in "${k0_variants[@]}"; do
  IFS='|' read -r name k0_options <<< "$variant"
  echo "== k0-$name" | tee -a "$out/summary.txt"
  read -r -a k0_cmd <<< "./beam-k0 -r 4 -n 5 --device cuda $k0_options"
  (cd "$B" && env OMP_NUM_THREADS=1 timeout 1800 "${k0_cmd[@]}") \
    < /dev/null > "$out/k0-$name.txt" 2>&1
  echo "   status $?" | tee -a "$out/summary.txt"
  /usr/bin/grep -E "^dofs|resolution" "$out/k0-$name.txt" |
    tee -a "$out/summary.txt"
done
# 2. complete computation
linear_solver=pcg-amg
linear_tolerance=1e-4
base="-nls fixed-point -fpm newton-k0 -k0 -pk0 -alst EisenstatWalker1 -alsti 0.1"
# name:options of the beam
variants=(
  "amg:-els amg"
  "amg-aggressive:-els amg -amgo AggressiveCoarseningLevels=1"
  "amgx-default:-els amgx"
  "amgx-tuned:-els amgx -amgxo $amgx_tuned"
  "amgx-tuned-sweeps2:-els amgx -amgxo $amgx_tuned,presweeps=2,postsweeps=2"
  "amgx-tuned-simple:-els amgx -amgxs -amgxo $amgx_tuned")
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
    IFS=: read -r name extra <<< "$variant"
    options="$base $extra"
    run "beam-$where-$name" 1 $np "$library" "$batch" "${kernels[@]}"
  done
done
echo "results in $out"
