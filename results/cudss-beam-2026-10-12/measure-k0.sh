#!/bin/bash
# beam r4 of the prototype (111 843 unknowns, 10 steps), all by kernels, one
# process on the GPU (host cores given by CORES, 0-3 by default), NewtonK0
# (fixed point solver, elastic operator K0 assembled once, prediction by
# K0): K0 preconditioned by one cycle of BoomerAMG (best variant measured
# on 10/10, with the tolerance of Eisenstat and Walker) or factorized by
# cuDSS, fixed tolerance (1e-4) or tolerance of Eisenstat and Walker; two
# passes
P=$(cd "$(dirname "$0")/../.." && pwd)
source $P/scripts/env.sh > /dev/null
out=$P/results/cudss-beam-2026-10-12
L=$P/build/cases/beam/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
cd $P/build/cases-cudss/beam
export OMP_NUM_THREADS=1 MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch \
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels \
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels \
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
k0="-nls fixed-point -fpm newton-k0 -k0 -pk0 -ls pcg-amg -lt 1e-4"
for pass in 1 2; do
  for variant in k0-amg-ew1 k0-cudss k0-cudss-ew1; do
    case $variant in
      k0-amg-ew1) extra="$k0 -els amg -alst EisenstatWalker1" ;;
      k0-cudss) extra="$k0 -els cudss" ;;
      k0-cudss-ew1) extra="$k0 -els cudss -alst EisenstatWalker1" ;;
    esac
    echo "== pass $pass, $variant, load $(cut -d ' ' -f 1 /proc/loadavg)"
    taskset -c ${CORES:-0-3} timeout 900 mpirun -n 1 ./beam -r 4 -n 10 \
      -mi 3000 -l $L --profiling --device cuda $extra \
      > $out/beam-gpu-$variant-$pass.txt 2>&1
    echo "   exit $?"
  done
done
