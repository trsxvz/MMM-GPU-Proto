#!/bin/bash
# beam r4 of the prototype (111 843 unknowns, 10 steps), all by kernels, one
# process on the GPU (host cores given by CORES, 0-3 by default), Newton:
# PCG of hypre and BoomerAMG (tolerance 1e-4) or cuDSS (symmetric positive
# definite, analysis reused); hypre on the GPU; two passes
P=$(cd "$(dirname "$0")/../.." && pwd)
source $P/scripts/env.sh > /dev/null
out=$P/results/cudss-beam-2026-10-12
L=$P/build/cases/beam/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
cd $P/build/cases-cudss/beam
export OMP_NUM_THREADS=1 MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch \
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels \
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels \
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
for pass in 1 2; do
  for variant in pcg-amg cudss; do
    case $variant in
      pcg-amg) extra="-ls pcg-amg -lt 1e-4" ;;
      cudss) extra="-ls cudss -mra" ;;
    esac
    echo "== pass $pass, $variant, load $(cut -d ' ' -f 1 /proc/loadavg)"
    taskset -c ${CORES:-0-3} timeout 900 mpirun -n 1 ./beam -r 4 -n 10 \
      -mi 3000 -l $L --profiling --device cuda $extra \
      > $out/beam-gpu-$variant-$pass.txt 2>&1
    echo "   exit $?"
  done
done
