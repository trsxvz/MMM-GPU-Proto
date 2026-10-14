#!/bin/bash
# beam r4 of the prototype (111 843 unknowns, 10 steps), all by kernels, one
# process on the GPU (host cores given by CORES, 0-3 by default), hypre on
# the GPU: stiffness matrices added by atomic additions (Atomic) or row by
# row in a fixed order (Ordered), with Newton (PCG and BoomerAMG, tangent
# assembled at each iteration) and with NewtonK0 (K0 assembled once, one
# cycle of BoomerAMG, tolerance of Eisenstat and Walker); two passes
P=$(cd "$(dirname "$0")/../.." && pwd)
source $P/scripts/env.sh > /dev/null
out=$P/results/ordered-assembly-2026-10-12
L=$P/build/cas/poutre/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
cd $P/build/cas-ordered/poutre
export OMP_NUM_THREADS=1 MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch \
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels \
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels \
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
k0="-nls fixed-point -fpm newton-k0 -k0 -pk0 -ls pcg-amg -lt 1e-4 -els amg -alst EisenstatWalker1"
for pass in 1 2; do
  for variant in newton k0; do
    for assembly in Atomic Ordered; do
      case $variant in
        newton) extra="-ls pcg-amg -lt 1e-4" ;;
        k0) extra="$k0" ;;
      esac
      echo "== pass $pass, $variant, $assembly, load $(cut -d ' ' -f 1 /proc/loadavg)"
      MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=$assembly \
        taskset -c ${CORES:-0-3} timeout 900 mpirun -n 1 ./poutre -r 4 \
        -n 10 -mi 3000 -l $L --profiling --device cuda $extra \
        > $out/poutre-gpu-$variant-$assembly-$pass.txt 2>&1
      echo "   exit $?"
    done
  done
done
