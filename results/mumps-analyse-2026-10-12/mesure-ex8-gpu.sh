#!/bin/bash
# ex8 on the GPU (one process, host on cores 4-7, all by kernels, cuda
# laws, hypre on the host), with and without the reuse of the analysis of
# MUMPS (option -mra), two interleaved passes
P=$(cd "$(dirname "$0")/../.." && pwd)
source $P/scripts/env.sh > /dev/null
out=$P/results/mumps-analyse-2026-10-12
L=$P/build/cas-memoire/ex8/generic-parallel-cuda/src
cd $P/build/examples-mumps/ex8
export OMP_NUM_THREADS=1 MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch \
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels \
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels \
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
for pass in 1 2; do
  for variant in sans avec; do
    extra=""
    [ $variant = avec ] && extra=-mra
    echo "== pass $pass, $variant, load $(cut -d ' ' -f 1 /proc/loadavg)"
    taskset -c 4-7 timeout 900 mpirun -n 1 ./rjh_plate -m assemblage_hexa.msh \
      -r 0 -ns 5 -lsMc MUMPSSolver -msMc SymmetricPositiveDefinite -no-pp \
      --device cuda --hypre-on-host \
      -lU $L/libGenericParallelU3SI2-cuda.so \
      -lA $L/libGenericParallelALFENI-cuda.so $extra \
      > $out/ex8-gpu-noyaux-$variant-$pass.txt 2>&1
    echo "   exit $?"
  done
done
