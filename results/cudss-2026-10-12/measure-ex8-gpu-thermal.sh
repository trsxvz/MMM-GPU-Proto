#!/bin/bash
# ex8 on the GPU, thermal and mechanical problems solved by cuDSS with the
# reuse of the analysis (thermal matrices given whole, mechanical ones by
# their lower triangle), hypre on the host or on the GPU; one process,
# host on the cores given by CORES (0-3 by default); two passes; timings in
# summary-thermal.txt
P=$(cd "$(dirname "$0")/../.." && pwd)
source $P/scripts/env.sh > /dev/null
out=$P/results/cudss-2026-10-12
L=$P/build/cas-memoire/ex8/generic-parallel-cuda/src
cd $P/build/examples-mumps/ex8
export OMP_NUM_THREADS=1 MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch \
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels \
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels \
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
for pass in 1 2; do
  for variant in hypre-hote hypre-gpu; do
    case $variant in
      hypre-hote) extra="--hypre-on-host" ;;
      hypre-gpu) extra="" ;;
    esac
    echo "== pass $pass, $variant, load $(cut -d ' ' -f 1 /proc/loadavg)"
    taskset -c ${CORES:-0-3} timeout 900 mpirun -n 1 ./rjh_plate \
      -m assemblage_hexa.msh -r 0 -ns 5 -msMc SymmetricPositiveDefinite \
      -lsMc CuDSSSolver -lsTh CuDSSSolver -mra -no-pp --device cuda \
      -lU $L/libGenericParallelU3SI2-cuda.so \
      -lA $L/libGenericParallelALFENI-cuda.so $extra \
      > $out/ex8-gpu-thermique-cudss-$variant-$pass.txt 2>&1
    echo "   exit $?"
  done
done
