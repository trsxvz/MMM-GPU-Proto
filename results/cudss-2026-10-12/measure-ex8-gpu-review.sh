#!/bin/bash
# ex8 on the GPU after the independent review of the cuDSS solver (solutions
# checked by their residual, pattern compared without the values): cuDSS
# with the reuse of the analysis, hypre on the host or on the GPU, and
# without the reuse, hypre on the GPU; one process, host on the cores given
# by CORES (0-3 by default); two passes; timings in summary-review.txt
P=$(cd "$(dirname "$0")/../.." && pwd)
source $P/scripts/env.sh > /dev/null
out=$P/results/cudss-2026-10-12
L=$P/build/cases-memory/ex8/generic-parallel-cuda/src
cd $P/build/examples-mumps/ex8
export OMP_NUM_THREADS=1 MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch \
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels \
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels \
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
for pass in 1 2; do
  for variant in cudss-with-hypre-host cudss-with-hypre-gpu cudss-without-hypre-gpu; do
    case $variant in
      cudss-with-hypre-host) extra="-mra --hypre-on-host" ;;
      cudss-with-hypre-gpu) extra="-mra" ;;
      cudss-without-hypre-gpu) extra="" ;;
    esac
    echo "== pass $pass, $variant, load $(cut -d ' ' -f 1 /proc/loadavg)"
    taskset -c ${CORES:-0-3} timeout 900 mpirun -n 1 ./rjh_plate \
      -m assemblage_hexa.msh -r 0 -ns 5 -msMc SymmetricPositiveDefinite \
      -lsMc CuDSSSolver -no-pp --device cuda \
      -lU $L/libGenericParallelU3SI2-cuda.so \
      -lA $L/libGenericParallelALFENI-cuda.so $extra \
      > $out/ex8-gpu-review-$variant-$pass.txt 2>&1
    echo "   exit $?"
  done
done
