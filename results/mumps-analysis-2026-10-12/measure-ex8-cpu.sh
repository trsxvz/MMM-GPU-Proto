#!/bin/bash
# ex8 on the CPU (2 processes of one thread, cores 0-3, all by kernels,
# generic laws), with and without the reuse of the analysis of MUMPS
# (option -mra), two interleaved passes; timings in resume-*.txt
P=$(cd "$(dirname "$0")/../.." && pwd)
source $P/scripts/env.sh > /dev/null
out=$P/results/mumps-analysis-2026-10-12
L=$P/build/mfem-mgis-examples/ex8/src
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
    taskset -c 0-3 timeout 900 mpirun -n 2 ./rjh_plate -m assemblage_hexa.msh \
      -r 0 -ns 5 -lsMc MUMPSSolver -msMc SymmetricPositiveDefinite -no-pp \
      -lU $L/libU3SI2-generic.so -lA $L/libALFENI-generic.so $extra \
      > $out/ex8-np2-noyaux-$variant-$pass.txt 2>&1
    echo "   exit $?"
  done
done
