#!/bin/bash
# ex7 complete (40 steps, 2634 unknowns), all by kernels: on the CPU (two
# processes, cores 0-3, generic laws) with GMRES and BoomerAMG or MUMPS, and
# on the GPU (one process, cuda laws) with GMRES and BoomerAMG (hypre on the
# host) or cuDSS (hypre on the host or on the GPU); direct solvers with the
# reuse of the analysis, symmetric positive definite matrices; two passes.
# usage: measure-ex7.sh [cpu|gpu]; durations in durations-ex7.txt
P=$(cd "$(dirname "$0")/../.." && pwd)
source $P/scripts/env.sh > /dev/null
out=$P/results/cudss-2026-10-12
L=$P/build/cas-post/ex7/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
cd $P/build/examples-mumps/ex7
export OMP_NUM_THREADS=1 MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch \
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels \
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels \
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
EX7="./mox2 --no-post-processing --order 2 --mesh mesh/OneSphere.msh
  --nbsteps 40 --reference-file OneSphere-avgStress-40steps.ref"
DIRECT="-mra -ms SymmetricPositiveDefinite"
run() {
  local name=$1
  shift
  local start=$(date +%s.%N)
  "$@" > $out/ex7-$name.txt 2>&1
  local status=$?
  local end=$(date +%s.%N)
  LC_NUMERIC=C printf "%-28s %6.2f s, exit %s\n" "$name" "$(echo "$end - $start" | bc)" \
    "$status" | tee -a $out/durations-ex7.txt
}
for pass in 1 2; do
  if [ "${1:-cpu}" = cpu ]; then
    run np2-gmres-$pass taskset -c 0-3 mpirun -n 2 $EX7
    run np2-mumps-$pass taskset -c 0-3 mpirun -n 2 $EX7 -ls MUMPSSolver $DIRECT
  else
    run gpu-gmres-hypre-hote-$pass taskset -c 0-3 mpirun -n 1 $EX7 -l $L \
      --device cuda --hypre-on-host
    run gpu-cudss-hypre-hote-$pass taskset -c 0-3 mpirun -n 1 $EX7 -l $L \
      --device cuda --hypre-on-host -ls CuDSSSolver $DIRECT
    run gpu-cudss-hypre-gpu-$pass taskset -c 0-3 mpirun -n 1 $EX7 -l $L \
      --device cuda -ls CuDSSSolver $DIRECT
  fi
done
