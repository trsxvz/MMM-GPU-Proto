#!/bin/bash
# nsys profile of the transfers between the host and the GPU (work package
# L7.3): thermal case r5 and beam r4 by kernels on the GPU. The call stacks
# of the copies (--cudabacktrace) are not recorded on this machine (no
# sampling of the processor): see gdb-copies.txt in the results.
# usage: cases/transfers-profile.sh directory
# Analysis: nsys export -t sqlite, then grouping of the copies
# (CUPTI_ACTIVITY_KIND_MEMCPY) by direction and by size.
set -u
P=$(cd "$(dirname "$0")/.." && pwd)
out=$(mkdir -p "$1" && cd "$1" && pwd)
source $P/scripts/env.sh > /dev/null
unset OMP_PROC_BIND OMP_PLACES
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
nsys_cmd=(nsys profile -t cuda --cudabacktrace=memory -f true)
B=$P/build/cases/thermal
(cd "$B" && env OMP_NUM_THREADS=1 "${kernels[@]}" timeout -k 10 600 \
  "${nsys_cmd[@]}" -o "$out/thermal-r5" ./thermal -r 5 -n 5 \
  -ls gmres-amg -lt 1e-8 \
  -l $B/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so \
  --profiling --device cuda) < /dev/null > "$out/thermal-r5.txt" 2>&1
echo "thermal status $?"
B=$P/build/cases/beam
(cd "$B" && env OMP_NUM_THREADS=1 "${kernels[@]}" timeout -k 10 600 \
  "${nsys_cmd[@]}" -o "$out/beam-r4" ./beam -r 4 -n 10 -mi 3000 \
  -ls pcg-amg -lt 1e-4 -nls fixed-point -fpm newton-k0 -k0 -pk0 \
  -alst EisenstatWalker1 -alsti 0.1 -els amg -amgo AggressiveCoarseningLevels=1 \
  -l $B/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so \
  --profiling --device cuda) < /dev/null > "$out/beam-r4.txt" 2>&1
echo "beam status $?"
