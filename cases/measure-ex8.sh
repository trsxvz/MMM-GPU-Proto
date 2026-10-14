#!/bin/bash
# End to end measurements of ex8 of mfem-mgis-examples (thermomechanical
# plate: Robin boundary condition, iterative coupling, swelling), point by
# point and everything by kernels, original mesh (8424 hexahedra): on the
# CPU, 2 processes of one thread, like the test of the example; on the GPU,
# one process (MUMPS solves the mechanics on the host), with the behaviours
# of the example compiled by the `generic` interface (integrated on the
# host) or by the `generic-parallel` interface for the GPU (cases/ex8, built
# in `$cases/ex8`). The example is the one of the build given by the
# variable `ex8` (by default build/examples-pointwise-ex8, linked to
# install/mfem-mgis-pointwise). The mesh refined once is out of reach:
# more than 20 min per step with MUMPS on 8 processes, and 1,100 s of PCG
# with BoomerAMG for 5 steps on 2 processes. To be run on a quiet machine,
# with the GPU free.
#
# usage: cases/measure-ex8.sh [results directory]
work_package=ex8
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/ex8-$(date +%Y-%m-%d-%H%M)}
source "$P/cases/measure-common.sh" ex8 "$out"
B=${ex8:-$P/build/examples-pointwise-ex8}
cases=${cases:-$P/build/cases-others}
if [ ! -x "$B/rjh_plate" ]; then
  echo "$B/rjh_plate missing" >&2
  exit 1
fi
cuda_libraries="-lU $cases/ex8/generic-parallel-cuda/src/libGenericParallelU3SI2-cuda.so
  -lA $cases/ex8/generic-parallel-cuda/src/libGenericParallelALFENI-cuda.so"
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
refinement=0
run ex8-np2-point-by-point 1 2 ""
run ex8-np2-kernels 1 2 "" "${kernels[@]}"
device=cuda
run ex8-np1-gpu-point-by-point 1 1 ""
run ex8-np1-gpu-kernels 1 1 "" "${kernels[@]}"
options=$cuda_libraries
run ex8-np1-gpu-kernels-behaviours-cuda 1 1 "" "${kernels[@]}"
echo "results in $out"
