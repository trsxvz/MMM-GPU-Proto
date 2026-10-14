#!/bin/bash
# Effect of the memory types of MFEM on the GPU. MFEM compiled with Umpire
# also allocates the memory of the host with Umpire as soon as a device is
# configured (HOST_UMPIRE): each temporary of the loops made on the host
# (shape functions, element matrices) then goes through the allocation
# table of Umpire. For one case, computations on the GPU (one process)
# with three choices of memory:
# - umpire: MFEM_MEMORY=umpire, Umpire for the host and the GPU (choice of
#   MFEM);
# - host: MFEM_MEMORY=host, standard memory of the host and of the GPU
#   (cudaMalloc);
# - default: variable not set, choice of mfem-mgis: standard memory of the
#   host and memory of the GPU from Umpire since the branch host-memory of
#   mfem-mgis, choice of MFEM (umpire) before.
# The outputs are compared to the ones of the first variant for the
# bubble, the cermet and the polycrystal, not for the beam nor for ex8. The
# polycrystal, the cermet and ex8 run without post-processings, the bubble
# with them. To be run on a quiet machine, with the GPU free.
#
# usage: cases/measure-host-memory.sh case [results directory]
# case: bubble, cermet, polycrystal, beam or ex8. Variables:
# - `exe`: directory of the executable (by default the one of
#   build/mm-opera-hpc-memory, build/cases-memory/beam or
#   build/examples-pointwise-ex8);
# - `cases`: build directory of cases/, which provides the cuda behaviours
#   (by default build/cases-memory);
# - `variants`: measured choices of memory (by default "umpire host
#   default");
# - `kernels`: 1 (default) to compute everything by kernels with the cuda
#   behaviours, 0 for the integration point by point with the `generic`
#   behaviours;
# - `mfem_mgis_sources`: see cases/measure-common.sh; to be given to record
#   the version of mfem-mgis of the executables, by default the one of the
#   submodule mfem-mgis.
# The relative paths are taken from the current directory.
work_package=memory
P=$(cd "$(dirname "$0")/.." && pwd)
case_name=${1:?usage: $0 case [results directory]}
out=${2:-$P/results/host-memory-$case_name-$(date +%Y-%m-%d-%H%M)}
cases=$(cd "${cases:-$P/build/cases-memory}" && pwd) || exit 1
source "$P/cases/measure-common.sh" "$case_name" "$out"
# executable and cuda behaviours of the case
case "$case_name" in
  bubble | cermet | polycrystal)
    B=$(cd "${exe:-$P/build/mm-opera-hpc-memory/$case_name}" && pwd) || exit 1
    cuda_library=$cases/$case_name/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
    ;;
  beam)
    B=$(cd "${exe:-$cases/beam}" && pwd) || exit 1
    cuda_library=$cases/beam/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
    ;;
  ex8)
    B=$(cd "${exe:-$P/build/examples-pointwise-ex8}" && pwd) || exit 1
    # behaviours of the example, given by its options
    cuda_library=""
    ;;
esac
generic=$B/src/libBehaviour.so
device=cuda
settings=()
library=$generic
if [ "${kernels:-1}" = 1 ]; then
  settings=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
    MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
    MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
    MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
  if [ "$case_name" = ex8 ]; then
    options="-lU $cases/ex8/generic-parallel-cuda/src/libGenericParallelU3SI2-cuda.so
      -lA $cases/ex8/generic-parallel-cuda/src/libGenericParallelALFENI-cuda.so"
  else
    library=$cuda_library
  fi
else
  settings=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=PointByPoint)
fi
[ "$case_name" = beam ] && executable=./beam
first=""
for variant in ${variants:-umpire host default}; do
  memory=()
  [ "$variant" != default ] && memory=(MFEM_MEMORY=$variant)
  name=$case_name-np1-gpu-$variant
  # reference of the comparisons: the first computation
  reference=${first:-$name}
  # no outputs of a previous computation if this one fails before writing
  # them
  [ "$case_name" = bubble ] &&
    rm -f "$B/sig_hydro.txt" "$B/bubbles_and_stresses_selected.txt"
  run "$name" 1 1 "$library" "${settings[@]}" "${memory[@]}"
  if [ "$case_name" = bubble ]; then
    # outputs written in the directory of the bubble
    for f in sig_hydro bubbles_and_stresses_selected; do
      cp "$B/$f.txt" "$out/$name-$f.txt" || continue
      if [ -n "$first" ]; then
        python3 "$P/cases/compare-outputs.py" --csv "$out/$first-$f.txt" \
          "$out/$name-$f.txt" --rtol 1e-4 | sed "s/^/   $f: /" |
          tee -a "$out/summary.txt"
      fi
    done
  fi
  first=${first:-$name}
done
echo "results in $out"
