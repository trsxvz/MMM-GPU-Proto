#!/bin/bash
# Comparison of the GPU backends of the generic-parallel interface on the
# polycrystal of mm-opera-hpc (single crystal behaviour of UO2, expensive):
# one process and the GPU, everything by kernels, behaviour integrated on
# the host (generic interface, reference of the macroscopic stresses) or on
# the GPU with the library of each backend which works without changing the
# program (cuda by nvcc and by clang, hip by clang and by hipcc on the
# NVIDIA platform), built in `$cases/polycrystal/generic-parallel-<directory>`
# (by default build/cases-others). To be run on a quiet machine, with the GPU
# free.
#
# Variables:
# - `polycrystal`: directory of the executable of mm-opera-hpc (by default
#   build/mm-opera-hpc/polycrystal);
# - `backends`: measured backends, in the form <directory>:<suffix> (by
#   default the ones above); for kokkos (kokkos:kokkos), the executable
#   must be linked to a mfem-mgis built with Kokkos (option
#   enable-kokkos), which initializes it (variable MFEM_MGIS_KOKKOS);
# - `cases`: build directory of cases/ (by default build/cases-others);
# - `mfem_mgis_sources`: see cases/measure-common.sh.
# The relative paths are taken from the current directory.
#
# usage: cases/measure-backends-polycrystal.sh [results directory]
work_package=backends
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/backends-polycrystal-$(date +%Y-%m-%d-%H%M)}
source "$P/cases/measure-common.sh" polycrystal "$out"
# absolute paths: the computations run in the directory of the executable
cases=$(cd "${cases:-$P/build/cases-others}" && pwd)
B=$(cd "${polycrystal:-$B}" && pwd)
generic=$B/src/libBehaviour.so
backends=${backends:-cuda:cuda cuda-clang:cuda hip-clang-nvidia-platform:hip hip-hipcc:hip}
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
device=cuda
reference=polycrystal-np1-gpu-generic-1
for pass in 1 2; do
  run "polycrystal-np1-gpu-generic-$pass" 1 1 "$generic" "${kernels[@]}"
  for backend in $backends; do
    IFS=: read -r directory suffix <<< "$backend"
    settings=("${kernels[@]}")
    if [ "$suffix" = kokkos ]; then
      # Kokkos initialized by mfem-mgis
      settings+=(MFEM_MGIS_KOKKOS=1)
    fi
    run "polycrystal-np1-gpu-$directory-$pass" 1 1 \
      "$cases/polycrystal/generic-parallel-$directory/src/libGenericParallelBehaviour-$suffix.so" \
      "${settings[@]}"
  done
done
echo "results in $out"
