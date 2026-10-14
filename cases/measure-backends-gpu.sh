#!/bin/bash
# Comparison of the GPU backends of the generic-parallel interface in
# mfem-mgis: beam r4 on one process and the GPU, everything by kernels,
# assembled tangent operator, Newton and PCG + BoomerAMG, with the same
# behaviour generated for each backend which works in mfem-mgis (cuda by
# nvcc and by clang, hip by clang and by hipcc on the NVIDIA platform,
# kokkos with Kokkos initialized by the option --kokkos, by the case or by
# mfem-mgis if it is built with Kokkos). The libraries are the ones of the
# directories generic-parallel-<directory> of build/cases-others/beam (cuda
# for nvcc, otherwise the name of the configuration), or the ones of
# `$libs/<configuration>/src` if the variable `libs` is given. The
# executable is the one of the directory given by the variable `beam` (by
# default build/cases-others/beam), built with Kokkos for the option --kokkos
# (Kokkos found by the configuration of cases/, or mfem-mgis built with
# Kokkos). To be run on a quiet machine, with the GPU free.
#
# usage: cases/measure-backends-gpu.sh [results directory]
work_package=backends
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/backends-gpu-$(date +%Y-%m-%d-%H%M)}
source "$P/cases/measure-common.sh" beam "$out"
B=${beam:-$P/build/cases-others/beam}
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
device=cuda
library_of() {
  local configuration=$1 directory=$2 suffix=$3
  if [ -n "${libs:-}" ]; then
    echo "$libs/$configuration/src/libGenericParallelBehaviour-$suffix.so"
  else
    echo "$B/generic-parallel-$directory/src/libGenericParallelBehaviour-$suffix.so"
  fi
}
for pass in 1 2; do
  for backend in cuda-nvcc:cuda:cuda cuda-clang:cuda-clang:cuda \
    hip-clang-nvidia-platform:hip-clang-nvidia-platform:hip \
    hip-hipcc:hip-hipcc:hip kokkos:kokkos:kokkos; do
    IFS=: read -r configuration directory suffix <<< "$backend"
    options=""
    if [ "$suffix" = kokkos ]; then
      options="--kokkos"
    fi
    run "beam-np1-gpu-$configuration-$pass" 1 1 \
      "$(library_of "$configuration" "$directory" "$suffix")" "${kernels[@]}"
  done
done
echo "results in $out"
