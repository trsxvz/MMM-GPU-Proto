#!/bin/bash
# Comparaison des backends GPU de l'interface generic-parallel dans
# mfem-mgis : poutre r4 sur un processus et le GPU, tout par noyaux,
# tangente assemblée, Newton et PCG + BoomerAMG, avec la même loi générée
# pour chaque backend qui fonctionne dans mfem-mgis (cuda par nvcc et par
# clang, hip par clang et par hipcc sur la plateforme NVIDIA, kokkos avec
# Kokkos initialisé à l'option --kokkos, par le cas ou par mfem-mgis s'il
# est construit avec Kokkos). Les bibliothèques sont
# celles des dossiers generic-parallel-<dossier> de build/cas-autres/poutre
# (cuda pour nvcc, sinon le nom de la configuration), ou celles de
# `$libs/<configuration>/src` si la variable `libs` est donnée.
# L'exécutable est celui du dossier donné par la variable `poutre` (par
# défaut build/cas-autres/poutre), construit avec Kokkos pour l'option
# --kokkos (Kokkos trouvé par la configuration de cases/, ou mfem-mgis
# construit avec Kokkos). À lancer sur une machine calme, GPU libre.
#
# usage : cases/measure-backends-gpu.sh [dossier de résultats]
lot=backends
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/backends-gpu-$(date +%Y-%m-%d-%H%M)}
source "$P/cases/measure-common.sh" poutre "$out"
B=${poutre:-$P/build/cas-autres/poutre}
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
    run "poutre-np1-gpu-$configuration-$pass" 1 1 \
      "$(library_of "$configuration" "$directory" "$suffix")" "${kernels[@]}"
  done
done
echo "résultats dans $out"
