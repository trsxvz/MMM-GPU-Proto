#!/bin/bash
# Comparaison des backends GPU de l'interface generic-parallel sur le
# polycristal de mm-opera-hpc (loi de monocristal de l'UO2, coûteuse) : un
# processus et le GPU, tout par noyaux, loi intégrée sur l'hôte (interface
# generic, référence des contraintes macroscopiques) ou sur le GPU avec la
# bibliothèque de chaque backend qui fonctionne sans modifier le programme
# (cuda par nvcc et par clang, hip par clang et par hipcc sur la plateforme
# NVIDIA), construites dans `$cas/polycristal/generic-parallel-<dossier>`
# (par défaut build/cas-autres). À lancer sur une machine calme, GPU libre.
#
# usage : cas/mesures-backends-polycristal.sh [dossier de résultats]
lot=backends
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/backends-polycristal-$(date +%Y-%m-%d-%H%M)}
source "$P/cas/mesures-commun.sh" polycristal "$out"
cas=${cas:-$P/build/cas-autres}
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
device=cuda
reference=polycristal-np1-gpu-generic-1
for pass in 1 2; do
  run "polycristal-np1-gpu-generic-$pass" 1 1 "$generic" "${kernels[@]}"
  for backend in cuda:cuda cuda-clang:cuda hip-clang-nvidia-platform:hip \
    hip-hipcc:hip; do
    IFS=: read -r directory suffix <<< "$backend"
    run "polycristal-np1-gpu-$directory-$pass" 1 1 \
      "$cas/polycristal/generic-parallel-$directory/src/libGenericParallelBehaviour-$suffix.so" \
      "${kernels[@]}"
  done
done
echo "résultats dans $out"
