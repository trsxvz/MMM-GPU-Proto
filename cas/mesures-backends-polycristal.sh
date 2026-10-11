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
# Variables :
# - `polycristal` : dossier de l'exécutable de mm-opera-hpc (par défaut
#   build/mm-opera-hpc/polycrystal) ;
# - `backends` : backends mesurés, sous la forme <dossier>:<suffixe> (par
#   défaut ceux qui précèdent) ; pour kokkos (kokkos:kokkos), l'exécutable
#   doit être lié à un mfem-mgis construit avec Kokkos (option
#   enable-kokkos), qui l'initialise (variable MFEM_MGIS_KOKKOS) ;
# - `cas` : dossier de construction de cas/ (par défaut build/cas-autres) ;
# - `mfem_mgis_sources` : voir cas/mesures-commun.sh.
# Les chemins relatifs sont pris depuis le dossier courant.
#
# usage : cas/mesures-backends-polycristal.sh [dossier de résultats]
lot=backends
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/backends-polycristal-$(date +%Y-%m-%d-%H%M)}
source "$P/cas/mesures-commun.sh" polycristal "$out"
# chemins absolus : les calculs tournent dans le dossier de l'exécutable
cas=$(cd "${cas:-$P/build/cas-autres}" && pwd)
B=$(cd "${polycristal:-$B}" && pwd)
generic=$B/src/libBehaviour.so
backends=${backends:-cuda:cuda cuda-clang:cuda hip-clang-nvidia-platform:hip hip-hipcc:hip}
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
device=cuda
reference=polycristal-np1-gpu-generic-1
for pass in 1 2; do
  run "polycristal-np1-gpu-generic-$pass" 1 1 "$generic" "${kernels[@]}"
  for backend in $backends; do
    IFS=: read -r directory suffix <<< "$backend"
    settings=("${kernels[@]}")
    if [ "$suffix" = kokkos ]; then
      # Kokkos initialisé par mfem-mgis
      settings+=(MFEM_MGIS_KOKKOS=1)
    fi
    run "polycristal-np1-gpu-$directory-$pass" 1 1 \
      "$cas/polycristal/generic-parallel-$directory/src/libGenericParallelBehaviour-$suffix.so" \
      "${settings[@]}"
  done
done
echo "résultats dans $out"
