#!/bin/bash
# Mesures du barreau entaillé ssna303 en 3D (campagne de mesures, cas 2) :
# exemple ex4 de mfem-mgis-examples, variante résolue par hypre (Newton,
# FGMRES de hypre préconditionné par BoomerAMG), plasticité en grandes
# transformations (déformation logarithmique de Hencky), chargement complet
# en 10 pas, maillage raffiné zéro, une et deux fois (28 000 à 1,7 million
# d'inconnues) :
# - point par point puis par lots et noyaux sur 8 processus CPU ;
# - par noyaux sur le GPU (loi `cuda`), avec la réserve de mémoire de hypre
#   réduite à 0,5 Gio ;
# puis les mêmes calculs avec la formulation FBar (`--use-fbar`), par
# noyaux depuis le lot L7.9.
# Avant chaque calcul sur le GPU, le script attend au plus 10 minutes qu'un
# autre processus libère le GPU, puis le signale dans le résumé.
# À lancer sur une machine calme, après build-mfem-mgis.sh (examples et cas).
#
# usage : cas/mesures-ssna303.sh [dossier de résultats]
lot=ssna303
source "$(dirname "$0")/mesures-commun.sh" ssna303 "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so

# waits at most 10 minutes for the GPU to be free of other processes
wait_for_gpu() {
  local i
  for ((i = 0; i != 120; ++i)); do
    if [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; then
      return 0
    fi
    sleep 5
  done
  echo "   GPU partagé avec un autre processus" | tee -a "$out/resume.txt"
}

for refinement in 0 1 2; do
  for fbar in "" "--use-fbar"; do
    suffix=${fbar:+-fbar}
    device=""
    options="$fbar"
    run "ssna303-r$refinement$suffix-np8-point-par-point" 1 8 "$generic"
    run "ssna303-r$refinement$suffix-np8-noyaux" 1 8 "$generic" "$batch" \
      "${kernels[@]}"
    device=cuda
    options="$fbar -hdps 0.5"
    wait_for_gpu
    run "ssna303-r$refinement$suffix-gpu-noyaux" 1 1 "$cuda" "$batch" \
      "${kernels[@]}"
  done
done
echo "résultats dans $out"
