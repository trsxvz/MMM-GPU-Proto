#!/bin/bash
# Passage à l'échelle (campagne de mesures, cas 5) : poutre raffinée trois,
# quatre puis cinq fois (16 000 à 840 000 inconnues, 10 pas), meilleure
# méthode du lot L5 (NewtonK0, gradient conjugué préconditionné par le
# multigrille de K0, prédiction, tolérance d'Eisenstat-Walker 1), tous les
# noyaux :
# - sur le GPU, multigrille BoomerAMG agressif ou AmgX réglé en simple
#   précision, tangente assemblée, et BoomerAMG agressif avec la tangente
#   appliquée sans assemblage (seul K0 est assemblé), avec la mémoire
#   maximale occupée sur le GPU (relevée par nvidia-smi toutes les 200 ms,
#   autres processus compris) ;
# - sur 8 processus CPU, BoomerAMG agressif.
# Sur le GPU, la réserve de mémoire de hypre (Umpire, 4 Gio par défaut,
# allouée d'un coup) est réduite à 0,5 Gio (option -hdps de mfem-mgis) :
# temps à quelques % près, 3,5 Gio de moins. Le raffinement suivant
# (6,5 millions d'inconnues) ne tient pas dans les 16 Go du GPU, même ainsi
# (allocation refusée vers 15,3 Gio, tangente assemblée ou non).
# Avant chaque calcul sur le GPU, le script attend au plus 10 minutes qu'un
# autre processus libère le GPU, puis le signale dans le résumé.
# À lancer sur une machine calme, après build-mfem-mgis.sh.
#
# usage : cas/mesures-echelle-poutre.sh [dossier de résultats]
lot=echelle
source "$(dirname "$0")/mesures-commun.sh" poutre "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
matrix_free=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
linear_tolerance=1e-4
base="-nls fixed-point -fpm newton-k0 -k0 -pk0 -alst EisenstatWalker1 -alsti 0.1"
pool="-hdps 0.5"
amg="-els amg -amgo AggressiveCoarseningLevels=1"
amgx="-els amgx -amgxs -amgxo strength_threshold=0.8,interpolator=MULTIPASS"

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

# run on the GPU with the peak of the memory used on the GPU
# usage: run_gpu name settings...
run_gpu() {
  local name=$1
  shift
  wait_for_gpu
  nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -lms 200 \
    > "$out/$name.memoire.txt" &
  local sampler=$!
  run "$name" 1 1 "$cuda" "$@"
  kill "$sampler"
  wait "$sampler" 2> /dev/null
  echo "   mémoire du GPU au plus $(sort -n "$out/$name.memoire.txt" |
    tail -1) Mio" | tee -a "$out/resume.txt"
}

for refinement in 3 4 5; do
  device=cuda
  options="$base $pool $amg"
  run_gpu "poutre-r$refinement-gpu-amg-agressif" "$batch" "${kernels[@]}"
  options="$base $pool $amgx"
  run_gpu "poutre-r$refinement-gpu-amgx-simple" "$batch" "${kernels[@]}"
  options="$base $pool $amg"
  run_gpu "poutre-r$refinement-gpu-amg-agressif-sans-matrice" "$batch" \
    "${matrix_free[@]}"
  device=""
  options="$base $amg"
  run "poutre-r$refinement-np8-amg-agressif" 1 8 "$generic" "$batch" \
    "${kernels[@]}"
done
echo "résultats dans $out"
