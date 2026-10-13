#!/bin/bash
# Pression imposée sur l'extrémité de la poutre (lot L7.6). Avec les forces
# intérieures calculées par noyaux, la pression, indépendante des inconnues,
# est assemblée une fois pour une pression unité et ajoutée au résidu sur le
# GPU, sans intégrateur dans la forme non linéaire : la tangente est
# assemblée par noyaux sans passer par l'assemblage de MFEM, et peut être
# appliquée sans assemblage. Comparaison avec un exécutable construit avec
# la version précédente de mfem-mgis (premier argument), où la pression est
# un intégrateur de bord évalué sur l'hôte à chaque résidu et où les
# matrices élémentaires calculées par noyaux sont assemblées par MFEM sur
# l'hôte. Poutre r4, 10 pas, NewtonK0 (meilleure méthode du lot L5), tous
# les noyaux, sur le GPU et sur 8 processus CPU, avec et sans pression.
# Avant chaque calcul sur le GPU, le script attend au plus 10 minutes qu'un
# autre processus libère le GPU, puis le signale dans le résumé.
#
# La variable APRES donne l'exécutable de la nouvelle version (par défaut
# celui du dossier de construction du cas).
#
# usage : cas/mesures-pression-poutre.sh exécutable-avant [dossier de résultats]
lot=pression
before=$(realpath "$1")
shift
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
pressure=${PRESSION:-1e8}
base="-nls fixed-point -fpm newton-k0 -k0 -pk0 -alst EisenstatWalker1 -alsti 0.1 -els amg -amgo AggressiveCoarseningLevels=1"
echo "pression : $pressure Pa, exécutable avant : $before" | tee -a "$out/resume.txt"

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

for version in avant apres; do
  if [ "$version" = avant ]; then
    executable=$before
  else
    executable=${APRES:-./poutre}
  fi
  device=cuda
  options="$base -hdps 0.5"
  wait_for_gpu
  run "poutre-gpu-sans-pression-$version" 1 1 "$cuda" "$batch" "${kernels[@]}"
  options="$base -hdps 0.5 -pa $pressure"
  wait_for_gpu
  run "poutre-gpu-pression-$version" 1 1 "$cuda" "$batch" "${kernels[@]}"
  if [ "$version" = apres ]; then
    wait_for_gpu
    run "poutre-gpu-pression-sans-matrice-$version" 1 1 "$cuda" "$batch" \
      "${matrix_free[@]}"
  fi
  device=""
  options="$base"
  run "poutre-np8-sans-pression-$version" 1 8 "$generic" "$batch" \
    "${kernels[@]}"
  options="$base -pa $pressure"
  run "poutre-np8-pression-$version" 1 8 "$generic" "$batch" "${kernels[@]}"
done
python3 "$P/cas/resume-profilage.py" "$out"/poutre-*.txt > "$out/profilage.txt" 2>&1
echo "fin : $(date '+%F %T')" | tee -a "$out/resume.txt"
