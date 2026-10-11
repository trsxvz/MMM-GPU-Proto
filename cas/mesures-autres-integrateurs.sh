#!/bin/bash
# Mesures de bout en bout des autres intégrateurs de la forme ajoutés aux
# calculs par noyaux (lot 7m) : thermique r5 sur 8 processus, échange
# convectif à la seconde extrémité (option -hc, intégrateur de faces de bord
# comme la condition de Robin de l'exemple ex8), tout par noyaux, tangente
# assemblée, avec la construction des cas liée à mfem-mgis installé dans
# install/mfem-mgis-globale (avant : boucles de MFEM sur tous les éléments,
# build/cas-globale) et dans install/mfem-mgis-autres (après,
# build/cas-autres), et, pour comparaison, le même calcul avec la
# température imposée (sans autre intégrateur). À lancer sur une machine
# calme, rien d'autre ne tournant sur les cœurs (les processeurs N et N+8
# partagent un cœur).
#
# usage : cas/mesures-autres-integrateurs.sh [dossier de résultats]
lot=autres
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/autres-integrateurs-$(date +%Y-%m-%d-%H%M)}
source "$P/cas/mesures-commun.sh" thermique "$out"
# la loi est celle de build/cas/thermique pour les deux versions (elle ne
# dépend pas de mfem-mgis) ; versions de mfem-mgis mesurées
for version in globale autres; do
  if [ ! -x "$P/build/cas-$version/thermique/thermique" ]; then
    echo "build/cas-$version/thermique/thermique absent" >&2
    exit 1
  fi
  source_dir=$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' \
    "$P/build/mfem-mgis-$version/CMakeCache.txt")
  commit=$(git -C "$source_dir" log -1 --format=%h 2> /dev/null ||
    echo "inconnu, $source_dir absent")
  echo "mfem-mgis-$version : $commit" | tee -a "$out/resume.txt"
done
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
linear_solver=gmres-amg
linear_tolerance=1e-8
for version in globale autres; do
  B=$P/build/cas-$version/thermique
  options="-hc 1"
  run "thermique-np8-robin-$version" 1 8 "$generic" "${kernels[@]}"
  options=""
  run "thermique-np8-dirichlet-$version" 1 8 "$generic" "${kernels[@]}"
done
echo "résultats dans $out"
