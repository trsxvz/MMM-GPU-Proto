#!/bin/bash
# Mesures de bout en bout des positions des entrées dans la matrice creuse
# (mfem-mgis, branche positions) : thermique r5 et poutre r4 sur 8
# processus, tangente assemblée par noyaux, avec la construction des cas
# liée à mfem-mgis installé dans install/mfem-mgis (avant) et à celui
# installé dans install/mfem-mgis-positions (après, build/cas-positions),
# plus le calcul point par point comme référence. À lancer sur une machine
# calme, rien d'autre ne tournant sur les cœurs (les processeurs N et N+8
# partagent un cœur).
#
# usage : cas/mesures-positions.sh [dossier de résultats]
lot=positions
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/positions-$(date +%Y-%m-%d-%H%M)}
batch=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
for case_name in thermique poutre; do
  source "$P/cas/mesures-commun.sh" "$case_name" "$out"
  if [ "$case_name" = thermique ]; then
    linear_solver=gmres-amg
    linear_tolerance=1e-8
  fi
  run "$case_name-np8-point-par-point" 1 8 "$generic"
  run "$case_name-np8-noyaux-avant" 1 8 "$generic" "${batch[@]}"
  B=$P/build/cas-positions/$case_name
  executable=./poutre
  run "$case_name-np8-noyaux-positions" 1 8 "$generic" "${batch[@]}"
done
echo "résultats dans $out"
