#!/bin/bash
# Mesures de bout en bout de la matrice globale avec motif et échanges
# gardés (lot 7k) : thermique r5 et poutre r4 sur 8 processus, tangente
# assemblée par noyaux, avec la construction des cas liée à mfem-mgis
# installé dans install/mfem-mgis-jacobiennes (avant, build/cas-jacobiennes)
# et dans install/mfem-mgis-globale (après, build/cas-globale). À lancer sur
# une machine calme, rien d'autre ne tournant sur les cœurs (les
# processeurs N et N+8 partagent un cœur).
#
# usage : cases/measure-global-matrix.sh [dossier de résultats]
lot=globale
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/matrice-globale-bout-en-bout-$(date +%Y-%m-%d-%H%M)}
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
for case_name in thermique poutre; do
  source "$P/cases/measure-common.sh" "$case_name" "$out"
  if [ "$case_name" = thermique ]; then
    linear_solver=gmres-amg
    linear_tolerance=1e-8
  fi
  for version in jacobiennes globale; do
    B=$P/build/cas-$version/$case_name
    run "$case_name-np8-noyaux-$version" 1 8 "$generic" "${kernels[@]}"
  done
done
echo "résultats dans $out"
