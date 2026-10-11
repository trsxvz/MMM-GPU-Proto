#!/bin/bash
# Mesures de bout en bout sur le GPU des inverses des jacobiennes calculées
# une fois (lot 7j) : poutre r4 et thermique r5 sur un processus et le GPU,
# tout par noyaux (loi `cuda`), avec la construction des cas liée à
# mfem-mgis installé dans install/mfem-mgis-registres (avant,
# build/cas-registres) et dans install/mfem-mgis-jacobiennes (après,
# build/cas-jacobiennes). Les bibliothèques des lois sont celles de
# build/cas. À lancer sur une machine calme, GPU libre.
#
# usage : cas/mesures-jacobiennes-gpu.sh [dossier de résultats]
lot=jacobiennes
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/jacobiennes-gpu-$(date +%Y-%m-%d-%H%M)}
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
for case_name in poutre thermique; do
  source "$P/cas/mesures-commun.sh" "$case_name" "$out"
  if [ "$case_name" = thermique ]; then
    linear_solver=gmres-amg
    linear_tolerance=1e-8
  fi
  cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
  device=cuda
  for version in registres jacobiennes; do
    B=$P/build/cas-$version/$case_name
    run "$case_name-np1-gpu-noyaux-$version" 1 1 "$cuda" "${kernels[@]}"
  done
done
echo "résultats dans $out"
