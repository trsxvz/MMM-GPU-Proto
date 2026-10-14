#!/bin/bash
# Bilan sur le GPU des lots du 11/10/2026 (positions, matrices et action sans
# former B, inverses des jacobiennes, matrice globale) : poutre r4 sur un
# processus et le GPU, tout par noyaux (loi `cuda`), tangente assemblée ou
# appliquée sans assemblage, Newton ou NewtonK0 avec Eisenstat-Walker, avec
# la construction des cas liée à mfem-mgis installé dans install/mfem-mgis
# (avant, build/cas) et dans install/mfem-mgis-globale (après,
# build/cas-globale). À lancer sur une machine calme, GPU libre.
#
# usage : cases/measure-gpu-summary.sh [dossier de résultats]
lot=bilan-gpu
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/bilan-gpu-$(date +%Y-%m-%d-%H%M)}
source "$P/cases/measure-common.sh" poutre "$out"
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels)
assembled=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
matrix_free=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
device=cuda
linear_tolerance=1e-4
k0="-nls fixed-point -fpm newton-k0 -els amg -k0 -pk0"
ew="-alst EisenstatWalker1 -alsti 0.1"
variants=(
  "assemblee-newton:pcg-amg:$assembled:"
  "assemblee-newton-k0-ew1:pcg-amg:$assembled:$k0 $ew"
  "sans-matrice-newton-k0-ew1:pcg-amg:$matrix_free:$k0 $ew")
for version in avant apres; do
  if [ "$version" = avant ]; then
    B=$P/build/cas/poutre
  else
    B=$P/build/cas-globale/poutre
  fi
  for variant in "${variants[@]}"; do
    IFS=: read -r name linear_solver tangent options <<< "$variant"
    run "poutre-np1-gpu-$name-$version" 1 1 "$cuda" "${kernels[@]}" \
      "$tangent"
  done
done
echo "résultats dans $out"
