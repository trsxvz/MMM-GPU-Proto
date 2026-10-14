#!/bin/bash
# Mesures de la tangente appliquée sans assemblage (lot L5) sur la poutre,
# sur le GPU (loi `cuda`) : gradient conjugué de MFEM préconditionné par le
# lisseur de Jacobi ou de Chebyshev (ordres 2 à 5), et NewtonK0 (gradient
# conjugué préconditionné par BoomerAMG construit une fois sur K0, avec ou
# sans Eisenstat-Walker 1), comparés à NewtonK0 avec la tangente assemblée.
# À lancer sur une machine calme, après build-mfem-mgis.sh.
#
# usage : cases/measure-L5-matrix-free.sh [dossier de résultats]
lot=L5
source "$(dirname "$0")/measure-common.sh" poutre "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
common=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels)
assembled=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
matrix_free=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
device=cuda
linear_tolerance=1e-4
k0="-nls fixed-point -fpm newton-k0 -els amg -k0 -pk0"
ew="-alst EisenstatWalker1 -alsti 0.1"
# nom:solveur linéaire:tangente:options
variants=(
  "sans-matrice-cg-jacobi:cg-jacobi:$matrix_free:"
  "sans-matrice-cg-chebyshev-2:cg-chebyshev:$matrix_free:-cho 2"
  "sans-matrice-cg-chebyshev-3:cg-chebyshev:$matrix_free:-cho 3"
  "sans-matrice-cg-chebyshev-4:cg-chebyshev:$matrix_free:-cho 4"
  "sans-matrice-cg-chebyshev-5:cg-chebyshev:$matrix_free:-cho 5"
  "sans-matrice-newton-k0:pcg-amg:$matrix_free:$k0"
  "sans-matrice-newton-k0-ew1:pcg-amg:$matrix_free:$k0 $ew"
  "assemblee-newton-k0-ew1:pcg-amg:$assembled:$k0 $ew")
for variant in "${variants[@]}"; do
  IFS=: read -r name linear_solver tangent options <<< "$variant"
  run "poutre-np1-gpu-$name" 1 1 "$cuda" "$batch" "${common[@]}" "$tangent"
done
echo "résultats dans $out"
