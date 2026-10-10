#!/bin/bash
# Mesures de la thermique stationnaire non linéaire par noyaux (lot L7) sur
# la poutre (cas/thermique, r5, 280 000 inconnues) : sur 8 processus CPU,
# point par point puis par noyaux ; sur le GPU (loi `cuda`), par noyaux,
# tangente assemblée (GMRES de hypre et BoomerAMG, tolérance fixe ou
# Eisenstat-Walker 1) ou appliquée sans assemblage (GMRES de MFEM et
# Jacobi). La tangente thermique est définie négative et les ddl imposés ont
# une diagonale à un : GMRES au lieu du gradient conjugué. À lancer sur une
# machine calme, après build-mfem-mgis.sh.
#
# usage : cas/mesures-L7-thermique.sh [dossier de résultats]
lot=L7
source "$(dirname "$0")/mesures-commun.sh" thermique "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
common=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels)
assembled=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
matrix_free=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_tolerance=1e-8
linear_solver=gmres-amg
run thermique-np8-point-par-point 1 8 "$generic"
run thermique-np8-noyaux 1 8 "$generic" "$batch" "${common[@]}" "$assembled"
device=cuda
run thermique-np1-gpu-noyaux 1 1 "$cuda" "$batch" "${common[@]}" "$assembled"
options="-alst EisenstatWalker1"
run thermique-np1-gpu-noyaux-ew1 1 1 "$cuda" "$batch" "${common[@]}" \
  "$assembled"
options=""
linear_solver=gmres-jacobi
run thermique-np1-gpu-sans-matrice 1 1 "$cuda" "$batch" "${common[@]}" \
  "$matrix_free"
echo "résultats dans $out"
