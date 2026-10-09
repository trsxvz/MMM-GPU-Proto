#!/bin/bash
# Mesures de la tangente sans matrice (lot L4) sur la poutre : tangente
# assemblée par noyaux et PCG + BoomerAMG de hypre, contre la tangente
# appliquée sans assemblage par noyaux (MatrixFree) et le gradient
# conjugué de MFEM préconditionné par BoomerAMG construit sur la tangente
# assemblée, refaite toutes les 1, 2 ou 4 itérations de Newton (solveur
# cg-amg, option -pu). Sur 8 processus CPU puis sur le GPU (loi `cuda`).
# À lancer sur une machine calme, après build-mfem-mgis.sh.
#
# usage : cas/mesures-L4.sh [dossier de résultats]
lot=L4
source "$(dirname "$0")/mesures-commun.sh" poutre "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels)
assembled=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
matrix_free=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
for device in "" cuda; do
  if [ -z "$device" ]; then
    np=8
    library=$generic
    where=np8
  else
    np=1
    library=$cuda
    where=np1-gpu
  fi
  linear_solver=pcg-amg
  options=""
  run "poutre-$where-assemblee-pcg-amg" 1 $np "$library" "$batch" \
    "${kernels[@]}" "$assembled"
  linear_solver=cg-amg
  for period in 1 2 4; do
    options="-pu $period"
    run "poutre-$where-sans-matrice-cg-amg-pu$period" 1 $np "$library" \
      "$batch" "${kernels[@]}" "$matrix_free"
  done
done
echo "résultats dans $out"
