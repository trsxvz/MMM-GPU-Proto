#!/bin/bash
# Mesures de NewtonK0 (lot L5) sur le polycristal de mm-opera-hpc : Newton
# (PCG + BoomerAMG de hypre, nouvelle mise en place de BoomerAMG à chaque
# itération) contre NewtonK0 (méthode `NewtonK0` du solveur `FixedPoint` de
# mfem-mgis : gradient conjugué préconditionné par BoomerAMG construit une
# fois sur l'opérateur élastique, obtenu par une intégration de la loi).
# Tous les noyaux, sur 8 processus CPU (loi `generic`) puis sur le GPU (loi
# `cuda`). Les contraintes macroscopiques de NewtonK0 sont comparées à
# celles de Newton dans la même configuration. À lancer sur une machine
# calme, après build-mfem-mgis.sh.
#
# usage : cas/mesures-L5-polycristal.sh [dossier de résultats]
lot=L5
source "$(dirname "$0")/mesures-commun.sh" polycristal "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
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
  options=""
  reference=""
  run "polycristal-$where-newton" 1 $np "$library" "$batch" "${kernels[@]}"
  options="--non-linear-solver NewtonK0"
  reference="polycristal-$where-newton"
  run "polycristal-$where-newton-k0" 1 $np "$library" "$batch" \
    "${kernels[@]}"
done
echo "résultats dans $out"
