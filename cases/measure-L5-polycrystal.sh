#!/bin/bash
# Mesures de NewtonK0 et du Newton inexact (lot L5) sur le polycristal de
# mm-opera-hpc : Newton (PCG + BoomerAMG de hypre, nouvelle mise en place de
# BoomerAMG à chaque itération) contre NewtonK0 (méthode `NewtonK0` du
# solveur `FixedPoint` de mfem-mgis : gradient conjugué préconditionné par
# BoomerAMG construit une fois sur l'opérateur élastique, obtenu par une
# intégration de la loi), avec la tolérance du solveur linéaire du
# polycristal (1e-12) ou adaptée à la décroissance du résidu
# (Eisenstat-Walker 1, tolérance initiale 0,1). Tous les noyaux, sur 8
# processus CPU (loi `generic`) puis sur le GPU (loi `cuda`). Les
# contraintes macroscopiques sont comparées à celles de Newton dans la même
# configuration. À lancer sur une machine calme, après build-mfem-mgis.sh.
#
# usage : cases/measure-L5-polycrystal.sh [dossier de résultats]
lot=L5
source "$(dirname "$0")/measure-common.sh" polycristal "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
ew="--adaptive-linear-solver-tolerance EisenstatWalker1"
ew="$ew --adaptive-linear-solver-initial-tolerance 0.1"
# nom:options
variants=(
  "newton:"
  "newton-ew1:$ew"
  "newton-k0:--non-linear-solver NewtonK0"
  "newton-k0-ew1:--non-linear-solver NewtonK0 $ew")
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
  for variant in "${variants[@]}"; do
    IFS=: read -r name options <<< "$variant"
    reference=""
    if [ "$name" != newton ]; then
      reference="polycristal-$where-newton"
    fi
    run "polycristal-$where-$name" 1 $np "$library" "$batch" "${kernels[@]}"
  done
done
echo "résultats dans $out"
