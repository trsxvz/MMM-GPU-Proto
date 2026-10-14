#!/bin/bash
# Mesures des options de BoomerAMG (lot L5) sur la poutre, pour la meilleure
# méthode : NewtonK0 (gradient conjugué préconditionné par BoomerAMG
# construit une fois sur K0, prédiction par ce gradient conjugué) avec la
# tolérance d'Eisenstat-Walker 1 (tolérance initiale 0,1). Options de hypre
# passées par -amgo (seuil de connexion, grossissement agressif, lisseur,
# interpolation, nombre de balayages). Sur le GPU (loi `cuda`) puis sur 8
# processus CPU. À lancer sur une machine calme, après build-mfem-mgis.sh.
#
# usage : cases/measure-L5-amg.sh [dossier de résultats]
lot=L5
source "$(dirname "$0")/measure-common.sh" poutre "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
linear_tolerance=1e-4
base="-nls fixed-point -fpm newton-k0 -els amg -k0 -pk0 -alst EisenstatWalker1 -alsti 0.1"
# nom:options de BoomerAMG (vide : options par défaut de mfem-mgis)
variants=(
  "defaut:"
  "seuil-05:StrengthThreshold=0.5"
  "seuil-07:StrengthThreshold=0.7"
  "agressif-1:AggressiveCoarseningLevels=1"
  "hmis:CoarseningType=10"
  "jacobi:RelaxationType=7"
  "chebyshev:RelaxationType=16"
  "balayages-2:RelaxationSweeps=2"
  "interp-ext:InterpolationType=14"
  "seuil-05-agressif-1:StrengthThreshold=0.5,AggressiveCoarseningLevels=1")
for device in cuda ""; do
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
    IFS=: read -r name amg <<< "$variant"
    options="$base${amg:+ -amgo $amg}"
    run "poutre-$where-amg-$name" 1 $np "$library" "$batch" "${kernels[@]}"
  done
done
echo "résultats dans $out"
