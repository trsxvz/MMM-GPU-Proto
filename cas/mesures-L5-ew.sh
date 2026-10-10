#!/bin/bash
# Mesures du Newton inexact (lot L5) sur la poutre : tolérance du solveur
# linéaire fixe (1e-4) ou adaptée à la décroissance du résidu (choix 1 et 2
# d'Eisenstat et Walker, option -alst), pour Newton (PCG + BoomerAMG
# reconstruit à chaque itération) et pour NewtonK0 (gradient conjugué
# préconditionné par BoomerAMG construit une fois sur K0, prédiction par ce
# gradient conjugué). Sur 8 processus CPU puis sur le GPU (loi `cuda`). La
# tolérance fixe sert aussi aux prédictions. À lancer sur une machine calme,
# après build-mfem-mgis.sh.
#
# usage : cas/mesures-L5-ew.sh [dossier de résultats]
lot=L5
source "$(dirname "$0")/mesures-commun.sh" poutre "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
k0="-nls fixed-point -fpm newton-k0 -els amg -k0 -pk0"
# nom:tolérance du solveur linéaire:options
variants=(
  "newton:1e-4:"
  "newton-ew1:1e-4:-alst EisenstatWalker1"
  "newton-ew2:1e-4:-alst EisenstatWalker2"
  "newton-k0:1e-4:$k0"
  "newton-k0-ew1:1e-4:$k0 -alst EisenstatWalker1"
  "newton-k0-ew2:1e-4:$k0 -alst EisenstatWalker2")
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
    IFS=: read -r name linear_tolerance options <<< "$variant"
    run "poutre-$where-$name" 1 $np "$library" "$batch" "${kernels[@]}"
  done
done
echo "résultats dans $out"
