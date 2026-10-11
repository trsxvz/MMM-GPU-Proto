#!/bin/bash
# Mesures d'AmgX (lot L5) contre BoomerAMG sur le polycristal, sur le GPU :
# NewtonK0 (gradient conjugué préconditionné par un cycle du multigrille
# construit une fois sur l'opérateur élastique) avec la tolérance
# d'Eisenstat-Walker 1, le multigrille étant BoomerAMG ou AmgX (réglages
# par défaut de mfem-mgis, seuil 0,8 et interpolation MULTIPASS, en double
# ou en simple précision). Le calcul de Newton sert de référence aux
# contraintes macroscopiques. À lancer sur une machine calme, après
# build-mfem-mgis.sh.
#
# usage : cas/mesures-L5-polycristal-amgx.sh [dossier de résultats]
lot=L5
source "$(dirname "$0")/mesures-commun.sh" polycristal "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
device=cuda
k0="--non-linear-solver NewtonK0"
k0="$k0 --adaptive-linear-solver-tolerance EisenstatWalker1"
k0="$k0 --adaptive-linear-solver-initial-tolerance 0.1"
amgx="--elastic-operator-linear-solver AmgX"
tuned="--amgx-parameters strength_threshold=0.8,interpolator=MULTIPASS"
variants=(
  "newton:"
  "k0-amg:$k0"
  "k0-amgx-defaut:$k0 $amgx"
  "k0-amgx-regle:$k0 $amgx $tuned"
  "k0-amgx-regle-simple:$k0 $amgx $tuned --amgx-single-precision")
for variant in "${variants[@]}"; do
  IFS=: read -r name options <<< "$variant"
  reference=""
  if [ "$name" != newton ]; then
    reference="polycristal-np1-gpu-newton"
  fi
  run "polycristal-np1-gpu-$name" 1 1 "$cuda" "$batch" "${kernels[@]}"
done
echo "résultats dans $out"
