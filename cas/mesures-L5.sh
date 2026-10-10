#!/bin/bash
# Mesures du point fixe (lot L5) sur la poutre : meilleur Newton (tangente
# assemblée par noyaux, PCG + BoomerAMG de hypre, tolérance 1e-4) contre le
# point fixe sur l'opérateur élastique K0, accéléré (Anderson de profondeur
# 8, Cast3M), K0 résolu par un cycle de BoomerAMG gardé pour tout le calcul
# (lisseur par défaut ou de Chebyshev), ou par PCG + BoomerAMG gardé, et
# contre le Newton modifié accéléré. Pour le point fixe, la prédiction est
# résolue à 1e-8 près (sinon les pas élastiques demandent des centaines
# d'itérations). Sur 8 processus CPU puis sur le GPU (loi `cuda`). À lancer
# sur une machine calme, après build-mfem-mgis.sh.
#
# usage : cas/mesures-L5.sh [dossier de résultats]
lot=L5
source "$(dirname "$0")/mesures-commun.sh" poutre "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
# nom:tolérance du solveur linéaire:options
fp="-nls fixed-point"
variants=(
  "newton:1e-4:"
  "point-fixe-anderson8-amg-k0:1e-8:$fp -acc anderson -ad 8 -els amg -k0"
  "point-fixe-anderson8-amg-chebyshev-k0:1e-8:$fp -acc anderson -ad 8 -els amg -k0 -amgo RelaxationType=16"
  "point-fixe-cast3m-amg-k0:1e-8:$fp -acc cast3m -els amg -k0"
  "point-fixe-anderson8-pcg-amg-k0:1e-8:$fp -acc anderson -ad 8 -els pcg-amg -k0"
  "newton-modifie-anderson8:1e-8:$fp -fpm modified-newton -acc anderson -ad 8")
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
