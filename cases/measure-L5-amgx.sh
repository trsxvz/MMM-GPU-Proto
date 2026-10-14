#!/bin/bash
# Mesures d'AmgX (lot L5) contre BoomerAMG sur la poutre :
# 1. opérateur élastique K0 seul (poutre-k0) : mise en place et
#    résolutions du gradient conjugué préconditionné par chaque
#    multigrille, sur le GPU ;
# 2. calcul complet par la meilleure méthode, NewtonK0 (gradient conjugué
#    préconditionné par un cycle du multigrille construit une fois sur K0,
#    prédiction par ce gradient conjugué) avec la tolérance
#    d'Eisenstat-Walker 1, le multigrille de K0 étant BoomerAMG ou AmgX,
#    sur le GPU puis sur 8 processus CPU (AmgX calcule alors sur le GPU,
#    partagé par les 8 processus).
# Paramètres d'AmgX passés par -amgxo (noms d'AmgX) : réglages par défaut
# de mfem-mgis (ceux de MFEM : seuil de connexion 0,25, interpolation D2),
# seuil 0,8 et interpolation MULTIPASS (meilleur réglage du balayage sur
# K0), avec deux balayages du lisseur, et en simple précision (-amgxs).
# Sur 8 processus, chaque processus a ses ressources AmgX sur le même GPU.
# À lancer sur une machine calme, après build-mfem-mgis.sh.
#
# usage : cases/measure-L5-amgx.sh [dossier de résultats]
lot=L5
source "$(dirname "$0")/measure-common.sh" poutre "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
amgx_tuned=strength_threshold=0.8,interpolator=MULTIPASS
# 1. opérateur K0 seul, sur le GPU (nom|options de poutre-k0)
k0_variants=(
  "amg|-s cg-amg"
  "amg-agressif|-s cg-amg -amgo AggressiveCoarseningLevels=1"
  "pcg-amg|-s pcg-amg"
  "amgx-defaut|-s cg-amgx"
  "amgx-regle|-s cg-amgx -amgxo $amgx_tuned"
  "amgx-regle-balayages2|-s cg-amgx -amgxo $amgx_tuned,presweeps=2,postsweeps=2"
  "amgx-regle-simple|-s cg-amgx -amgxs -amgxo $amgx_tuned"
  "amgx-defaut-simple|-s cg-amgx -amgxs"
  "amgx-pcg-regle|-s amgx-pcg -amgxo $amgx_tuned"
  "jacobi|-s cg-jacobi")
for variant in "${k0_variants[@]}"; do
  IFS='|' read -r name k0_options <<< "$variant"
  echo "== k0-$name" | tee -a "$out/summary.txt"
  read -r -a k0_cmd <<< "./poutre-k0 -r 4 -n 5 --device cuda $k0_options"
  (cd "$B" && env OMP_NUM_THREADS=1 timeout 1800 "${k0_cmd[@]}") \
    < /dev/null > "$out/k0-$name.txt" 2>&1
  echo "   statut $?" | tee -a "$out/summary.txt"
  /usr/bin/grep -E "^dofs|resolution" "$out/k0-$name.txt" |
    tee -a "$out/summary.txt"
done
# 2. calcul complet
linear_solver=pcg-amg
linear_tolerance=1e-4
base="-nls fixed-point -fpm newton-k0 -k0 -pk0 -alst EisenstatWalker1 -alsti 0.1"
# nom:options de la poutre
variants=(
  "amg:-els amg"
  "amg-agressif:-els amg -amgo AggressiveCoarseningLevels=1"
  "amgx-defaut:-els amgx"
  "amgx-regle:-els amgx -amgxo $amgx_tuned"
  "amgx-regle-balayages2:-els amgx -amgxo $amgx_tuned,presweeps=2,postsweeps=2"
  "amgx-regle-simple:-els amgx -amgxs -amgxo $amgx_tuned")
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
    IFS=: read -r name extra <<< "$variant"
    options="$base $extra"
    run "poutre-$where-$name" 1 $np "$library" "$batch" "${kernels[@]}"
  done
done
echo "résultats dans $out"
