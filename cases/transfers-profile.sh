#!/bin/bash
# Profil nsys des transferts entre l'hôte et le GPU (lot L7.3) : thermique
# r5 et poutre r4 par noyaux sur le GPU. Les piles d'appel des copies
# (--cudabacktrace) ne sont pas enregistrées sur cette machine (pas
# d'échantillonnage du processeur) : voir gdb-copies.txt dans les résultats.
# usage : cases/transfers-profile.sh dossier
# Analyse : nsys export -t sqlite, puis regroupement des copies
# (CUPTI_ACTIVITY_KIND_MEMCPY) par sens et par taille.
set -u
P=$(cd "$(dirname "$0")/.." && pwd)
out=$(mkdir -p "$1" && cd "$1" && pwd)
source $P/scripts/env.sh > /dev/null
unset OMP_PROC_BIND OMP_PLACES
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
nsys_cmd=(nsys profile -t cuda --cudabacktrace=memory -f true)
B=$P/build/cas/thermique
(cd "$B" && env OMP_NUM_THREADS=1 "${kernels[@]}" timeout -k 10 600 \
  "${nsys_cmd[@]}" -o "$out/thermique-r5" ./thermique -r 5 -n 5 \
  -ls gmres-amg -lt 1e-8 \
  -l $B/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so \
  --profiling --device cuda) < /dev/null > "$out/thermal-r5.txt" 2>&1
echo "thermique statut $?"
B=$P/build/cas/poutre
(cd "$B" && env OMP_NUM_THREADS=1 "${kernels[@]}" timeout -k 10 600 \
  "${nsys_cmd[@]}" -o "$out/poutre-r4" ./poutre -r 4 -n 10 -mi 3000 \
  -ls pcg-amg -lt 1e-4 -nls fixed-point -fpm newton-k0 -k0 -pk0 \
  -alst EisenstatWalker1 -alsti 0.1 -els amg -amgo AggressiveCoarseningLevels=1 \
  -l $B/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so \
  --profiling --device cuda) < /dev/null > "$out/beam-r4.txt" 2>&1
echo "poutre statut $?"
