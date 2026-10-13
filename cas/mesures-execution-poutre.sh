#!/bin/bash
# Surcoût de chaque placement des calculs donné par l'option --execution de
# mfem-mgis (tâche 9a) : poutre r4 (111 843 inconnues), 10 pas, Newton,
# PCG de hypre et BoomerAMG, tout sur le CPU (`PROCESSUS_CPU` processus, 8
# par défaut, sur une machine calme), tout sur le GPU
# (un processus), et les huit placements des éléments finis, des lois et
# des solveurs linéaires en mode mixte (un processus et le GPU ; la poutre
# n'a pas de post-traitements). La loi est celle du cas, pour le CPU
# (interface generic) ou pour le GPU (interface generic-parallel, backend
# cuda), selon le placement des lois. Les copies des tableaux des matériaux
# entre l'hôte et le GPU sont relevées sous le tableau de profilage.
#
# La variable `cas_build` donne la construction des cas, construite contre
# un mfem-mgis qui a l'option --execution (par défaut build/cas-9a).
#
# usage : cas/mesures-execution-poutre.sh [dossier de résultats]
lot=execution
cas_build=${cas_build:-$(cd "$(dirname "$0")/.." && pwd)/build/cas-9a}
source "$(dirname "$0")/mesures-commun.sh" poutre "$@"
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
linear_tolerance=1e-4
summary="^dofs|^time|host to device|device to host"
device=""
# tout sur le CPU
options="--execution cpu"
np=${PROCESSUS_CPU:-8}
run "execution-cpu-${np}p" 1 "$np" "$generic"
# tout sur le GPU
options="--execution gpu"
run execution-gpu 1 1 "$cuda"
# placements mixtes
for fe in cpu gpu; do
  for bh in cpu gpu; do
    for lsd in cpu gpu; do
      library=$generic
      if [ "$bh" = gpu ]; then
        library=$cuda
      fi
      options="--execution mixed --finite-elements $fe --behaviours $bh --linear-solver-device $lsd"
      run "execution-mixte-fe-$fe-lois-$bh-solveur-$lsd" 1 1 "$library"
    done
  done
done
echo "résultats dans $out"
