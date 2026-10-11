#!/bin/bash
# Mesures de bout en bout d'ex8 de mfem-mgis-examples (plaque
# thermomécanique : condition de Robin, couplage itératif, gonflement),
# point par point et tout par noyaux, maillage d'origine (8424 hexaèdres) :
# sur CPU, 2 processus d'un fil, comme le test de l'exemple ; sur le GPU, un
# processus (MUMPS résout la mécanique sur l'hôte), avec les lois de
# l'exemple compilées par l'interface `generic` (intégrées sur l'hôte) ou
# par l'interface `generic-parallel` pour le GPU (cas/ex8, construites dans
# `$cas/ex8`). L'exemple est celui de la construction donnée par la
# variable `ex8` (par défaut build/examples-pointwise-ex8, liée à
# install/mfem-mgis-pointwise). Le maillage raffiné une fois est hors de
# portée : plus de 20 min par pas avec MUMPS sur 8 processus, et 1 100 s de
# PCG avec BoomerAMG pour 5 pas sur 2 processus. À lancer sur une machine
# calme, GPU libre.
#
# usage : cas/mesures-ex8.sh [dossier de résultats]
lot=ex8
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/ex8-$(date +%Y-%m-%d-%H%M)}
source "$P/cas/mesures-commun.sh" ex8 "$out"
B=${ex8:-$P/build/examples-pointwise-ex8}
cas=${cas:-$P/build/cas-autres}
if [ ! -x "$B/rjh_plate" ]; then
  echo "$B/rjh_plate absent" >&2
  exit 1
fi
cuda_libraries="-lU $cas/ex8/generic-parallel-cuda/src/libGenericParallelU3SI2-cuda.so
  -lA $cas/ex8/generic-parallel-cuda/src/libGenericParallelALFENI-cuda.so"
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
refinement=0
run ex8-np2-point-par-point 1 2 ""
run ex8-np2-noyaux 1 2 "" "${kernels[@]}"
device=cuda
run ex8-np1-gpu-point-par-point 1 1 ""
run ex8-np1-gpu-noyaux 1 1 "" "${kernels[@]}"
options=$cuda_libraries
run ex8-np1-gpu-noyaux-lois-cuda 1 1 "" "${kernels[@]}"
echo "résultats dans $out"
