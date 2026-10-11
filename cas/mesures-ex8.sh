#!/bin/bash
# Mesures de bout en bout d'ex8 de mfem-mgis-examples (plaque
# thermomécanique : condition de Robin, couplage itératif, gonflement),
# point par point et tout par noyaux, avec la construction séparée
# build/examples-autres liée à install/mfem-mgis-autres, maillage d'origine
# (8424 hexaèdres) : sur CPU, 2 processus d'un fil, comme le test de
# l'exemple ; sur le GPU, un processus (MUMPS résout la mécanique sur
# l'hôte). Le maillage raffiné une fois est hors de portée : plus de 20 min
# par pas avec MUMPS sur 8 processus, et 1 100 s de PCG avec BoomerAMG pour
# 5 pas sur 2 processus. À lancer sur une machine calme, GPU libre.
#
# usage : cas/mesures-ex8.sh [dossier de résultats]
lot=ex8
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/ex8-$(date +%Y-%m-%d-%H%M)}
source "$P/cas/mesures-commun.sh" ex8 "$out"
B=$P/build/examples-autres/ex8
if [ ! -x "$B/rjh_plate" ]; then
  echo "$B/rjh_plate absent" >&2
  exit 1
fi
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
echo "résultats dans $out"
