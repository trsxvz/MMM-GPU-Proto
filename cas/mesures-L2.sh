#!/bin/bash
# Mesures du lot L2 (lois intégrées par lots sur CPU) : intégration point
# par point (bibliothèque `generic`) contre intégration par lots
# (bibliothèques `generic` et `generic-parallel`, backends `openmp` et
# `stlpar`), à 1 et à 8 processus, puis avec plusieurs fils par processus.
# À lancer sur une machine calme, après build-mfem-mgis.sh.
#
# usage : cas/mesures-L2.sh [poutre|polycristal] [dossier de résultats]
# (voir cas/mesures-commun.sh)
lot=L2
source "$(dirname "$0")/mesures-commun.sh" "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
for np in 8 1; do
  run "$case_name-np$np-point-par-point" 1 "$np" "$generic"
  run "$case_name-np$np-lots-generic" 1 "$np" "$generic" "$batch"
  run "$case_name-np$np-lots-openmp-1fil" 1 "$np" "$openmp" "$batch"
done
# un processus sur les 8 cœurs physiques (le backend stlpar de g++ utilise
# TBB, qui prend tous les fils matériels des cœurs réservés)
run "$case_name-np1-lots-openmp-8fils" 8 1 "$openmp" "$batch"
run "$case_name-np1-lots-stlpar-8coeurs" 8 1 "$stlpar" "$batch"
# deux processus de quatre cœurs (et la référence point par point sur deux
# processus d'un cœur)
run "$case_name-np2-point-par-point" 1 2 "$generic"
run "$case_name-np2-lots-openmp-4fils" 4 2 "$openmp" "$batch"
echo "résultats dans $out"
