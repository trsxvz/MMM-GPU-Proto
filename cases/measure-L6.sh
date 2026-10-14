#!/bin/bash
# Mesures de la tranche GPU (lots L3 à L6) : calcul de référence sur CPU
# (8 processus, point par point et tout par noyaux), puis sur le GPU avec
# tous les noyaux, la loi `generic` (intégrée sur l'hôte) ou la loi `cuda`
# (intégrée sur le GPU), à 1 et 2 processus. À lancer sur une machine
# calme, après build-mfem-mgis.sh.
#
# usage : cases/measure-L6.sh [poutre|polycristal] [dossier de résultats]
# (voir cases/measure-common.sh)
lot=L6
source "$(dirname "$0")/measure-common.sh" "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
device=""
# référence des comparaisons du polycristal : calcul point par point sur
# le même nombre de processus (voir measure-common.sh)
run "$case_name-np8-point-par-point" 1 8 "$generic"
run "$case_name-np1-point-par-point" 1 1 "$generic"
run "$case_name-np2-point-par-point" 1 2 "$generic"
run "$case_name-np8-noyaux" 1 8 "$generic" "$batch" "${kernels[@]}"
device=cuda
run "$case_name-np1-gpu-noyaux-loi-generic" 1 1 "$generic" "$batch" \
  "${kernels[@]}"
run "$case_name-np1-gpu-noyaux-loi-cuda" 1 1 "$cuda" "$batch" "${kernels[@]}"
run "$case_name-np2-gpu-noyaux-loi-cuda" 1 2 "$cuda" "$batch" "${kernels[@]}"
echo "résultats dans $out"
