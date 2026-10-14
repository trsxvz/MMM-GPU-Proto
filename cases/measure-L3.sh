#!/bin/bash
# Mesures du lot L3 (calculs par noyaux) : intégration par lots avec les
# gradients, les forces intérieures et les matrices élémentaires de la
# tangente calculés élément par élément ou par noyaux, à 1 et à 8
# processus, sur CPU. À lancer sur une machine calme, après
# build-mfem-mgis.sh.
#
# usage : cases/measure-L3.sh [poutre|polycristal] [dossier de résultats]
# (voir cases/measure-common.sh)
lot=L3
source "$(dirname "$0")/measure-common.sh" "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
gradients=MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
forces=MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
tangent=MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
for np in 8 1; do
  run "$case_name-np$np-point-par-point" 1 "$np" "$generic"
  run "$case_name-np$np-lots" 1 "$np" "$generic" "$batch"
  run "$case_name-np$np-lots-noyaux-gradients" 1 "$np" "$generic" \
    "$batch" "$gradients"
  run "$case_name-np$np-lots-noyaux-forces" 1 "$np" "$generic" \
    "$batch" "$forces"
  run "$case_name-np$np-lots-noyaux-tangente" 1 "$np" "$generic" \
    "$batch" "$tangent"
  run "$case_name-np$np-lots-noyaux" 1 "$np" "$generic" \
    "$batch" "$gradients" "$forces" "$tangent"
done
echo "résultats dans $out"
