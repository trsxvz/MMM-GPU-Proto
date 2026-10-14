#!/bin/bash
# Coût des post-traitements (lot L8) : polycristal de mm-opera-hpc, 10 pas,
# tous les noyaux, Newton, sur le GPU et sur 8 processus CPU, sans puis
# avec les post-traitements du cas (export ParaView des résultats et
# moyennes des contraintes à chaque pas), calculés sur l'hôte puis par
# noyaux (MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels). La variable
# `executable` donne un autre exécutable que celui du dossier de
# construction (chemin absolu). Les sections de profilage de
# mfem-mgis séparent la recopie des matériaux vers l'hôte
# (NLEP::PP::synchronizeOnHost) de l'exécution des post-traitements
# (NLEP::PP::execute). Les fichiers ParaView sont écrits dans le dossier de
# construction du cas, puis effacés.
#
# usage : cases/measure-post-processings.sh [dossier de résultats]
lot=L8
source "$(dirname "$0")/measure-common.sh" polycristal "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
options=""

# waits at most 10 minutes for the GPU to be free of other processes
wait_for_gpu() {
  local i
  for ((i = 0; i != 120; ++i)); do
    if [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; then
      return 0
    fi
    sleep 5
  done
  echo "   GPU partagé avec un autre processus" | tee -a "$out/summary.txt"
}

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
  post_processing=--disable-post-processing
  reference=""
  [ -n "$device" ] && wait_for_gpu
  run "polycristal-$where-sans-post-traitements" 1 $np "$library" "$batch" \
    "${kernels[@]}"
  post_processing=""
  reference="polycristal-$where-sans-post-traitements"
  [ -n "$device" ] && wait_for_gpu
  run "polycristal-$where-post-traitements" 1 $np "$library" "$batch" \
    "${kernels[@]}"
  [ -n "$device" ] && wait_for_gpu
  run "polycristal-$where-post-traitements-noyaux" 1 $np "$library" \
    "$batch" "${kernels[@]}" MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels
  rm -rf "$B"/OutputFile-Uniaxial-polycrystal* "$B"/PolycrystalParaviewOutput
done
python3 "$P/cases/summary-profiling.py" "$out"/polycristal-*.txt > "$out/profiling.txt" 2>&1
echo "fin : $(date '+%F %T')" | tee -a "$out/summary.txt"
