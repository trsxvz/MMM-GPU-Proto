#!/bin/bash
# Coût des post-traitements (lot L8) : polycristal de mm-opera-hpc, 10 pas,
# tous les noyaux, Newton, sur le GPU et sur 8 processus CPU, sans puis
# avec les post-traitements du cas (export ParaView des résultats et
# moyennes des contraintes à chaque pas). Les sections de profilage de
# mfem-mgis séparent la recopie des matériaux vers l'hôte
# (NLEP::PP::synchronizeOnHost) de l'exécution des post-traitements
# (NLEP::PP::execute). Les fichiers ParaView sont écrits dans le dossier de
# construction du cas, puis effacés.
#
# usage : cas/mesures-post-traitements.sh [dossier de résultats]
lot=L8
source "$(dirname "$0")/mesures-commun.sh" polycristal "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
options=""
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
  run "polycristal-$where-sans-post-traitements" 1 $np "$library" "$batch" \
    "${kernels[@]}"
  post_processing=""
  reference="polycristal-$where-sans-post-traitements"
  run "polycristal-$where-post-traitements" 1 $np "$library" "$batch" \
    "${kernels[@]}"
  rm -rf "$B"/OutputFile-Uniaxial-polycrystal* "$B"/PolycrystalParaviewOutput
done
python3 "$P/cas/resume-profilage.py" "$out"/polycristal-*.txt > "$out/profilage.txt" 2>&1
echo "fin : $(date '+%F %T')" | tee -a "$out/resume.txt"
