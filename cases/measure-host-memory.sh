#!/bin/bash
# Effet des types de mémoire de MFEM sur le GPU. MFEM compilé avec Umpire
# alloue aussi la mémoire de l'hôte avec Umpire dès qu'un périphérique est
# configuré (HOST_UMPIRE) : chaque temporaire des boucles faites sur l'hôte
# (fonctions de forme, matrices élémentaires) passe alors par la table des
# allocations d'Umpire. Pour un cas, calculs sur le GPU (un processus)
# avec trois choix de mémoire :
# - umpire : MFEM_MEMORY=umpire, Umpire pour l'hôte et le GPU (choix de
#   MFEM) ;
# - host : MFEM_MEMORY=host, mémoire standard de l'hôte et du GPU
#   (cudaMalloc) ;
# - defaut : variable absente, choix de mfem-mgis : mémoire standard de
#   l'hôte et mémoire du GPU d'Umpire depuis la branche memoire-hote de
#   mfem-mgis, choix de MFEM (umpire) avant.
# Les sorties sont comparées à celles de la première variante pour la
# bulle, le cermet et le polycristal, pas pour la poutre ni pour ex8. Le
# polycristal, le cermet et ex8 tournent sans post-traitements, la bulle
# avec. À lancer sur une machine calme, GPU libre.
#
# usage : cases/measure-host-memory.sh cas [dossier de résultats]
# cas : bulle, cermet, polycristal, poutre ou ex8. Variables :
# - `exe` : dossier de l'exécutable (par défaut celui de
#   build/mm-opera-hpc-memoire, build/cas-memoire/poutre ou
#   build/examples-pointwise-ex8) ;
# - `cas` : dossier de construction de cases/, qui fournit les lois cuda
#   (par défaut build/cas-memoire) ;
# - `variantes` : choix de mémoire mesurés (par défaut « umpire host
#   defaut ») ;
# - `noyaux` : 1 (par défaut) pour tout calculer par noyaux avec les lois
#   cuda, 0 pour l'intégration point par point avec les lois `generic` ;
# - `mfem_mgis_sources` : voir cases/measure-common.sh ; à donner pour noter
#   la version de mfem-mgis des exécutables, par défaut celle du
#   sous-module mfem-mgis.
# Les chemins relatifs sont pris depuis le dossier courant.
lot=memoire
P=$(cd "$(dirname "$0")/.." && pwd)
case_name=${1:?usage : $0 cas [dossier de résultats]}
out=${2:-$P/results/memoire-hote-$case_name-$(date +%Y-%m-%d-%H%M)}
cas=$(cd "${cas:-$P/build/cas-memoire}" && pwd) || exit 1
source "$P/cases/measure-common.sh" "$case_name" "$out"
# exécutable et lois cuda du cas
case "$case_name" in
  bulle | cermet | polycristal)
    directory=$case_name
    [ "$case_name" = bulle ] && directory=bubble
    [ "$case_name" = polycristal ] && directory=polycrystal
    B=$(cd "${exe:-$P/build/mm-opera-hpc-memoire/$directory}" && pwd) || exit 1
    cuda_library=$cas/$case_name/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
    ;;
  poutre)
    B=$(cd "${exe:-$cas/poutre}" && pwd) || exit 1
    cuda_library=$cas/poutre/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
    ;;
  ex8)
    B=$(cd "${exe:-$P/build/examples-pointwise-ex8}" && pwd) || exit 1
    # lois de l'exemple, données par ses options
    cuda_library=""
    ;;
esac
generic=$B/src/libBehaviour.so
device=cuda
settings=()
library=$generic
if [ "${noyaux:-1}" = 1 ]; then
  settings=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
    MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
    MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
    MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
  if [ "$case_name" = ex8 ]; then
    options="-lU $cas/ex8/generic-parallel-cuda/src/libGenericParallelU3SI2-cuda.so
      -lA $cas/ex8/generic-parallel-cuda/src/libGenericParallelALFENI-cuda.so"
  else
    library=$cuda_library
  fi
else
  settings=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=PointByPoint)
fi
[ "$case_name" = poutre ] && executable=./poutre
first=""
for variant in ${variantes:-umpire host defaut}; do
  memory=()
  [ "$variant" != defaut ] && memory=(MFEM_MEMORY=$variant)
  name=$case_name-np1-gpu-$variant
  # référence des comparaisons : le premier calcul
  reference=${first:-$name}
  # pas de sorties d'un calcul précédent si celui-ci échoue avant de les
  # écrire
  [ "$case_name" = bulle ] &&
    rm -f "$B/sig_hydro.txt" "$B/bubbles_and_stresses_selected.txt"
  run "$name" 1 1 "$library" "${settings[@]}" "${memory[@]}"
  if [ "$case_name" = bulle ]; then
    # sorties écrites dans le dossier de la bulle
    for f in sig_hydro bubbles_and_stresses_selected; do
      cp "$B/$f.txt" "$out/$name-$f.txt" || continue
      if [ -n "$first" ]; then
        python3 "$P/cases/compare-outputs.py" --csv "$out/$first-$f.txt" \
          "$out/$name-$f.txt" --rtol 1e-4 | sed "s/^/   $f : /" |
          tee -a "$out/summary.txt"
      fi
    done
  fi
  first=${first:-$name}
done
echo "résultats dans $out"
