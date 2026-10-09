#!/bin/bash
# Fonctions communes aux scripts de mesures (cas/mesures-L*.sh), à sourcer
# après avoir donné le nom du lot dans la variable `lot` :
#
#   lot=L2
#   source "$(dirname "$0")/mesures-commun.sh" "$@"
#
# Arguments : cas (poutre ou polycristal), dossier de résultats (par défaut
# results/<lot>-<cas>-<date>). Définit `run` (voir plus bas), `generic`,
# `openmp` et `stlpar` (bibliothèques de la loi du cas).
#
# - poutre : poutre r4 du prototype (cas/poutre), von Mises à écrouissage
#   linéaire ;
# - polycristal : polycristal de mm-opera-hpc (5 grains, loi de monocristal
#   de l'UO2, orthotrope, en grandes transformations), 10 pas.
#
# Sorties brutes (avec les tableaux de profilage) dans le dossier de
# résultats, résumé dans resume.txt ; temps par phase :
# cas/resume-profilage.py <dossier>/*.txt. Pour le polycristal, les
# contraintes macroscopiques de chaque calcul sont comparées à celles du
# calcul point par point sur le même nombre de processus
# (cas/comparer-sorties.py).
set -u
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$P/scripts/env.sh" > /dev/null
# env.sh lie les fils d'OpenMP (pour Kokkos et AdaptiveCpp) : libgomp,
# chargée par les bibliothèques des cas, attache alors le fil principal au
# premier processeur du processus, si bien que TBB (backend stlpar) ne voit
# qu'un processeur, et `spread` sur les fils matériels n'occupe que la
# moitié des cœurs. Le placement est laissé à OpenMPI.
unset OMP_PROC_BIND OMP_PLACES
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
MPI=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
  location -i openmpi)
export PATH=$MPI/bin:$PATH
case_name=${1:-poutre}
out=${2:-$P/results/$lot-$case_name-$(date +%Y-%m-%d-%H%M)}
case "$case_name" in
  poutre)
    B=$P/build/cas/poutre
    generic=$B/src/libBehaviour.so
    gp=$B
    ;;
  polycristal)
    B=$P/build/mm-opera-hpc/polycrystal
    generic=$B/src/libBehaviour.so
    gp=$P/build/cas/polycristal
    ;;
  *)
    echo "usage : $0 [poutre|polycristal] [dossier de résultats]" >&2
    exit 1
    ;;
esac
openmp=$gp/generic-parallel-openmp/src/libGenericParallelBehaviour-openmp.so
stlpar=$gp/generic-parallel-stlpar/src/libGenericParallelBehaviour-stlpar.so
mkdir -p "$out"
# absolute path: the cases run in their build directories
out=$(cd "$out" && pwd)
{
  echo "date : $(date '+%F %T')"
  echo "cas : $case_name"
  echo "machine : $(nproc --all) cœurs, $(free -g | awk '/^Mem/ {print $2}') Go"
  echo "charge : $(cut -d ' ' -f 1-3 /proc/loadavg)"
  echo "mfem-mgis : $(git -C "$P/mfem-mgis" log -1 --format=%h)"
} | tee -a "$out/resume.txt"
# périphérique de MFEM des calculs (option --device), vide pour le CPU
device=""
# usage : command_of bibliothèque nom
# affiche la commande du cas
command_of() {
  local library=$1 name=$2
  case "$case_name" in
    poutre)
      echo ./poutre -r 4 -n 10 -mi 3000 -ls pcg-amg -lt 1e-4 -l "$library" \
        --profiling ${device:+--device $device}
      ;;
    polycristal)
      echo ./uniaxial-polycrystal --disable-post-processing --nstep 10 \
        --duration 10 -l "$library" \
        --macroscopic-stress-output-file "$out/$name.res"
      ;;
  esac
}
# usage : run nom fils processus bibliothèque [VARIABLE=valeur...]
# Les variables d'environnement données choisissent les options de
# mfem-mgis (MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch, ...).
run() {
  local name=$1 threads=$2 np=$3 library=$4
  shift 4
  local settings=("$@")
  echo "== $name : $np processus, $threads fils, ${device:-cpu}, ${settings[*]:-défaut}" |
    tee -a "$out/resume.txt"
  # un fil : placement par défaut d'OpenMPI (un cœur par processus) ;
  # plusieurs fils : `threads` cœurs réservés à chaque processus (avec
  # `--bind-to none`, OpenMPI 5 attache tous les processus au premier cœur)
  local placement=()
  if [ "$threads" != 1 ]; then
    placement=(--map-by "slot:PE=$threads" --bind-to core)
  fi
  local cmd
  read -r -a cmd <<< "$(command_of "$library" "$name")"
  (cd "$B" && env OMP_NUM_THREADS="$threads" "${settings[@]}" \
    timeout 3600 mpirun "${placement[@]}" -n "$np" "${cmd[@]}") \
    < /dev/null > "$out/$name.txt" 2>&1
  echo "   statut $?" | tee -a "$out/resume.txt"
  /usr/bin/grep -E "^dofs|^time" "$out/$name.txt" | tee -a "$out/resume.txt"
  # comparaison au calcul point par point sur le même nombre de processus
  # (la partition change les résultats du solveur itératif)
  local reference=$out/$case_name-np$np-point-par-point.res
  if [ -f "$out/$name.res" ] && [ -f "$reference" ] &&
    [ "$out/$name.res" != "$reference" ]; then
    python3 "$P/cas/comparer-sorties.py" "$reference" \
      "$out/$name.res" --rtol 1e-8 --groupes 0 1-3 4-6 | sed 's/^/   /' |
      tee -a "$out/resume.txt"
  fi
}
