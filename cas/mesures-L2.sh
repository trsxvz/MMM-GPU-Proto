#!/bin/bash
# Mesures du lot L2 (lois intégrées par lots sur CPU) : intégration point
# par point (bibliothèque `generic`) contre intégration par lots
# (bibliothèques `generic` et `generic-parallel`, backends `openmp` et
# `stlpar`), à 1 et à 8 processus, puis avec plusieurs fils par processus.
# À lancer sur une machine calme, après build-mfem-mgis.sh.
#
# usage : cas/mesures-L2.sh [poutre|polycristal] [dossier de résultats]
#
# - poutre : poutre r4 du prototype (cas/poutre), von Mises à écrouissage
#   linéaire ;
# - polycristal : polycristal de mm-opera-hpc (5 grains, loi de monocristal
#   de l'UO2, orthotrope, en grandes transformations), 10 pas.
#
# Sorties brutes (avec les tableaux de profilage) dans le dossier de
# résultats (par défaut results/L2-<cas>-<date>), résumé dans resume.txt ;
# temps par phase : cas/resume-profilage.py <dossier>/*.txt. Pour le
# polycristal, les contraintes macroscopiques de chaque calcul sont
# comparées à celles du calcul point par point (cas/comparer-sorties.py).
set -u
P=$(cd "$(dirname "$0")/.." && pwd)
source "$P/scripts/env.sh" > /dev/null
# env.sh binds the threads of OpenMP (for Kokkos and AdaptiveCpp): libgomp,
# loaded by the libraries of the cases, then binds the main thread to the
# first CPU of the process, so that TBB (stlpar backend) sees one CPU, and
# `spread` over the hardware threads only uses half of the cores. The
# processes are placed by OpenMPI instead.
unset OMP_PROC_BIND OMP_PLACES
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
MPI=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
  location -i openmpi)
export PATH=$MPI/bin:$PATH
case_name=${1:-poutre}
out=${2:-$P/results/L2-$case_name-$(date +%Y-%m-%d-%H%M)}
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
# usage: command library integration name
# prints the command running the case
command_of() {
  local library=$1 integration=$2 name=$3
  case "$case_name" in
    poutre)
      echo ./poutre -r 4 -n 10 -mi 3000 -ls pcg-amg -lt 1e-4 -l "$library" \
        -bi "$integration" --profiling
      ;;
    polycristal)
      echo ./uniaxial-polycrystal --disable-post-processing --nstep 10 \
        --duration 10 -l "$library" \
        --macroscopic-stress-output-file "$out/$name.res"
      ;;
  esac
}
# usage: run name threads processes library integration
run() {
  local name=$1 threads=$2 np=$3 library=$4 integration=$5
  echo "== $name : $np processus, $threads fils, $integration" |
    tee -a "$out/resume.txt"
  # one thread: the default placement of OpenMPI (one core per process);
  # several threads: `threads` cores reserved for each process (with
  # `--bind-to none`, OpenMPI 5 binds all the processes to the first core)
  local placement=()
  if [ "$threads" != 1 ]; then
    placement=(--map-by "slot:PE=$threads" --bind-to core)
  fi
  local cmd
  read -r -a cmd <<< "$(command_of "$library" "$integration" "$name")"
  (cd "$B" && OMP_NUM_THREADS=$threads \
    MFEM_MGIS_BEHAVIOUR_INTEGRATION=$integration \
    timeout 3600 mpirun "${placement[@]}" -n "$np" "${cmd[@]}") \
    < /dev/null > "$out/$name.txt" 2>&1
  echo "   statut $?" | tee -a "$out/resume.txt"
  /usr/bin/grep -E "^dofs|^time" "$out/$name.txt" | tee -a "$out/resume.txt"
  # comparison with the integration point by point on the same number of
  # processes (the partition changes the results of the iterative solver)
  local reference=$out/$case_name-np$np-point-par-point.res
  if [ -f "$out/$name.res" ] && [ -f "$reference" ] &&
    [ "$out/$name.res" != "$reference" ]; then
    python3 "$P/cas/comparer-sorties.py" "$reference" \
      "$out/$name.res" --rtol 1e-8 | sed 's/^/   /' | tee -a "$out/resume.txt"
  fi
}
for np in 8 1; do
  run "$case_name-np$np-point-par-point" 1 "$np" "$generic" PointByPoint
  run "$case_name-np$np-lots-generic" 1 "$np" "$generic" Batch
  run "$case_name-np$np-lots-openmp-1fil" 1 "$np" "$openmp" Batch
done
# un processus sur les 8 cœurs physiques (le backend stlpar de g++ utilise
# TBB, qui prend tous les fils matériels des cœurs réservés)
run "$case_name-np1-lots-openmp-8fils" 8 1 "$openmp" Batch
run "$case_name-np1-lots-stlpar-8coeurs" 8 1 "$stlpar" Batch
# deux processus de quatre cœurs (et la référence point par point sur deux
# processus d'un cœur)
run "$case_name-np2-point-par-point" 1 2 "$generic" PointByPoint
run "$case_name-np2-lots-openmp-4fils" 4 2 "$openmp" Batch
echo "résultats dans $out"
