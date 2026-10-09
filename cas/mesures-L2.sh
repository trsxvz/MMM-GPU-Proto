#!/bin/bash
# Mesures du lot L2 (lois intégrées par lots sur CPU) sur la poutre :
# intégration point par point (bibliothèque `generic`) contre intégration
# par lots (bibliothèques `generic` et `generic-parallel`, backends
# `openmp` et `stlpar`), à 1 et à 8 processus. À lancer sur une machine
# calme, après build-mfem-mgis.sh.
#
# usage : cas/mesures-L2.sh [dossier de résultats]
#
# Sorties brutes (avec les tableaux de profilage) dans le dossier de
# résultats (par défaut results/L2-<date>), résumé dans resume.txt ; temps
# par phase : cas/resume-profilage.py <dossier>/*.txt
set -u
P=$(cd "$(dirname "$0")/.." && pwd)
source "$P/scripts/env.sh" > /dev/null
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
MPI=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
  location -i openmpi)
export PATH=$MPI/bin:$PATH
out=${1:-$P/results/L2-$(date +%Y-%m-%d-%H%M)}
mkdir -p "$out"
B=$P/build/cas/poutre
generic=$B/src/libBehaviour.so
openmp=$B/generic-parallel-openmp/src/libGenericParallelBehaviour-openmp.so
stlpar=$B/generic-parallel-stlpar/src/libGenericParallelBehaviour-stlpar.so
{
  echo "date : $(date '+%F %T')"
  echo "machine : $(nproc --all) cœurs, $(free -g | awk '/^Mem/ {print $2}') Go"
  echo "mfem-mgis : $(git -C "$P/mfem-mgis" log -1 --format=%h)"
} | tee -a "$out/resume.txt"
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
  (cd "$B" && OMP_NUM_THREADS=$threads \
    timeout 3600 mpirun "${placement[@]}" -n "$np" ./poutre -r 4 -n 10 \
    -mi 3000 -ls pcg-amg -lt 1e-4 -l "$library" -bi "$integration" \
    --profiling) < /dev/null > "$out/$name.txt" 2>&1
  echo "   statut $?" | tee -a "$out/resume.txt"
  /usr/bin/grep -E "^dofs|^time" "$out/$name.txt" | tee -a "$out/resume.txt"
}
for np in 8 1; do
  run poutre-np$np-point-par-point 1 "$np" "$generic" PointByPoint
  run poutre-np$np-lots-generic 1 "$np" "$generic" Batch
  run poutre-np$np-lots-openmp-1fil 1 "$np" "$openmp" Batch
done
# un processus sur les 8 cœurs physiques (le backend stlpar de g++ utilise
# TBB, qui prend tous les fils matériels des cœurs réservés)
run poutre-np1-lots-openmp-8fils 8 1 "$openmp" Batch
run poutre-np1-lots-stlpar-8coeurs 8 1 "$stlpar" Batch
# deux processus de quatre cœurs
run poutre-np2-lots-openmp-4fils 4 2 "$openmp" Batch
echo "résultats dans $out"
