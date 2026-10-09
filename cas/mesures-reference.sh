#!/bin/bash
# Temps de référence du chemin actuel de mfem-mgis (lot L0) : cas de la
# campagne de mesures, sur CPU (pas de périphérique), un fil OpenMP par
# processus. À lancer sur une machine calme, après build-mfem-mgis.sh.
#
# usage : cas/mesures-reference.sh [dossier de résultats]
#
# Sorties brutes (avec les tableaux de profilage) dans le dossier de
# résultats (par défaut results/reference-L0-<date>), une par cas, et un
# résumé dans resume.txt : commande, statut, durée totale du processus,
# lignes de bilan du cas.
set -u
P=$(cd "$(dirname "$0")/.." && pwd)
source "$P/scripts/env.sh" > /dev/null
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
MPI=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
  location -i openmpi)
export PATH=$MPI/bin:$PATH
export OMP_NUM_THREADS=1
out=${1:-$P/results/reference-L0-$(date +%Y-%m-%d-%H%M)}
mkdir -p "$out"
B=$P/build
# usage: run name directory command...
run() {
  local name=$1 dir=$2
  shift 2
  echo "== $name : $*" | tee -a "$out/resume.txt"
  local start end status
  start=$(date +%s.%N)
  (cd "$dir" && timeout 3600 "$@") < /dev/null > "$out/$name.txt" 2>&1
  status=$?
  end=$(date +%s.%N)
  echo "   statut $status, $(echo "$end - $start" | bc) s" |
    tee -a "$out/resume.txt"
  /usr/bin/grep -E "^dofs|^time" "$out/$name.txt" | tee -a "$out/resume.txt"
}
{
  echo "date : $(date '+%F %T')"
  echo "machine : $(nproc) cœurs, $(free -g | awk '/^Mem/ {print $2}') Go"
  echo "mfem-mgis : $(git -C "$P/mfem-mgis" log -1 --format=%h)"
} | tee -a "$out/resume.txt"

# 1. poutre en plasticité J2 (petites déformations, 3D), Newton de
# mfem-mgis : 15 795 inconnues (r3) à 4 processus, 111 843 (r4) à 8
poutre=$B/cas/poutre
for ls in "pcg-amg" "pcg-amg -lt 1e-4" "mumps"; do
  name=poutre-r3-np4-$(echo "$ls" | tr ' ' '_')
  run "$name" "$poutre" mpirun -n 4 ./poutre -r 3 -n 10 -mi 3000 -ls $ls \
    --profiling
done
for ls in "pcg-amg -lt 1e-4" "mumps"; do
  name=poutre-r4-np8-$(echo "$ls" | tr ' ' '_')
  run "$name" "$poutre" mpirun -n 8 ./poutre -r 4 -n 10 -mi 3000 -ls $ls \
    --profiling
done

# 2. ssna303 en 3D (grandes déformations), FGMRES et BoomerAMG de hypre,
# 10 pas jusqu'à la fin du chargement
ex4=$B/mfem-mgis-examples/ex4
run ssna303-3d-hypre-r0-np4 "$ex4" mpirun -n 4 ./ssna303_3d_hypre \
  -ns 10 -et 1 --profiling
run ssna303-3d-hypre-r1-np8 "$ex4" mpirun -n 8 ./ssna303_3d_hypre \
  -ns 10 -et 1 -r 1 --profiling

# 3. polycristal (5 grains, loi de monocristal), sans post-traitement
poly=$B/mm-opera-hpc/polycrystal
run polycristal-r0-np4 "$poly" mpirun -n 4 ./uniaxial-polycrystal \
  --disable-post-processing --nstep 10 --duration 10
run polycristal-r1-np8 "$poly" mpirun -n 8 ./uniaxial-polycrystal \
  --disable-post-processing --nstep 10 --duration 10 -r 1

# 4. thermomécanique (thermique non linéaire couplée), plaque de l'exemple
# ex8, maillage complet, MUMPS pour la mécanique
ex8=$B/mfem-mgis-examples/ex8
run thermomecanique-np4 "$ex8" mpirun -n 4 ./rjh_plate \
  -m assemblage_hexa.msh -lsMc MUMPSSolver
echo "résultats dans $out"
