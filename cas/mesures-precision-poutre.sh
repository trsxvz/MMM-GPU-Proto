#!/bin/bash
# Précision float contre double sur le CPU (tâche 9b) : poutre du prototype
# (cas/poutre), `RAFFINEMENTS` raffinements (4 par défaut, 111 843
# inconnues), `PAS` pas (10), Newton, PCG de hypre et BoomerAMG (tolérance
# 1e-4), `PROCESSUS` processus (2) attachés aux cœurs `CPUS` (0,1), chaque
# calcul répété `REPETITIONS` fois (3).
#
# Les deux constructions ne diffèrent que par la précision : environnements
# spack jumeaux deps/spack/mmm-double-cpu et deps/spack/mmm-float-cpu (MFEM,
# hypre et MUMPS en double ou en simple, sans CUDA), MGIS en double
# (install/mgis) ou en float (install/mgis-float), mfem-mgis construit
# contre chacun, cas construits contre chaque mfem-mgis (build/cas-double-cpu
# et build/cas-float, ou les dossiers donnés par `CAS_DOUBLE` et
# `CAS_FLOAT`).
#
# Critères de Newton : celui du cas (1e-8) en double seulement, comme
# référence, puis `CRITERE` (1e-4 par défaut) pour les deux précisions, pour
# comparer les temps à critère égal. En float, le critère du cas, borné à
# 100 ε = 1,2e-5 par mfem-mgis, n'est pas atteint à r4 : le résidu relatif
# plafonne vers 1,3e-5 et Newton s'arrête au bout de 50 itérations.
#
# usage : cas/mesures-precision-poutre.sh [dossier de résultats]
set -u
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/precision-poutre-$(date +%Y-%m-%d-%H%M)}
mkdir -p "$out"
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
np=${PROCESSUS:-2}
cpus=${CPUS:-0,1}
repetitions=${REPETITIONS:-3}
refinements=${RAFFINEMENTS:-4}
steps=${PAS:-10}
criterion=${CRITERE:-1e-4}
declare -A cas=([double]=${CAS_DOUBLE:-$P/build/cas-double-cpu}
                [float]=${CAS_FLOAT:-$P/build/cas-float})
declare -A mgis=([double]=$P/install/mgis [float]=$P/install/mgis-float)
declare -A environment=([double]=$P/deps/spack/mmm-double-cpu
                        [float]=$P/deps/spack/mmm-float-cpu)
{
  echo "poutre r$refinements, $steps pas, $np processus sur les cœurs $cpus,"
  echo "$repetitions répétitions ; $(date)"
  echo "machine : $(nproc) fils ; charge : $(cut -d' ' -f1-3 /proc/loadavg)"
} | tee "$out/resume.txt"
# run <précision> <nom> <répétition> <options...>
run() {
  local precision=$1 name=$2 i=$3
  shift 3
  local mpi
  mpi=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "${environment[$precision]}" \
    location -i mpi)
  local B=${cas[$precision]}/poutre
  if [ ! -x "$B/poutre" ]; then
    echo "== $name : exécutable $B/poutre absent" | tee -a "$out/resume.txt"
    return 1
  fi
  echo "== $name-$i : $*" | tee -a "$out/resume.txt"
  # environnement propre : seuls MPI, TFEL et le MGIS de la précision ;
  # taskset autour de mpirun, qui attache chaque processus à un cœur de
  # l'ensemble autorisé
  (cd "$B" && env -i HOME="$HOME" PATH="$mpi/bin:/usr/bin:/bin" \
    LD_LIBRARY_PATH="$P/install/tfel-mmm-gpu/lib:${mgis[$precision]}/lib" \
    taskset -c "$cpus" timeout 3600 mpirun -n "$np" ./poutre \
    -r "$refinements" -n "$steps" -ls pcg-amg -lt 1e-4 -prof "$@") \
    < /dev/null > "$out/$name-$i.txt" 2>&1
  echo "   statut $?" | tee -a "$out/resume.txt"
  /usr/bin/grep -E "^dofs|^time" "$out/$name-$i.txt" | sed 's/^/   /' |
    tee -a "$out/resume.txt"
}
# les précisions alternent à chaque répétition, pour qu'une dérive de la
# machine touche les deux de la même façon
for i in $(seq 1 "$repetitions"); do
  run double double-critere-du-cas "$i"
  for precision in double float; do
    run "$precision" "$precision-critere-$criterion" "$i" -rt "$criterion"
  done
done
echo "résultats dans $out"
