#!/bin/bash
# Fonctions communes aux scripts de mesures (cas/mesures-L*.sh), à sourcer
# après avoir donné le nom du lot dans la variable `lot` :
#
#   lot=L2
#   source "$(dirname "$0")/mesures-commun.sh" "$@"
#
# Arguments : cas (poutre, polycristal, thermique, ssna303, ex8, bulle ou
# cermet), dossier
# de résultats (par défaut results/<lot>-<cas>-<date>). Définit `run` (voir
# plus bas), `generic`, `openmp` et `stlpar` (bibliothèques de la loi du
# cas), et les variables `device`, `linear_solver`, `linear_tolerance`,
# `options` et `reference` lues par `run`, et `refinement` (nombre de
# raffinements du maillage : 4 pour la poutre, 5 pour la thermique, 0 pour
# ssna303 et ex8 par défaut). La variable `mfem_mgis_sources`, lue avant,
# donne les sources de mfem-mgis dont la version est notée dans resume.txt
# (par défaut le sous-module mfem-mgis), pour un exécutable construit
# contre une autre installation.
#
# - poutre : poutre r4 du prototype (cas/poutre), von Mises à écrouissage
#   linéaire ;
# - polycristal : polycristal de mm-opera-hpc (5 grains, loi de monocristal
#   de l'UO2, orthotrope, en grandes transformations), 10 pas ;
# - thermique : thermique stationnaire non linéaire dans la poutre
#   (cas/thermique), r5, 5 pas ;
# - ssna303 : barreau entaillé ssna303 en 3D (exemple ex4 de
#   mfem-mgis-examples, variante résolue par hypre), plasticité en grandes
#   transformations (déformation logarithmique de Hencky), 10 pas ;
# - ex8 : plaque thermomécanique (exemple ex8 de mfem-mgis-examples :
#   thermique avec condition de Robin, mécanique avec gonflement, couplage
#   itératif), 5 pas, MUMPS pour la mécanique, lois de l'exemple (la
#   bibliothèque donnée à `run` est ignorée) ;
# - bulle : bulle de mm-opera-hpc (cavité sphérique sous pression dans une
#   cellule périodique, élasticité, éléments quadratiques, maillage par
#   défaut single_sphere.msh, 467 109 inconnues, PCG et BoomerAMG), un pas,
#   post-traitements compris (export des contraintes aux nœuds, projection
#   L2) ; les sorties de la bulle sont écrites dans son dossier ;
# - cermet : cermet de mm-opera-hpc (polycristal de l'UO2 et interface
#   métallique), `nsteps` pas (1 par défaut) sur `duration` (0,4 s par
#   défaut).
# Les exécutables de la bulle et du cermet doivent avoir les options de
# mfem-mgis et l'option -l (mm-opera-hpc gpu-proto 88be07a ou après).
#
# Sorties brutes (avec les tableaux de profilage) dans le dossier de
# résultats, résumé dans resume.txt ; temps par phase :
# cas/resume-profilage.py <dossier>/*.txt. Pour le polycristal, les
# contraintes macroscopiques de chaque calcul sont comparées à celles du
# calcul de référence, par défaut le calcul point par point sur le même
# nombre de processus (cas/comparer-sorties.py).
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
  thermique)
    B=$P/build/cas/thermique
    generic=$B/src/libBehaviour.so
    gp=$B
    ;;
  ssna303)
    B=$P/build/mfem-mgis-examples/ex4
    generic=$B/src/libBehaviour.so
    gp=$P/build/cas/ssna303
    ;;
  ex8)
    B=$P/build/mfem-mgis-examples/ex8
    generic=""
    gp=$B
    ;;
  bulle)
    B=$P/build/mm-opera-hpc/bubble
    generic=$B/src/libBehaviour.so
    gp=$P/build/cas/bulle
    ;;
  cermet)
    B=$P/build/mm-opera-hpc/cermet
    generic=$B/src/libBehaviour.so
    gp=$P/build/cas/cermet
    ;;
  *)
    echo "usage : $0 [poutre|polycristal|thermique|ssna303|ex8|bulle|cermet] [dossier de résultats]" >&2
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
  echo "mfem-mgis : $(git -C "${mfem_mgis_sources:-$P/mfem-mgis}" log -1 --format=%h)"
} | tee -a "$out/resume.txt"
# périphérique de MFEM des calculs (option --device), vide pour le CPU
device=""
# poutre : solveur linéaire (option -ls), sa tolérance (option -lt) ;
# options supplémentaires du cas (par exemple `-pu 2` pour la poutre, voir
# cas/poutre/poutre.cxx, ou `--non-linear-solver NewtonK0` pour le
# polycristal)
linear_solver=pcg-amg
linear_tolerance=1e-4
options=""
# polycristal : nom du calcul de référence des comparaisons (par défaut, le
# calcul point par point sur le même nombre de processus)
reference=""
# lignes des sorties recopiées dans resume.txt (expression régulière
# étendue) : inconnues et temps de la boucle sur les pas, ou, pour ex8,
# statistiques des champs et temps du couplage
summary="^dofs|^time"
if [ "$case_name" = ex8 ]; then
  summary="^Temperature|^DisplacementNorm|IterativeCouplingScheme::computeNextState"
elif [ "$case_name" = bulle ]; then
  summary="^Newton iteration|^PCG Iterations"
fi
# poutre : exécutable (par défaut celui du dossier de construction, chemin
# absolu pour un autre exécutable, par exemple construit avec une version
# précédente de mfem-mgis)
executable=./poutre
# polycristal : option des post-traitements (vide pour les exécuter : export
# ParaView et moyennes des contraintes à chaque pas)
post_processing=--disable-post-processing
# usage : command_of bibliothèque nom
# affiche la commande du cas
command_of() {
  local library=$1 name=$2
  case "$case_name" in
    poutre)
      echo "$executable" -r "${refinement:-4}" -n 10 -mi 3000 -ls "$linear_solver" \
        -lt "$linear_tolerance" \
        -l "$library" --profiling ${device:+--device $device} $options
      ;;
    thermique)
      echo ./thermique -r "${refinement:-5}" -n 5 -ls "$linear_solver" \
        -lt "$linear_tolerance" \
        -l "$library" --profiling ${device:+--device $device} $options
      ;;
    ssna303)
      echo ./ssna303_3d_hypre --no-paraview --nbsteps 10 --end-time 1 \
        -r "${refinement:-0}" -l "$library" --profiling \
        ${device:+--device $device} $options
      ;;
    ex8)
      echo ./rjh_plate -m assemblage_hexa.msh -r "${refinement:-0}" -ns 5 \
        -lsMc MUMPSSolver -no-pp ${device:+--device $device} $options
      ;;
    bulle)
      echo ./test-bubble --mesh mesh/single_sphere.msh \
        --bubble-file mesh/single_bubble.txt \
        --linear-solver-preconditioner HypreBoomerAMG -v 1 \
        -l "$library" ${device:+--device $device} $options
      ;;
    cermet)
      echo ./cermet --duration "${duration:-0.4}" --nstep "${nsteps:-1}" \
        --postprocessing 0 -l "$library" \
        --macroscopic-stress-output-file "$out/$name.res" \
        ${device:+--device $device} $options
      ;;
    polycristal)
      echo ./uniaxial-polycrystal $post_processing --nstep 10 \
        --duration 10 -l "$library" \
        --macroscopic-stress-output-file "$out/$name.res" \
        ${device:+--device $device} $options
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
  if ! (cd "$B" && [ -x "${cmd[0]}" ]); then
    echo "   exécutable ${cmd[0]} absent de $B" | tee -a "$out/resume.txt"
    return 1
  fi
  (cd "$B" && env OMP_NUM_THREADS="$threads" "${settings[@]}" \
    timeout 3600 mpirun "${placement[@]}" -n "$np" "${cmd[@]}") \
    < /dev/null > "$out/$name.txt" 2>&1
  echo "   statut $?" | tee -a "$out/resume.txt"
  /usr/bin/grep -E "$summary" "$out/$name.txt" | tee -a "$out/resume.txt"
  # comparaison au calcul de référence, par défaut le calcul point par
  # point sur le même nombre de processus (la partition change les
  # résultats du solveur itératif)
  local reference_file=$out/${reference:-$case_name-np$np-point-par-point}.res
  if [ -f "$out/$name.res" ] && [ -f "$reference_file" ] &&
    [ "$out/$name.res" != "$reference_file" ]; then
    python3 "$P/cas/comparer-sorties.py" "$reference_file" \
      "$out/$name.res" --rtol 1e-8 --groupes 0 1-3 4-6 | sed 's/^/   /' |
      tee -a "$out/resume.txt"
  fi
}
