#!/bin/bash
# Functions shared by the measurement scripts (cases/measure-L*.sh), to be
# sourced after giving the name of the work package in the variable
# `work_package`:
#
#   work_package=L2
#   source "$(dirname "$0")/measure-common.sh" "$@"
#
# Arguments: case (beam, polycrystal, thermal, ssna303, ex8, bubble or
# cermet), results directory (by default
# results/<work_package>-<case>-<date>). Defines `run`
# (see below), `generic`, `openmp` and `stlpar` (libraries of the behaviour
# of the case), and the variables `device`, `linear_solver`,
# `linear_tolerance`, `options` and `reference` read by `run`, and
# `refinement` (number of refinements of the mesh: 4 for the beam, 5 for the
# thermal case, 0 for ssna303 and ex8 by default). The variable
# `mfem_mgis_sources`, read before, gives the sources of mfem-mgis whose
# version is written in summary.txt (by default the submodule mfem-mgis),
# for an executable built against another installation. The variable
# `examples_build`, read before, gives the build of mfem-mgis-examples used
# for ssna303 and ex8 (by default build/mfem-mgis-examples). The variable
# `cases_build`, read before, gives the build of the cases of the prototype
# used for the beam and the thermal case (by default build/cases).
#
# - beam: beam r4 of the prototype (cases/beam), von Mises with linear
#   hardening;
# - polycrystal: polycrystal of mm-opera-hpc (5 grains, single crystal
#   behaviour of UO2, orthotropic, finite strain), 10 steps;
# - thermal: stationary non linear heat transfer in the beam
#   (cases/thermal), r5, 5 steps;
# - ssna303: notched bar ssna303 in 3D (example ex4 of
#   mfem-mgis-examples, variant solved by hypre), finite strain plasticity
#   (logarithmic strain of Hencky), 10 steps;
# - ex8: thermomechanical plate (example ex8 of mfem-mgis-examples: heat
#   transfer with a Robin boundary condition, mechanics with swelling,
#   iterative coupling), 5 steps, MUMPS for the mechanics, behaviours of the
#   example (the library given to `run` is ignored);
# - bubble: bubble of mm-opera-hpc (pressurized spherical cavity in a
#   periodic cell, elasticity, quadratic elements, default mesh
#   single_sphere.msh, 467,109 unknowns, PCG and BoomerAMG), one step,
#   post-processings included (export of the stresses at the nodes, L2
#   projection); the outputs of the bubble are written in its directory;
# - cermet: cermet of mm-opera-hpc (polycrystal of UO2 and metallic
#   interface), `nsteps` steps (1 by default) over `duration` (0.4 s by
#   default).
# The executables of the bubble and of the cermet must have the options of
# mfem-mgis and the option -l (mm-opera-hpc gpu-proto 88be07a or later).
#
# Raw outputs (with the profiling tables) in the results directory, summary
# in summary.txt; time per phase:
# cases/summary-profiling.py <directory>/*.txt. For the polycrystal, the
# macroscopic stresses of each computation are compared to the ones of the
# reference computation, by default the computation point by point on the
# same number of processes (cases/compare-outputs.py).
set -u
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$P/scripts/env.sh" > /dev/null
# env.sh binds the OpenMP threads (for Kokkos and AdaptiveCpp): libgomp,
# loaded by the libraries of the cases, then binds the main thread to the
# first processor of the process, so that TBB (stlpar backend) only sees one
# processor, and `spread` over the hardware threads only uses half of the
# cores. The placement is left to OpenMPI.
unset OMP_PROC_BIND OMP_PLACES
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
MPI=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
  location -i openmpi)
export PATH=$MPI/bin:$PATH
case_name=${1:-beam}
out=${2:-$P/results/$work_package-$case_name-$(date +%Y-%m-%d-%H%M)}
case "$case_name" in
  beam)
    B=${cases_build:-$P/build/cases}/beam
    generic=$B/src/libBehaviour.so
    gp=$B
    ;;
  polycrystal)
    B=$P/build/mm-opera-hpc/polycrystal
    generic=$B/src/libBehaviour.so
    gp=$P/build/cases/polycrystal
    ;;
  thermal)
    B=${cases_build:-$P/build/cases}/thermal
    generic=$B/src/libBehaviour.so
    gp=$B
    ;;
  ssna303)
    B=${examples_build:-$P/build/mfem-mgis-examples}/ex4
    generic=$B/src/libBehaviour.so
    gp=$P/build/cases/ssna303
    ;;
  ex8)
    B=${examples_build:-$P/build/mfem-mgis-examples}/ex8
    generic=""
    gp=$B
    ;;
  bubble)
    B=$P/build/mm-opera-hpc/bubble
    generic=$B/src/libBehaviour.so
    gp=$P/build/cases/bubble
    ;;
  cermet)
    B=$P/build/mm-opera-hpc/cermet
    generic=$B/src/libBehaviour.so
    gp=$P/build/cases/cermet
    ;;
  *)
    echo "usage: $0 [beam|polycrystal|thermal|ssna303|ex8|bubble|cermet] [results directory]" >&2
    exit 1
    ;;
esac
openmp=$gp/generic-parallel-openmp/src/libGenericParallelBehaviour-openmp.so
stlpar=$gp/generic-parallel-stlpar/src/libGenericParallelBehaviour-stlpar.so
mkdir -p "$out"
# absolute path: the cases run in their build directories
out=$(cd "$out" && pwd)
{
  echo "date: $(date '+%F %T')"
  echo "case: $case_name"
  echo "machine: $(nproc --all) cores, $(free -g | awk '/^Mem/ {print $2}') GB"
  echo "load: $(cut -d ' ' -f 1-3 /proc/loadavg)"
  echo "mfem-mgis: $(git -C "${mfem_mgis_sources:-$P/mfem-mgis}" log -1 --format=%h)"
} | tee -a "$out/summary.txt"
# device of MFEM for the computations (option --device), empty for the CPU
device=""
# beam: linear solver (option -ls), its tolerance (option -lt);
# additional options of the case (for instance `-pu 2` for the beam, see
# cases/beam/beam.cxx, or `--non-linear-solver NewtonK0` for the
# polycrystal)
linear_solver=pcg-amg
linear_tolerance=1e-4
options=""
# polycrystal: name of the reference computation of the comparisons (by
# default, the computation point by point on the same number of processes)
reference=""
# lines of the outputs copied into summary.txt (extended regular
# expression): unknowns and time of the loop over the steps, or, for ex8,
# statistics of the fields and time of the coupling
summary="^dofs|^time"
if [ "$case_name" = ex8 ]; then
  summary="^Temperature|^DisplacementNorm|IterativeCouplingScheme::computeNextState"
elif [ "$case_name" = bubble ]; then
  summary="^Newton iteration|^PCG Iterations"
fi
# beam and polycrystal: executable (by default the one of the build
# directory, absolute path for another executable, for instance built with
# another version of mfem-mgis), read in the variable `executable` if it is
# defined before
executable=${executable:-}
# polycrystal: option of the post-processings (empty to run them: ParaView
# export and averages of the stresses at each step)
post_processing=--disable-post-processing
# usage: command_of library name
# prints the command of the case
command_of() {
  local library=$1 name=$2
  case "$case_name" in
    beam)
      echo "${executable:-./beam}" -r "${refinement:-4}" -n 10 -mi 3000 \
        -ls "$linear_solver" \
        -lt "$linear_tolerance" \
        -l "$library" --profiling ${device:+--device $device} $options
      ;;
    thermal)
      echo ./thermal -r "${refinement:-5}" -n 5 -ls "$linear_solver" \
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
    bubble)
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
    polycrystal)
      echo "${executable:-./uniaxial-polycrystal}" $post_processing --nstep 10 \
        --duration 10 -l "$library" \
        --macroscopic-stress-output-file "$out/$name.res" \
        ${device:+--device $device} $options
      ;;
  esac
}
# usage: run name threads processes library [VARIABLE=value...]
# The given environment variables select the options of mfem-mgis
# (MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch, ...).
run() {
  local name=$1 threads=$2 np=$3 library=$4
  shift 4
  local settings=("$@")
  echo "== $name: $np processes, $threads threads, ${device:-cpu}, ${settings[*]:-default}" |
    tee -a "$out/summary.txt"
  # one thread: default placement of OpenMPI (one core per process);
  # several threads: `threads` cores reserved for each process (with
  # `--bind-to none`, OpenMPI 5 binds all the processes to the first core)
  local placement=()
  if [ "$threads" != 1 ]; then
    placement=(--map-by "slot:PE=$threads" --bind-to core)
  fi
  local cmd
  read -r -a cmd <<< "$(command_of "$library" "$name")"
  if ! (cd "$B" && [ -x "${cmd[0]}" ]); then
    echo "   executable ${cmd[0]} missing from $B" | tee -a "$out/summary.txt"
    return 1
  fi
  (cd "$B" && env OMP_NUM_THREADS="$threads" "${settings[@]}" \
    timeout 3600 mpirun "${placement[@]}" -n "$np" "${cmd[@]}") \
    < /dev/null > "$out/$name.txt" 2>&1
  echo "   status $?" | tee -a "$out/summary.txt"
  /usr/bin/grep -E "$summary" "$out/$name.txt" | tee -a "$out/summary.txt"
  # comparison to the reference computation, by default the computation
  # point by point on the same number of processes (the partition changes
  # the results of the iterative solver)
  local reference_file=$out/${reference:-$case_name-np$np-point-by-point}.res
  if [ -f "$out/$name.res" ] && [ -f "$reference_file" ] &&
    [ "$out/$name.res" != "$reference_file" ]; then
    python3 "$P/cases/compare-outputs.py" "$reference_file" \
      "$out/$name.res" --rtol 1e-8 --groups 0 1-3 4-6 | sed 's/^/   /' |
      tee -a "$out/summary.txt"
  fi
}
