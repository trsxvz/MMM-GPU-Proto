#!/bin/bash
# Installation of the environment mmm-trilinos (phase 4 of the solver
# libraries, see docs/plan-solver-libraries.md) on a machine with 30 GB of
# memory, as deps/spack/install-mmm-solvers.sh: the packages of mmm-solvers
# are reused, only hwloc, Kokkos, Kokkos Kernels, PaStiX, PETSc, Trilinos
# and MFEM are built, on the cores 2-7 and 10-15 (the cores 0, 1, 8 and 9
# stay free for the tests), with the lowest priority. No other large
# compilation shall run at the same time: the compilation of MFEM needs
# about 22 GB for one file.
#
# usage: deps/spack/install-mmm-trilinos.sh [spack]
#
# Phases, each one under scripts/memory-governor.sh:
# 1. the dependencies of MFEM and Trilinos, two packages at a time, 4 jobs
#    each (Trilinos takes about 1 h 30 with nvcc);
# 2. MFEM, alone, 2 jobs (ptxas takes about 22 GB for
#    fem/lor/lor_batched.cpp with sm_120), resumed with one job where it
#    stopped if it fails.
# Each phase is tried again (up to 3 times) when it fails.
set -u
P=$(cd "$(dirname "$0")/../.." && pwd)
SPACK=${1:-$HOME/Bureau/spack-mmm/bin/spack}
export SPACK_DISABLE_LOCAL_CONFIG=1
E="$P/deps/spack/mmm-trilinos"
run() {
  echo "$(date '+%F %T') $*"
  taskset -c 2-7,10-15 nice -n 19 "$SPACK" -e "$E" "$@" &
  local pid=$!
  "$P/scripts/memory-governor.sh" $pid &
  wait $pid
  local status=$?
  echo "$(date '+%F %T') exit status $status"
  return $status
}
# usage: retry jobs-of-the-first-attempt jobs-of-the-next-ones arguments...
retry() {
  local first=$1 next=$2
  shift 2
  run install --dont-restage -j "$first" "$@" && return 0
  for attempt in 2 3; do
    echo "$(date '+%F %T') attempt $attempt"
    run install --dont-restage -j "$next" "$@" && return 0
  done
  return 1
}
status=0
retry 4 4 -p 2 --only dependencies mfem || status=1
retry 4 4 -p 2 trilinos || status=1
retry 2 1 -p 1 || status=1
echo "$(date '+%F %T') done (status $status)"
exit $status
