#!/bin/bash
# Installation of the environment mmm-trilinos-float (single precision twin
# of mmm-trilinos-teko, see its spack.yaml): only the packages built in
# single precision are new (Trilinos, PETSc, hypre, MFEM), the other ones
# are reused. Built with 4 jobs, the lowest priority, under
# scripts/memory-governor.sh, on the cores given by the variable CORES (by
# default all of them), and tried again with 2 jobs if it fails.
#
# usage: deps/spack/install-mmm-trilinos-float.sh [spack]
set -u
P=$(cd "$(dirname "$0")/../.." && pwd)
SPACK=${1:-$HOME/Bureau/spack-mmm/bin/spack}
export SPACK_DISABLE_LOCAL_CONFIG=1
E="$P/deps/spack/mmm-trilinos-float"
CORES=${CORES:-0-$(($(nproc) - 1))}
run() {
  echo "$(date '+%F %T') install -j $1"
  taskset -c "$CORES" nice -n 19 "$SPACK" -e "$E" install --dont-restage -j "$1" &
  local pid=$!
  "$P/scripts/memory-governor.sh" $pid &
  wait $pid
  local status=$?
  echo "$(date '+%F %T') exit status $status"
  return $status
}
run 4 || run 2 || run 2
