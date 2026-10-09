#!/bin/bash
# Installation of the environment mmm-gpu on a machine with 30 GB of memory.
#
# usage: deps/spack/install-mmm-gpu.sh [spack]
# where `spack` is the spack executable (default: ~/Bureau/spack-mmm/bin/spack).
#
# Phases, each one under scripts/memory-governor.sh, which pauses all the
# build but its largest process when the available memory is too low:
# 1. everything but MFEM, two packages at a time, 4 jobs each (the
#    dependencies of MFEM, then the other roots);
# 2. MFEM, alone, 2 jobs. ptxas takes about 22 GB for the CUDA kernels of
#    fem/lor/lor_batched.cpp (sm_120): if the build fails (lack of memory),
#    it is resumed with one job, where it stopped (--dont-restage, and the
#    package mfem of mmm_gpu keeps the objects of MFEM up to date).
# Each phase is tried again (up to 3 times) when it fails.
set -u
P=$(cd "$(dirname "$0")/../.." && pwd)
SPACK=${1:-$HOME/Bureau/spack-mmm/bin/spack}
export SPACK_DISABLE_LOCAL_CONFIG=1
E="$P/deps/spack/mmm-gpu"
run() {
  echo "$(date '+%F %T') $*"
  "$SPACK" -e "$E" "$@" &
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
retry 4 4 -p 2 kokkos kokkos-kernels cudss || status=1
retry 2 1 -p 1 || status=1
echo "$(date '+%F %T') done (status $status)"
exit $status
