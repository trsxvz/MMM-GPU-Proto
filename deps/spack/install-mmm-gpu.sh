#!/bin/bash
# Installation of the environment mmm-gpu on a machine with 30 GB of memory.
#
# usage: deps/spack/install-mmm-gpu.sh [spack]
# where `spack` is the spack executable (default: ~/Bureau/spack-mmm/bin/spack).
#
# Two phases:
# 1. everything but MFEM, two packages at a time, 4 jobs each (in two
#    steps: the dependencies of MFEM, then the other roots), under
#    scripts/memory-governor.sh, which pauses the build when the available
#    memory is too low;
# 2. MFEM, alone, one job, without the governor: ptxas takes about 22 GB
#    for the CUDA kernels of fem/lor/lor_batched.cpp (MFEM 4.9, sm_120).
#    Pausing it does not free its memory: the governor would wait forever
#    for the memory to come back.
# A failed package does not stop the other ones: the log tells which ones
# failed. A partial build is resumed where it stopped (--dont-restage).
set -u
P=$(cd "$(dirname "$0")/../.." && pwd)
SPACK=${1:-$HOME/Bureau/spack-mmm/bin/spack}
export SPACK_DISABLE_LOCAL_CONFIG=1
E="$P/deps/spack/mmm-gpu"
# usage: run governor|alone spack-arguments...
run() {
  local mode=$1
  shift
  echo "$(date '+%F %T') $*"
  "$SPACK" -e "$E" "$@" &
  local pid=$!
  if [ "$mode" = governor ]; then
    "$P/scripts/memory-governor.sh" $pid &
  fi
  wait $pid
  local status=$?
  echo "$(date '+%F %T') exit status $status"
  return $status
}
status=0
run governor install --dont-restage -p 2 -j 4 --only dependencies mfem || status=1
run governor install --dont-restage -p 2 -j 4 kokkos kokkos-kernels cudss || status=1
run alone install --dont-restage -p 1 -j 1 || status=1
echo "$(date '+%F %T') done (status $status)"
exit $status
