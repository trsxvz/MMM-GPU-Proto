#!/bin/bash
# Installation of the environment mmm-gpu on a machine with 30 GB of memory.
#
# usage: deps/spack/install-mmm-gpu.sh [spack]
# where `spack` is the spack executable (default: ~/Bureau/spack-mmm/bin/spack).
#
# Two phases:
# 1. everything but MFEM, two packages at a time, 4 jobs each;
# 2. MFEM, alone, 2 jobs: each ptxas process of the CUDA kernels of MFEM
#    can take about 10 GB.
# Each phase runs under scripts/memory-governor.sh, which pauses the build
# when the available memory is too low. A failed package does not stop the
# other ones: the log tells which ones failed.
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
status=0
run install -p 2 -j 4 --only dependencies mfem || status=1
run install -p 2 -j 4 kokkos kokkos-kernels cudss || status=1
run install -p 1 -j 2 || status=1
echo "$(date '+%F %T') done (status $status)"
exit $status
