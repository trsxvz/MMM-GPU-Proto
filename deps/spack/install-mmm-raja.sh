#!/bin/bash
# Installation of the environment mmm-raja (RAJA for the raja backend of TFEL
# and MGIS), as deps/spack/install-mmm-gpu.sh: on the cores 2-7 and 10-15
# (the cores 0, 1, 8 and 9 stay free for the tests), with the lowest
# priority, under scripts/memory-governor.sh.
#
# usage: deps/spack/install-mmm-raja.sh [spack]
set -u
P=$(cd "$(dirname "$0")/../.." && pwd)
SPACK=${1:-$HOME/Bureau/spack-mmm/bin/spack}
export SPACK_DISABLE_LOCAL_CONFIG=1
E="$P/deps/spack/mmm-raja"
echo "$(date '+%F %T') install"
taskset -c 2-7,10-15 nice -n 19 "$SPACK" -e "$E" install -j 4 &
pid=$!
"$P/scripts/memory-governor.sh" $pid &
wait $pid
status=$?
echo "$(date '+%F %T') exit status $status"
exit $status
