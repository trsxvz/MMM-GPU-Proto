#!/bin/bash
# Build the generic library of a test behaviour with the MFront of the
# development tree and check the rotation of its tangent operator blocks
# (check-rotations.cxx).
# usage: reproducers/rotations/run.sh [behaviour.mfront] [hypotheses...]
set -e
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
D=$P/tfel-dev
f=${1:-$D/test-suite/behaviours/generic-parallel/OrthotropicTangentOperatorBlocksTest.mfront}
shift || true
hypotheses=${*:-"Tridimensional PlaneStrain Axisymmetrical GeneralisedPlaneStrain AxisymmetricalGeneralisedPlaneStrain"}
source "$P/scripts/env.sh" > /dev/null 2>&1
W=$P/scratch/reproducers/rotations
rm -rf "$W" && mkdir -p "$W" && cd "$W"
export LD_LIBRARY_PATH=$(find "$P/scratch/build-dev2" -name "*.so" -printf "%h\n" | sort -u | paste -sd:):$LD_LIBRARY_PATH
INCLUDES="-I$D/mfront/include -I$D/include" "$P/scratch/build-dev2/mfront/src/mfront" \
  --obuild --interface=generic "$f" > build.log 2>&1
g++-14 -std=c++20 -O2 -Wall "$R/check-rotations.cxx" -o check-rotations -ldl
./check-rotations src/libBehaviour.so "$(basename "$f" .mfront)" $hypotheses
