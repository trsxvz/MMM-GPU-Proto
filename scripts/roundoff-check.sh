#!/bin/bash
# Sensitivity of a law to round-off errors, to check a ROUNDOFF status of
# scripts/runtime-sweep.py: the OpenMP library of the batch is compiled on
# the host with FMA contraction (-ffp-contract=fast, as nvcc in device code)
# and the reference without (-ffp-contract=off). Both are then compared with
# the options of the sweep. Differences of the same order as on the GPU
# show that the law amplifies round-off errors. Without FMA contraction, the
# differences are null.
#
# usage: scripts/roundoff-check.sh directory function options...
#
# The directory is the one of the law in the results of the sweep of the
# openmp-host-gcc configuration (for instance
# scratch/runtime-openmp-host-gcc/behaviours_AgeingBurger). The function and
# the options are the ones written on the first line of the compare.log file
# of the law in the results of the GPU configuration.
set -eu
P=$(cd "$(dirname "$0")/.." && pwd)
source=$(cd "$1" && pwd)
function=$2
shift 2
source "$P/scripts/env.sh" > /dev/null 2>&1
D=$P/tfel-dev
W=$P/scratch/roundoff-check/$(basename "$source")
rm -rf "$W"
mkdir -p "$(dirname "$W")"
cp -r "$source" "$W"
cd "$W/src"
flags="-I$D/mfront/include -I$D/include -O3 -fopenmp -std=c++20 -march=native
  -fPIC -I../include -I$D/mfront/tests/behaviours -I$D/mfront/tests/properties
  $(tfel-config --cppflags) -I$(tfel-config --include-path)"
libs="-L$(tfel-config --library-path) -fopenmp"
for l in $(tfel-config --library-dependency --material --mfront-profiling); do
  libs="$libs -l$l"
done
generic=$(ls ./*-generic.cxx)
batch=$(ls ./*-generic-parallel-openmp.cxx)
others=$(ls ./*.cxx | grep -v -- "-generic")
for c in off fast; do
  objects=""
  for s in $others; do
    g++-14 $flags -ffp-contract=$c "$s" -c -o "${s%.cxx}-$c.o" &
    objects="$objects ${s%.cxx}-$c.o"
  done
  g++-14 $flags -ffp-contract=$c "$generic" -c -o generic-$c.o &
  g++-14 $flags -ffp-contract=$c "$batch" -c -o batch-$c.o &
  wait
  g++-14 -shared generic-$c.o $objects -o libReference-$c.so $libs
  g++-14 -shared batch-$c.o $objects -o libBatch-$c.so $libs
done
cd "$W"
for c in off fast; do
  echo "batch compiled with -ffp-contract=$c:"
  MFRONT_GPB_COMPARE_NUMBER_OF_POINTS=1009 \
    LD_LIBRARY_PATH=src:${LD_LIBRARY_PATH:-} \
    "$P/scratch/sweep-env/openmp-host-gcc/compare" src/libReference-off.so \
    src/libBatch-$c.so "$function" "$@" | grep "by point" || true
done
