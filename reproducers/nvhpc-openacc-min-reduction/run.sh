#!/bin/bash
# Minimum reductions of OpenACC on 64 bits integers with nvc++, on the GPU
# (-acc=gpu) and on the cores of the host (-acc=multicore), at -O0 and -O2.
# The binaries are written in scratch/reproducers/nvhpc-openacc-min-reduction.
# usage: run.sh [nvc++]
set -u
here=$(cd "$(dirname "$0")" && pwd)
NVCXX=${1:-nvc++}
out="$here/../../scratch/reproducers/nvhpc-openacc-min-reduction"
mkdir -p "$out"
status=0
for target in gpu multicore; do
  for level in -O0 -O2; do
    exe="$out/min-reduction-$target$level"
    "$NVCXX" -std=c++20 "$level" -acc="$target" "$here/min-reduction.cxx" \
      -o "$exe" || exit 1
    echo "== -acc=$target $level"
    "$exe" || status=1
  done
done
exit $status
