#!/bin/bash
# Local stores of a kernel whose parameter, passed by value, is used through
# its address, without and with __grid_constant__ (see grid-constant.cu).
# usage: reproducers/grid-constant/run.sh
set -e
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
source "$P/scripts/env.sh" > /dev/null 2>&1
W=$P/scratch/reproducers/grid-constant
rm -rf "$W" && mkdir -p "$W" && cd "$W"
for v in copy grid-constant; do
  flag=""
  [ $v = grid-constant ] && flag="-DUSE_GRID_CONSTANT"
  nvcc -O3 -std=c++20 -arch=native $flag "$R/grid-constant.cu" -o $v
  echo "$v: $(cuobjdump -sass $v | grep -c 'STL') local stores in the" \
    "kernel, results $(./$v)"
done
