#!/bin/bash
# Time of a synchronous call of the generic-parallel interface (CUDA
# backend) on a batch of an elastic behaviour (Elasticity of the tests of
# MGIS), with the headers of tfel-dev at a given commit and at the working
# tree. The data are in managed memory.
# usage: reproducers/sync-call-overhead/run.sh [reference commit]
set -e
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
REFERENCE=${1:-1cf42bd07}
source "$P/scripts/env.sh" > /dev/null 2>&1
CUDA=$(dirname "$(dirname "$(which nvcc)")")
export LD_LIBRARY_PATH=$CUDA/lib64:$LD_LIBRARY_PATH
W=$P/scratch/reproducers/sync-call-overhead
rm -rf "$W" && mkdir -p "$W/reference" && cd "$W"
git -C "$P/tfel-dev" archive "$REFERENCE" mfront/include include | tar -x -C reference
for v in reference current; do
  if [ $v = reference ]; then
    I="-I$W/reference/mfront/include -I$W/reference/include"
  else
    I="-I$P/tfel-dev/mfront/include -I$P/tfel-dev/include"
  fi
  mkdir -p $v && cd $v
  cat > c.json << EOC
interfaces_options: {generic-parallel: {backend: "cuda"}},
compilation_options : {
  cuda : {
    compiler: "$(which nvcc)",
    compilation_flags: "$I -O3 -std=c++20 -arch=native -diag-suppress 20012,20208,2810,3189 -Xcompiler -fPIC --expt-relaxed-constexpr"
  },
  cxx : {compilation_flags: "$I -O3 -std=c++20"}
}
EOC
  mfront --obuild --interface=generic-parallel --configuration-file=c.json \
    "$P/mgis/tests/Elasticity.mfront" > build.log 2>&1
  g++-14 -O2 -std=c++20 $I -I"$CUDA/include" "$R/bench.cxx" -o bench -ldl \
    -L"$CUDA/lib64" -lcudart
  cd ..
done
for n in 1000 100000; do
  for v in reference current; do
    echo "$v $(./$v/bench $v/src/libGenericParallelBehaviour-cuda.so $n 2000)"
  done
done
