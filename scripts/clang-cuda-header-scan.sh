#!/bin/bash
# Compile each header of TFEL/Math and TFEL/Material alone with clang in
# CUDA mode (host and device passes) and report the inconsistencies of the
# host and device annotations that nvcc accepts but clang refuses:
# definitions without `__host__ __device__` whose declaration has them (or
# the opposite), C++11 attributes placed after `__host__ __device__`, and
# headers using `TFEL_HOST_DEVICE` without including the TFEL configuration.
# Errors due to headers which are not self-contained are not reported.
# CLANG defaults to the clang++ of the former spack environment
# tfelgpu_local (its PROTO_LLVM in scripts/env.sh), not to the LLVM of the
# default environment mmm-gpu: give CLANG=$PROTO_LLVM/bin/clang++ to use
# the latter.
# usage: scripts/clang-cuda-header-scan.sh [number of jobs]
set -u
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
JOBS=${1:-8}
INCLUDES=$PROTO/tfel/include
OUT=$PROTO/scratch/clang-cuda-header-scan
CLANG=${CLANG:-/home/t/Bureau/spack/opt/spack/linux-x86_64_v4/llvm-20.1.8-v3gg3rqravuex5zixx6xglsqudvhcnjw/bin/clang++}
CUDA=${CUDA:-/opt/nvidia/hpc_sdk/Linux_x86_64/25.11}
FLAGS="-x cuda --cuda-gpu-arch=sm_120 --cuda-path=$CUDA/cuda/12.9"
FLAGS="$FLAGS -I$CUDA/math_libs/12.9/include -Wno-unknown-cuda-version"
FLAGS="$FLAGS -std=c++20 -fsyntax-only -I$INCLUDES -DTFEL_ARCH64"
rm -rf "$OUT" && mkdir -p "$OUT"
scan() {
  h=$1
  d=$(mktemp -d -p "$OUT")
  echo "#include \"$h\"" > "$d/t.cu"
  $CLANG $FLAGS "$d/t.cu" 2>&1 | grep "error:" |
    grep -E "cannot overload|attribute list cannot appear|unknown type name 'TFEL_HOST_DEVICE'" |
    sed "s|$INCLUDES/||"
  rm -rf "$d"
}
export -f scan
export OUT CLANG FLAGS INCLUDES
(cd "$INCLUDES" && find TFEL/Math TFEL/Material -name "*.hxx" | sort) |
  xargs -P "$JOBS" -I{} bash -c 'scan {}' | sort -u > "$OUT/report.txt"
echo "$(wc -l < "$OUT/report.txt") inconsistencies (see $OUT/report.txt)"
[ ! -s "$OUT/report.txt" ]
