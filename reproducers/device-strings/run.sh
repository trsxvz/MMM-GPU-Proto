#!/bin/bash
# Compile and run the string constant reproducer with the GPU compilers of
# the prototype, for four variants of the constant used in device code:
# A: static constexpr std::string_view member, B: local constexpr
# std::string_view, C and D: same with a device compatible string type
# (constant-string.hxx). The results are 1 7 13 19 when the code is right.
# usage: reproducers/device-strings/run.sh
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
source "$P/scripts/env.sh" > /dev/null 2>&1
CLANG_CUDA=/home/t/Bureau/spack-mmm-externals/cuda-12.9.1
W=$P/scratch/reproducers/device-strings
rm -rf "$W" && mkdir -p "$W" && cd "$W" || exit 1
run() {
  local name=$1
  shift
  if "$@" > log 2>&1; then
    printf "%-14s %s exit=%s\n" "$name" "$(./t 2>&1 | tail -1)" \
      "$(./t > /dev/null 2>&1; echo $?)"
  else
    printf "%-14s COMPILE-ERROR %s\n" "$name" \
      "$(grep -m1 -iE 'error' log | cut -c1-90)"
  fi
  rm -f t
}
for v in A B C D; do
  case $v in
    A) d="-DMODEL_STATIC_MEMBER" ;;
    B) d="" ;;
    C) d="-DMODEL_STATIC_MEMBER -DMODEL_STRING=ConstantString" ;;
    D) d="-DMODEL_STRING=ConstantString" ;;
  esac
  echo "== variant $v"
  run nvcc nvcc $d -std=c++20 -arch=native --expt-relaxed-constexpr \
    "$R/test.cu" -o t
  run clang-cuda "$PROTO_LLVM/bin/clang++" $d -std=c++20 -x cuda \
    --cuda-path=$CLANG_CUDA --cuda-gpu-arch=sm_120 "$R/test.cu" -o t \
    -L$CLANG_CUDA/lib64 -lcudart
  run hipcc hipcc $d -x cu -arch=native --expt-relaxed-constexpr -std=c++20 \
    "$R/test.cu" -o t
  run nvc++-stdpar nvc++ $d -std=c++20 -O2 -stdpar=gpu -gpu=mem:managed \
    "$R/test-stdpar.cxx" -o t
  run nvc++-openmp nvc++ $d -std=c++20 -O2 -mp=gpu -gpu=mem:managed \
    "$R/test-omp.cxx" -o t
  run acpp acpp $d -std=c++20 -O2 "$R/test-sycl.cxx" -o t
done
