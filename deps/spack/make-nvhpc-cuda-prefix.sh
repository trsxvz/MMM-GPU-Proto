#!/bin/bash
# Builds a CUDA prefix with the layout of the standalone CUDA toolkit from a
# CUDA of the NVIDIA HPC SDK, whose mathematical libraries (cuBLAS,
# cuSPARSE, cuRAND, cuFFT, cuSOLVER) are installed in a separate directory
# (`math_libs`). The prefix is made of symbolic links. Build systems which
# do not know the layout of the HPC SDK (for instance the makefiles of MFEM,
# or clang in CUDA mode) find all the headers and libraries in it.
# usage: make-nvhpc-cuda-prefix.sh sdk cuda-version destination
#   sdk: directory of a version of the HPC SDK, for instance
#        /opt/nvidia/hpc_sdk/Linux_x86_64/25.11
#   cuda-version: major and minor versions, for instance 13.0
set -eu
sdk=$1
version=$2
dest=$3
cuda=$sdk/cuda/$version
math=$sdk/math_libs/$version
[ -d "$cuda" ] && [ -d "$math" ] || {
  echo "$cuda or $math not found" >&2
  exit 1
}
# link the entries of the source directories in a directory, merging the
# sub-directories found in several sources
merge() {
  local d=$1
  shift
  mkdir -p "$d"
  local s e n
  for s in "$@"; do
    for e in "$s"/*; do
      n=$d/$(basename "$e")
      if [ -L "$n" ] && [ -d "$n" ] && [ -d "$e" ]; then
        # directory already linked from another source: merge both
        local first
        first=$(readlink -f "$n")
        rm "$n"
        merge "$n" "$first" "$e"
      elif [ -d "$n" ] && [ ! -L "$n" ] && [ -d "$e" ]; then
        merge "$n" "$e"
      elif [ ! -e "$n" ]; then
        ln -s "$(readlink -f "$e")" "$n"
      fi
    done
  done
}
rm -rf "$dest"
mkdir -p "$dest"
for e in "$cuda"/*; do
  case $(basename "$e") in
    bin | include | lib64 | targets) ;;
    *) ln -s "$(readlink -f "$e")" "$dest/$(basename "$e")" ;;
  esac
done
# `bin` is a directory of links, not a link: nvcc finds its headers in
# `bin/../targets`, which must be the one of this prefix (with the headers of
# the mathematical libraries), not the one of the SDK
merge "$dest/bin" "$cuda/bin/"
merge "$dest/include" "$cuda/include/" "$math/include/"
merge "$dest/lib64" "$cuda/lib64/" "$math/lib64/"
mkdir -p "$dest/targets/x86_64-linux"
ln -s ../../include "$dest/targets/x86_64-linux/include"
ln -s ../../lib64 "$dest/targets/x86_64-linux/lib"
