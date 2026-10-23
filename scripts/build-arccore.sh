#!/bin/bash
# Builds and installs Arccore, the core component of the Arcane framework,
# for the arcane configuration of the generic-parallel interfaces of TFEL
# (see build-tfel.sh): accelerator mode CUDA, C, C++ and CUDA compiled by the
# clang of the mmm-gpu environment with CUDA 12.9 (clang 20 does not know
# CUDA 13), for the architecture of the GPU of the machine, without MPI, as
# shared libraries (required by Arcane). oneTBB (Thread queues) and GLib
# (libglib2.0-dev) are the ones of the system. The sources are a shallow
# clone of the arcane-v4.1.16.0 tag of the Arcane framework, made in
# deps/arcane-framework if no source directory is given.
# Build in build/arccore, installation in install/arccore.
# usage: scripts/build-arccore.sh [source directory]
set -e
P=$(cd "$(dirname "$0")/.." && pwd)
V=$P/deps/spack/mmm-gpu/.spack-env/view
CUDA=/home/t/Bureau/spack-mmm-externals/cuda-12.9.1
SRC=${1:-$P/deps/arcane-framework}
if [ ! -d "$SRC" ]; then
  git clone --depth 1 --branch arcane-v4.1.16.0 \
    https://github.com/arcaneframework/framework.git "$SRC"
fi
# architecture of the GPU of the machine (12.0 gives 120)
ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '.')
if [ -z "$ARCH" ]; then
  echo "unable to read the architecture of the GPU (nvidia-smi)" >&2
  exit 1
fi
env -u LD_LIBRARY_PATH cmake -S "$SRC" -B "$P/build/arccore" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$P/install/arccore" \
  -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arccore \
  -DARCCORE_ACCELERATOR_MODE=CUDA -DARCCORE_ENABLE_MPI=OFF \
  -DBUILD_SHARED_LIBS=TRUE \
  -DCMAKE_C_COMPILER="$V/bin/clang" -DCMAKE_CXX_COMPILER="$V/bin/clang++" \
  -DCMAKE_CUDA_COMPILER="$V/bin/clang++" \
  -DCUDAToolkit_ROOT="$CUDA" -DCMAKE_CUDA_ARCHITECTURES="$ARCH"
env -u LD_LIBRARY_PATH cmake --build "$P/build/arccore" -j8
env -u LD_LIBRARY_PATH cmake --install "$P/build/arccore"
