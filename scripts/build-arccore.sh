#!/bin/bash
# Builds and installs Arccore, the core component of the Arcane framework,
# for the arcane configuration of the generic-parallel interfaces of TFEL
# (see build-tfel.sh), without MPI, as shared libraries (required by
# Arcane). The accelerator mode of Arccore is fixed at its build and is
# given by ARCCORE_MODE:
# - cuda (default): C, C++ and CUDA compiled by the clang of the mmm-gpu
#   environment with CUDA 12.9 (clang 20 does not know CUDA 13), for the
#   architecture of the NVIDIA GPU of the machine. Build in build/arccore,
#   installation in install/arccore.
# - hip-amd: HIP on an AMD GPU, C, C++ and HIP compiled by the clang of
#   ROCm (Python packages of ROCm, TheRock, see build-hip-amd.sh), for the
#   architecture given by HIP_ARCH (gfx1036 by default, the integrated GPU
#   of the development machine). The headers of ROCm are given explicitly
#   (-isystem), since the old HIP headers of the system (/usr/include/hip)
#   would be read first. Build in build/arccore-hip-amd-<arch>,
#   installation in install/arccore-hip-amd-<arch>.
# HIP on the NVIDIA platform is not supported: Arccore builds, but the
# kernels of Arccore are only compiled for HIP by the HIP language of clang
# (`__HIP__`, see ARCCORE_COMPILING_HIP in arccore/base/ArccoreGlobal.h),
# which only targets AMD GPUs (nvcc and clang as a CUDA compiler do not
# define `__HIP__`), so the kernels of the behaviours would never run on the
# GPU.
# oneTBB (Thread queues) and GLib (libglib2.0-dev) are the ones of the
# system. The sources are a shallow clone of the arcane-v4.1.16.0 tag of the
# Arcane framework, made in deps/arcane-framework if no source directory is
# given.
# usage: [ARCCORE_MODE=cuda|hip-amd] [HIP_ARCH=gfx1036] \
#   scripts/build-arccore.sh [source directory]
set -e
P=$(cd "$(dirname "$0")/.." && pwd)
V=$P/deps/spack/mmm-gpu/.spack-env/view
CUDA=/home/t/Bureau/spack-mmm-externals/cuda-12.9.1
SRC=${1:-$P/deps/arcane-framework}
MODE=${ARCCORE_MODE:-cuda}
if [ ! -d "$SRC" ]; then
  git clone --depth 1 --branch arcane-v4.1.16.0 \
    https://github.com/arcaneframework/framework.git "$SRC"
fi
# architecture of the NVIDIA GPU of the machine (12.0 gives 120)
nvidia_architecture() {
  local arch
  arch=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '.')
  if [ -z "$arch" ]; then
    echo "unable to read the architecture of the GPU (nvidia-smi)" >&2
    exit 1
  fi
  echo "$arch"
}
# the installed libraries of Arccore find each other (`$ORIGIN`): a
# library or an executable linked with a run path of the new type
# (RUNPATH, given by hipcc) does not give it to the dependencies of Arccore
common=(-DCMAKE_BUILD_TYPE=Release
  -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arccore
  -DARCCORE_ENABLE_MPI=OFF -DBUILD_SHARED_LIBS=TRUE
  -DCMAKE_INSTALL_RPATH='$ORIGIN' -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON)
case "$MODE" in
  cuda)
    suffix=
    settings=()
    mode=(-DARCCORE_ACCELERATOR_MODE=CUDA
      -DCMAKE_C_COMPILER="$V/bin/clang" -DCMAKE_CXX_COMPILER="$V/bin/clang++"
      -DCMAKE_CUDA_COMPILER="$V/bin/clang++"
      -DCUDAToolkit_ROOT="$CUDA" -DCMAKE_CUDA_ARCHITECTURES="$(nvidia_architecture)")
    ;;
  hip-amd)
    ROCM_VENV=${ROCM_VENV:-$HOME/Bureau/rocm-sdk-therock/venv}
    ROCM=$("$ROCM_VENV/bin/rocm-sdk" path --root)
    HIP_ARCH=${HIP_ARCH:-gfx1036}
    suffix=-hip-amd-$HIP_ARCH
    settings=(HIP_PATH="$ROCM" ROCM_PATH="$ROCM" HIP_PLATFORM=amd)
    mode=(-DARCCORE_ACCELERATOR_MODE=HIP
      -DCMAKE_C_COMPILER="$ROCM/lib/llvm/bin/clang"
      -DCMAKE_CXX_COMPILER="$ROCM/lib/llvm/bin/clang++"
      -DCMAKE_HIP_COMPILER="$ROCM/lib/llvm/bin/clang++"
      -DCMAKE_HIP_ARCHITECTURES="$HIP_ARCH" -DCMAKE_PREFIX_PATH="$ROCM"
      -DCMAKE_CXX_FLAGS="-isystem $ROCM/include"
      -DCMAKE_HIP_FLAGS="-isystem $ROCM/include")
    ;;
  *)
    echo "unknown mode '$MODE' (cuda or hip-amd)" >&2
    exit 1
    ;;
esac
BUILD=$P/build/arccore$suffix
env -u LD_LIBRARY_PATH "${settings[@]}" cmake -S "$SRC" -B "$BUILD" \
  -DCMAKE_INSTALL_PREFIX="$P/install/arccore$suffix" "${common[@]}" "${mode[@]}"
env -u LD_LIBRARY_PATH "${settings[@]}" cmake --build "$BUILD" -j"${PROTO_JOBS:-8}"
env -u LD_LIBRARY_PATH "${settings[@]}" cmake --install "$BUILD"
