#!/bin/bash
# Builds and installs the TFEL fork of the prototype.
# PROTO_ENVIRONMENT (see scripts/env.sh) selects the environment: mmm-gpu
# by default (deps/spack/mmm-gpu: Kokkos, AdaptiveCpp, HIP 7.2.3 and LLVM
# from spack, CUDA 13.0.2 from the HPC SDK for nvcc and HIP; TFEL built in
# build/tfel-mmm-gpu and installed in install/tfel-mmm-gpu) or
# tfelgpu_local (old environment: Kokkos and AdaptiveCpp from install/,
# hipcc 5.7 with the CUDA 12.9 of the HPC SDK; build/tfel and install/tfel).
# generic-parallel configurations: cuda-nvcc, hip-hipcc (nvidia
# platform), stlpar-parunseq-gcc, stlpar-parunseq-nvhpc-gpu (nvc++ from the
# HPC SDK), kokkos (CUDA + OpenMP), raja (CUDA + OpenMP, RAJA of the
# environment mmm-raja, mmm-gpu only, if installed), sycl-default-acpp and
# sycl-gpu-acpp (AdaptiveCpp), openmp-host-gcc, openmp-offload-nvhpc-gpu
# (nvc++ -mp=gpu), openacc-nvhpc-gpu and openacc-nvhpc-multicore (nvc++ -acc=gpu
# and -acc=multicore), stlpar-parunseq-clang (clang 20 from the LLVM of spack:
# its library is in LD_LIBRARY_PATH, so the clang of the system cannot be used),
# cuda-clang and hip-clang (same clang, with CUDA 12.9: clang 20 does not
# know CUDA 13), stlpar-parunseq-icpx and sycl-default-icpx (oneAPI
# 2026.1, SYCL on CPU), arcane (accelerator API of Arccore built with CUDA
# by scripts/build-arccore.sh, mmm-gpu only, if installed: same clang and
# CUDA 12.9 as cuda-clang); TDLS enabled. To change the configurations of an
# existing build directory:
# cmake -Dgeneric-parallel-configurations="..." -U "tfel-substitutions*" .
# (with -DRAJA_DIR=<installation of RAJA>/lib/cmake/raja when adding raja,
# -DArccore_ROOT=<installation of Arccore> when adding arcane).
# The tests run in the environment of scripts/env.sh.
set -e
P=$(cd "$(dirname "$0")" && pwd)
case "${PROTO_ENVIRONMENT:-mmm-gpu}" in
  tfelgpu_local)
    V=/home/t/Bureau/spack/var/spack/environments/tfelgpu_local/.spack-env/view
    CUDA=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/12.9
    CUDA_TOOLKIT=$V
    KOKKOS=$P/install/kokkos
    ACPP=$P/install/adaptivecpp/bin/acpp
    CLANGXX=/home/t/Bureau/spack/opt/spack/linux-x86_64_v4/llvm-20.1.8-v3gg3rqravuex5zixx6xglsqudvhcnjw/bin/clang++
    BUILD=$P/build/tfel
    INSTALL=$P/install/tfel
    ;;
  mmm-gpu)
    V=$P/deps/spack/mmm-gpu/.spack-env/view
    CUDA=/home/t/Bureau/spack-mmm-externals/cuda-13.0.2
    CUDA_TOOLKIT=$CUDA
    KOKKOS=$V
    # RAJA and its camp built with OpenMP (deps/spack/mmm-raja), if
    # installed by deps/spack/install-mmm-raja.sh
    RAJA=$P/deps/spack/mmm-raja/.spack-env/view
    if [ ! -d "$RAJA/lib/cmake/raja" ]; then
      echo "no raja configuration: RAJA is not installed" \
        "(deps/spack/install-mmm-raja.sh)"
      RAJA=
    fi
    # Arccore built with CUDA, if installed by scripts/build-arccore.sh
    ARCCORE=$P/install/arccore
    if [ ! -f "$ARCCORE/include/arccore/arccore_config.h" ]; then
      echo "no arcane configuration: Arccore is not installed" \
        "(scripts/build-arccore.sh)"
      ARCCORE=
    fi
    ACPP=$V/bin/acpp
    CLANGXX=$V/bin/clang++
    BUILD=$P/build/tfel-mmm-gpu
    INSTALL=$P/install/tfel-mmm-gpu
    # spack environment possibly activated in the shell: made invisible
    # to CMake (see scripts/env.sh)
    unset SPACK_ENV SPACK_ENV_VIEW NVHPC_CUDA_HOME ACLOCAL_PATH
    export CMAKE_PREFIX_PATH=$V
    export PKG_CONFIG_PATH=$V/lib/pkgconfig:$V/share/pkgconfig
    export CUDA_HOME=$CUDA
    ;;
  *)
    echo "unknown environment '$PROTO_ENVIRONMENT'" >&2
    exit 1
    ;;
esac
export PATH=$V/bin:$CUDA/bin:$PATH
export HIP_PLATFORM=nvidia
export CUDA_PATH=$CUDA
CONFIGURATIONS="cuda-nvcc;cuda-clang;hip-hipcc;hip-clang"
CONFIGURATIONS="$CONFIGURATIONS;stlpar-parunseq-gcc;stlpar-parunseq-nvhpc-gpu"
CONFIGURATIONS="$CONFIGURATIONS;stlpar-parunseq-clang;stlpar-parunseq-icpx;kokkos"
CONFIGURATIONS="$CONFIGURATIONS;sycl-default-acpp;sycl-gpu-acpp;sycl-default-icpx"
CONFIGURATIONS="$CONFIGURATIONS;openmp-host-gcc;openmp-offload-nvhpc-gpu"
CONFIGURATIONS="$CONFIGURATIONS;openacc-nvhpc-gpu;openacc-nvhpc-multicore"
RAJA_OPTIONS=()
if [ -n "${RAJA:-}" ]; then
  CONFIGURATIONS="$CONFIGURATIONS;raja"
  # installation of RAJA itself, not the view: its include directory and the
  # run path written in the generated libraries only hold RAJA and camp
  RAJA_CONFIG=$(readlink -f "$RAJA/lib/cmake/raja/raja-config.cmake")
  RAJA_OPTIONS=(-DRAJA_DIR="$(dirname "$RAJA_CONFIG")")
fi
ARCANE_OPTIONS=()
if [ -n "${ARCCORE:-}" ]; then
  CONFIGURATIONS="$CONFIGURATIONS;arcane"
  ARCANE_OPTIONS=(-DArccore_ROOT="$ARCCORE")
fi
mkdir -p "$BUILD"
cd "$BUILD"
if [ ! -f CMakeCache.txt ]; then
  cmake "$P/tfel" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$INSTALL" \
    -DCMAKE_CXX_COMPILER=g++-14 -DCMAKE_C_COMPILER=gcc-14 \
    -Denable-website=OFF -Denable-python=OFF -Denable-python-bindings=OFF \
    -Denable-fortran=OFF -Denable-testing=ON \
    -Denable-tdls=ON -Dtdls_DIR=/home/t/Documents/GitHub/tdls/build-install/lib/cmake/tdls \
    -Dgeneric-parallel-configurations="$CONFIGURATIONS" \
    -DKokkos_ROOT="$KOKKOS" "${RAJA_OPTIONS[@]}" "${ARCANE_OPTIONS[@]}" \
    -DTFEL_NVHPC_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin/nvc++ \
    -DTFEL_ACPP_COMPILER="$ACPP" \
    -DTFEL_CLANGXX_COMPILER="$CLANGXX" \
    -DTFEL_ICPX_COMPILER=/opt/intel/oneapi/compiler/2026.1/bin/icpx \
    -DTFEL_CLANG_CUDA_PATH=/home/t/Bureau/spack-mmm-externals/cuda-12.9.1 \
    -DCUDAToolkit_ROOT=$CUDA_TOOLKIT
fi
make -j16 ${@:-install}
