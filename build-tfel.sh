#!/bin/bash
# Compile et installe le fork de TFEL du prototype dans install/tfel.
# Configurations generic-parallel : cuda-nvcc, hip-hipcc (plateforme nvidia,
# CUDA 12.9 du HPC SDK : hipcc 5.7 ne gère pas CUDA 13), stlpar-parunseq-gcc,
# stlpar-parunseq-nvhpc-gpu (nvc++ du HPC SDK), kokkos (Kokkos partagé de
# install/kokkos, CUDA + OpenMP), sycl-default-acpp (AdaptiveCpp de
# install/adaptivecpp), sycl-gpu-acpp, openmp-host-gcc, openmp-offload-nvhpc-gpu
# (nvc++ -mp=gpu), stlpar-parunseq-clang (clang 20 du LLVM de spack : sa
# bibliothèque est dans LD_LIBRARY_PATH, le clang du système ne peut donc
# pas servir), cuda-clang et hip-clang (même clang, avec le CUDA 12.9 du
# HPC SDK : clang 20 ne connaît pas CUDA 13), stlpar-parunseq-icpx et
# sycl-default-icpx (oneAPI 2026.1, SYCL sur CPU) ; TDLS activé. Pour changer
# les configurations d'un dossier de construction existant :
# cmake -Dgeneric-parallel-configurations="..." . Les tests s'exécutent dans
# l'environnement de scripts/env.sh.
# PROTO_ENVIRONMENT=mmm-gpu (voir scripts/env.sh) utilise l'environnement
# spack de deps/spack/mmm-gpu (Kokkos, AdaptiveCpp, HIP 7.2.3 et LLVM de
# spack, CUDA 13.0.2 du HPC SDK pour nvcc et HIP) et construit TFEL dans
# build/tfel-mmm-gpu, installé dans install/tfel-mmm-gpu.
set -e
P=$(cd "$(dirname "$0")" && pwd)
case "${PROTO_ENVIRONMENT:-tfelgpu_local}" in
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
    ACPP=$V/bin/acpp
    CLANGXX=$V/bin/clang++
    BUILD=$P/build/tfel-mmm-gpu
    INSTALL=$P/install/tfel-mmm-gpu
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
    -DKokkos_ROOT="$KOKKOS" \
    -DTFEL_NVHPC_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin/nvc++ \
    -DTFEL_ACPP_COMPILER="$ACPP" \
    -DTFEL_CLANGXX_COMPILER="$CLANGXX" \
    -DTFEL_ICPX_COMPILER=/opt/intel/oneapi/compiler/2026.1/bin/icpx \
    -DTFEL_CLANG_CUDA_PATH=/home/t/Bureau/spack-mmm-externals/cuda-12.9.1 \
    -DCUDAToolkit_ROOT=$CUDA_TOOLKIT
fi
make -j16 ${@:-install}
