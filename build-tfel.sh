#!/bin/bash
# Compile et installe le fork de TFEL du prototype dans install/tfel.
# Configurations generic-parallel : cuda-nvcc, hip-hipcc (plateforme nvidia,
# CUDA 12.9 du HPC SDK : hipcc 5.7 ne gère pas CUDA 13), stlpar-parunseq-gcc,
# stlpar-parunseq-nvhpc-gpu (nvc++ du HPC SDK), kokkos (Kokkos partagé de
# install/kokkos, CUDA + OpenMP), sycl-default-acpp (AdaptiveCpp de
# install/adaptivecpp), sycl-gpu-acpp, openmp-host-gcc, openmp-offload-nvhpc-gpu
# (nvc++ -mp=gpu), stlpar-parunseq-clang (clang 20), stlpar-parunseq-icpx et
# sycl-default-icpx (oneAPI 2026.1, SYCL sur CPU) ; TDLS activé. Pour changer
# les configurations d'un dossier de construction existant :
# cmake -Dgeneric-parallel-configurations="..." . Les tests s'exécutent dans
# l'environnement de scripts/env.sh.
set -e
P=$(cd "$(dirname "$0")" && pwd)
V=/home/t/Bureau/spack/var/spack/environments/tfelgpu_local/.spack-env/view
export PATH=$V/bin:$PATH
export HIP_PLATFORM=nvidia
export CUDA_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/12.9
CONFIGURATIONS="cuda-nvcc;hip-hipcc;stlpar-parunseq-gcc;stlpar-parunseq-nvhpc-gpu"
CONFIGURATIONS="$CONFIGURATIONS;stlpar-parunseq-clang;stlpar-parunseq-icpx;kokkos"
CONFIGURATIONS="$CONFIGURATIONS;sycl-default-acpp;sycl-gpu-acpp;sycl-default-icpx"
CONFIGURATIONS="$CONFIGURATIONS;openmp-host-gcc;openmp-offload-nvhpc-gpu"
mkdir -p "$P/build/tfel"
cd "$P/build/tfel"
if [ ! -f CMakeCache.txt ]; then
  cmake "$P/tfel" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$P/install/tfel" \
    -DCMAKE_CXX_COMPILER=g++-14 -DCMAKE_C_COMPILER=gcc-14 \
    -Denable-website=OFF -Denable-python=OFF -Denable-python-bindings=OFF \
    -Denable-fortran=OFF -Denable-testing=ON \
    -Denable-tdls=ON -Dtdls_DIR=/home/t/Documents/GitHub/tdls/build-install/lib/cmake/tdls \
    -Dgeneric-parallel-configurations="$CONFIGURATIONS" \
    -DKokkos_ROOT="$P/install/kokkos" \
    -DTFEL_NVHPC_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin/nvc++ \
    -DTFEL_ACPP_COMPILER="$P/install/adaptivecpp/bin/acpp" \
    -DTFEL_CLANGXX_COMPILER=/usr/bin/clang++-20 \
    -DTFEL_ICPX_COMPILER=/opt/intel/oneapi/compiler/2026.1/bin/icpx \
    -DCUDAToolkit_ROOT=$V
fi
make -j16 ${@:-install}
