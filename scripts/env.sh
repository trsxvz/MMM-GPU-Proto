# Environment of the prototype, to be sourced: `source scripts/env.sh`.
# The spack view tfelgpu_local of the shell contains an older TFEL: it must
# not appear in LD_LIBRARY_PATH or before install/tfel in PATH.
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PROTO=$P
# CUDA 13 (nvcc of the spack view, used by the cuda-nvcc configuration)
export PROTO_SPACK_VIEW=/home/t/Bureau/spack/var/spack/environments/tfelgpu_local/.spack-env/view
# NVIDIA HPC SDK: nvc++ and CUDA 12.9 (used by hipcc 5.7, which does not
# support CUDA 13, and by AdaptiveCpp)
export PROTO_NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11
export PROTO_CUDA129=$PROTO_NVHPC/cuda/12.9
# LLVM 20 built by spack (needed by AdaptiveCpp, which provides libomp)
export PROTO_LLVM=/home/t/Bureau/spack/opt/spack/linux-x86_64_v4/llvm-20.1.8-v3gg3rqravuex5zixx6xglsqudvhcnjw
# HIP on NVIDIA GPUs
# hipcc uses nvcc on the NVIDIA platform. HIP_COMPILER must not be set:
# MFront uses this variable to select the HIP compiler (instead of hipcc).
export HIP_PLATFORM=nvidia CUDA_PATH=$PROTO_CUDA129
unset HIP_COMPILER
unset LD_LIBRARY_PATH
export PATH=$P/install/tfel/bin:$P/install/adaptivecpp/bin:$PROTO_SPACK_VIEW/bin:$PROTO_NVHPC/compilers/bin:/usr/local/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=$P/install/tfel/lib:$P/install/kokkos/lib:$PROTO_LLVM/lib/x86_64-unknown-linux-gnu:$PROTO_LLVM/lib:$PROTO_CUDA129/lib64:$PROTO_NVHPC/compilers/lib
# OpenMP (Kokkos host execution space, AdaptiveCpp CPU device)
export OMP_PROC_BIND=spread OMP_PLACES=threads
