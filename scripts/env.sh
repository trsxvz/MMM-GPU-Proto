# Environment of the prototype, to be sourced: `source scripts/env.sh`.
#
# PROTO_ENVIRONMENT selects the spack environment providing the compilers
# and the libraries:
# - mmm-gpu (default): environment of deps/spack/mmm-gpu (Kokkos 5.2.2,
#   AdaptiveCpp 25.10, HIP 7.2.3 and LLVM 20.1.8 built by spack, CUDA
#   13.0.2 of the NVIDIA HPC SDK for nvcc and HIP);
# - tfelgpu_local: view of the spack environment of the shell, with Kokkos
#   and AdaptiveCpp built in install/ and HIP 5.7 (CUDA 12.9).
# TFEL is built in build/<tfel> and installed in install/<tfel>, where
# <tfel> is `tfel` for tfelgpu_local and `tfel-mmm-gpu` for mmm-gpu
# (PROTO_TFEL_BUILD and PROTO_TFEL_INSTALL).
#
# The spack view tfelgpu_local of the shell contains an older TFEL: it must
# not appear in LD_LIBRARY_PATH or before the installed TFEL in PATH.
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PROTO=$P
export PROTO_ENVIRONMENT=${PROTO_ENVIRONMENT:-mmm-gpu}
# NVIDIA HPC SDK: nvc++
export PROTO_NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11
case "$PROTO_ENVIRONMENT" in
  tfelgpu_local)
    # CUDA 13 (nvcc of the spack view, used by the cuda-nvcc configuration)
    export PROTO_SPACK_VIEW=/home/t/Bureau/spack/var/spack/environments/tfelgpu_local/.spack-env/view
    # CUDA 12.9 of the SDK, used by hipcc 5.7, which does not support
    # CUDA 13, and by AdaptiveCpp
    export PROTO_CUDA=$PROTO_NVHPC/cuda/12.9
    export PROTO_HIP_CUDA=$PROTO_CUDA
    # LLVM 20 built by spack (needed by AdaptiveCpp, which provides libomp)
    export PROTO_LLVM=/home/t/Bureau/spack/opt/spack/linux-x86_64_v4/llvm-20.1.8-v3gg3rqravuex5zixx6xglsqudvhcnjw
    export PROTO_KOKKOS=$P/install/kokkos
    export PROTO_ACPP=$P/install/adaptivecpp
    export PROTO_TFEL_BUILD=$P/build/tfel
    export PROTO_TFEL_INSTALL=$P/install/tfel
    ;;
  mmm-gpu)
    export PROTO_SPACK_VIEW=$P/deps/spack/mmm-gpu/.spack-env/view
    # CUDA 13.0.2 of the SDK, with the layout of the standalone toolkit
    # (deps/spack/make-nvhpc-cuda-prefix.sh), used by nvcc and HIP 7.2.3
    export PROTO_CUDA=/home/t/Bureau/spack-mmm-externals/cuda-13.0.2
    export PROTO_HIP_CUDA=$PROTO_CUDA
    export PROTO_LLVM=$PROTO_SPACK_VIEW
    export PROTO_KOKKOS=$PROTO_SPACK_VIEW
    export PROTO_ACPP=$PROTO_SPACK_VIEW
    export PROTO_TFEL_BUILD=$P/build/tfel-mmm-gpu
    export PROTO_TFEL_INSTALL=$P/install/tfel-mmm-gpu
    # a spack environment activated in the shell must not be seen by CMake:
    # its view may contain an older TFEL, whose headers come before the ones
    # of the installed TFEL as soon as a package of the view is found
    # (HDF5 for instance)
    unset SPACK_ENV SPACK_ENV_VIEW NVHPC_CUDA_HOME ACLOCAL_PATH
    export CMAKE_PREFIX_PATH=$PROTO_TFEL_INSTALL:$PROTO_SPACK_VIEW
    export PKG_CONFIG_PATH=$PROTO_SPACK_VIEW/lib/pkgconfig:$PROTO_SPACK_VIEW/share/pkgconfig
    export CUDA_HOME=$PROTO_CUDA
    # AmgX of mfem-mgis: the AmgX patched for CUDA 13 (classical multigrid,
    # recipe amgx of mmm_gpu), installed but not yet in the environment,
    # whose MFEM contains the unpatched AmgX (see deps/spack/README.md)
    export AMGX_DIR=${AMGX_DIR:-/home/t/Bureau/spack-mmm/opt/spack/linux-x86_64_v3/amgx-2.5.0.20260701-xupvzscun6yxw2kmubrm4eogjxtq6uns}
    ;;
  *)
    echo "env.sh: unknown environment '$PROTO_ENVIRONMENT'" >&2
    return 1 2> /dev/null || exit 1
    ;;
esac
# CUDA 12.9 of the SDK (AdaptiveCpp of install/, clang CUDA)
export PROTO_CUDA129=$PROTO_NVHPC/cuda/12.9
# HIP on NVIDIA GPUs
# hipcc uses nvcc on the NVIDIA platform. HIP_COMPILER must not be set:
# MFront uses this variable to select the HIP compiler (instead of hipcc).
export HIP_PLATFORM=nvidia CUDA_PATH=$PROTO_HIP_CUDA
unset HIP_COMPILER
unset LD_LIBRARY_PATH
export PATH=$PROTO_TFEL_INSTALL/bin:$PROTO_ACPP/bin:$PROTO_SPACK_VIEW/bin:$PROTO_CUDA/bin:$PROTO_NVHPC/compilers/bin:/usr/local/bin:/usr/bin:/bin
# the CUDA runtime 12.9 is also used by the programs compiled by clang
# (clang 20 does not support CUDA 13)
export LD_LIBRARY_PATH=$PROTO_TFEL_INSTALL/lib:$PROTO_KOKKOS/lib:$PROTO_LLVM/lib/x86_64-unknown-linux-gnu:$PROTO_LLVM/lib:$PROTO_CUDA/lib64:$PROTO_CUDA129/lib64:$PROTO_NVHPC/compilers/lib
# oneAPI 2026.1 (icpx, SYCL on the CPU): placed last. The OpenCL runtime of
# the CPU requires the TBB of oneAPI (2023.1), which is backward compatible
# with the TBB of the system.
PROTO_ONEAPI=/opt/intel/oneapi
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PROTO_ONEAPI/compiler/2026.1/lib:$PROTO_ONEAPI/umf/1.1/lib:$PROTO_ONEAPI/tbb/2023.1/lib/intel64/gcc4.8
# OpenMP (Kokkos host execution space, AdaptiveCpp CPU device)
export OMP_PROC_BIND=spread OMP_PLACES=threads
