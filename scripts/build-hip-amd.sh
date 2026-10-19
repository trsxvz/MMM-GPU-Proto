#!/usr/bin/env bash
# Build of an MFEM compiled with HIP for an AMD GPU, then of mfem-mgis
# and of the examples with it. The architecture is given by HIP_ARCH:
# - gfx90a (default, MI210 and MI250): compilation and link only, the
#   development machine not having this GPU;
# - gfx1036: the integrated GPU of the development machine (Ryzen), on
#   which the tests really run with GPU_TESTS=ON. /dev/kfd is reserved to
#   the group render: the tests go through `sg render` if the session was
#   opened before the addition to this group.
# The directories of an architecture other than gfx90a have its name as
# suffix (deps/hip-amd/mfem-src-gfx1036, build/mfem-mgis-hip-gfx1036...);
# hypre, compiled for the CPU, is common.
#
# ROCm comes from the Python packages of ROCm (TheRock), installed without
# administrator rights in a virtual environment:
#   python3 -m venv $ROCM_VENV
#   $ROCM_VENV/bin/pip install --index-url \
#     https://stable.repo.amd.com/rocm/whl-next/ \
#     "rocm[libraries,devel,device-gfx90a]==10.1.0"
#   $ROCM_VENV/bin/rocm-sdk init
# and, for another architecture, its package of library kernels:
#   $ROCM_VENV/bin/pip install --index-url \
#     https://stable.repo.amd.com/rocm/whl-next/ \
#     "rocm-sdk-device-gfx1036==10.1.0"
#
# Pitfalls:
# - a system package (libamdhip64-dev of Ubuntu) installs old HIP headers
#   in /usr/include/hip, which the compiler of ROCm reads before its own:
#   the headers of ROCm are given explicitly (-isystem) to hipcc by
#   HIPCC_COMPILE_FLAGS_APPEND and to the HIP compiler of CMake by
#   CMAKE_HIP_FLAGS;
# - scripts/env.sh exports HIP_PLATFORM=nvidia (HIP on NVIDIA GPU for
#   TFEL), which hipcc and hipconfig follow: MFEM is compiled with
#   HIP_PLATFORM=amd, and mfem-mgis chooses the AMD platform by itself when
#   MFEM uses HIP;
# - for gfx1036 only, the libraries of behaviours are compiled at -O1 (see
#   build_environment).
#
# hypre is compiled for the CPU (MFEM with HIP accepts a hypre on the host),
# or with HIP for the GPU if HYPRE_HIP=ON (needed by the tests on the GPU:
# the execution settings refuse a linear solver of hypre that computes on
# the host when it is placed on the GPU), from the archive cached by spack.
# With HYPRE_HIP=ON, MFEM and the rest have the additional suffix
# -hypre-hip. MFEM 4.10 comes from the official repository (tag v4.10).
# MPI, METIS, SuiteSparse and OpenBLAS are those of the spack environment
# given by SPACK_ENVIRONMENT (deps/spack/mmm-gpu by default). Result:
# deps/hip-amd/install/{hypre,mfem}, or HIP_DEPS_DIR/install/{hypre,mfem}.
# MFEM uses the hypre given by HYPRE_DIR, if any, instead of the one built
# by the step hypre (the continuous integration takes the hypre of its spack
# environment).
#
# The step mfem-mgis builds mfem-mgis (directory mfem-mgis of the prototype)
# with this MFEM, and the MGIS and TFEL of the prototype, in
# build/mfem-mgis-hip, with its tests (targets build-tests and
# BehaviourTest), those that need a GPU only with GPU_TESTS=ON, then runs
# them.
#
# The step examples installs this mfem-mgis in install/mfem-mgis-hip and
# builds the examples of mfem-mgis-examples with it (installed
# configuration, behaviours compiled for the AMD GPU), in
# build/examples-hip, then runs its tests (those on GPU with GPU_TESTS=ON,
# option MFEM_MGIS_EXAMPLES_GPU_TESTS).
#
# The compilations and the tests run on the cores given by CPUS (0,1,8,9 by
# default), or without restriction if CPUS is empty.
#
# Usage: [HIP_ARCH=gfx1036] [GPU_TESTS=ON] [HYPRE_HIP=ON] \
#   scripts/build-hip-amd.sh [hypre|mfem|mfem-mgis|examples|all]

set -euo pipefail

P=$(cd "$(dirname "$0")/.." && pwd)
ROCM_VENV=${ROCM_VENV:-$HOME/Bureau/rocm-sdk-therock/venv}
ROCM=$("$ROCM_VENV/bin/rocm-sdk" path --root)
HIP_ARCH=${HIP_ARCH:-gfx90a}
GPU_TESTS=${GPU_TESTS:-OFF}
HYPRE_HIP=${HYPRE_HIP:-OFF}
# suffix of the directories specific to the architecture
if [ "$HIP_ARCH" = gfx90a ]; then
  A=
else
  A=-$HIP_ARCH
fi
# suffixes of hypre and of the directories that depend on it
if [ "$HYPRE_HIP" = ON ]; then
  H=-hip-$HIP_ARCH
  A=$A-hypre-hip
else
  H=
fi
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
JOBS=${JOBS:-4}
SPACK_ENVIRONMENT=${SPACK_ENVIRONMENT:-$P/deps/spack/mmm-gpu}
# cores used by the compilation and by the tests, all if empty
CPUS=${CPUS-0,1,8,9}
on_cpus() {
  if [ -n "$CPUS" ]; then
    taskset -c "$CPUS" "$@"
  else
    "$@"
  fi
}

location() {
  SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$SPACK_ENVIRONMENT" \
    location -i "$1"
}
MPI=$(location mpi)
METIS=$(location metis)
SUITESPARSE=$(location suite-sparse)
OPENBLAS=$(location openblas)

D=${HIP_DEPS_DIR:-$P/deps/hip-amd}
mkdir -p "$D"
HYPRE_DIR=${HYPRE_DIR:-$D/install/hypre$H}

build_hypre() {
  local version=3.2.0
  # archive cached by spack (checksum of the recipe)
  local sha=5273205a310fb6aa3ae506ce216760fb67b30e02024874f3cdb8b811e4801de7
  local spack_root
  spack_root=$(dirname "$(dirname "$SPACK")")
  local archive=$spack_root/var/spack/cache/_source-cache/archive/${sha:0:2}/$sha.tar.gz
  echo "$sha  $archive" | sha256sum -c -
  local src=$D/hypre-src$H b=$D/build-hypre$H
  rm -rf "$src" "$b"
  mkdir -p "$src"
  tar -xf "$archive" -C "$src" --strip-components=1
  local gpu=()
  if [ "$HYPRE_HIP" = ON ]; then
    # hypre compiled by the HIP compiler of ROCm (which is not the MPI
    # wrapper: the headers of MPI are given to it), libraries rocsparse,
    # rocrand and rocblas of ROCm; without Umpire, which hypre enables by
    # default for a GPU (not installed here: hypre then allocates itself)
    export HIP_PATH=$ROCM ROCM_PATH=$ROCM HIP_PLATFORM=amd
    gpu=(-DHYPRE_ENABLE_HIP=ON -DCMAKE_HIP_ARCHITECTURES="$HIP_ARCH"
      -DHYPRE_ENABLE_UMPIRE=OFF
      -DCMAKE_HIP_COMPILER="$ROCM/lib/llvm/bin/clang++"
      -DCMAKE_HIP_FLAGS="-isystem $ROCM/include -isystem $MPI/include"
      -DCMAKE_PREFIX_PATH="$ROCM")
  fi
  cmake -S "$src/src" -B "$b" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$MPI/bin/mpicc" \
    -DCMAKE_CXX_COMPILER="$MPI/bin/mpic++" \
    -DHYPRE_ENABLE_MPI=ON -DBUILD_SHARED_LIBS=ON \
    "${gpu[@]}" \
    -DCMAKE_INSTALL_PREFIX="$D/install/hypre$H"
  on_cpus cmake --build "$b" -j "$JOBS"
  cmake --install "$b"
  echo "hypre $version installed in $D/install/hypre$H"
}

build_mfem() {
  # MFEM 4.10, as the environment mmm-gpu; the sources already extracted
  # are kept (the compilation resumes where it stopped)
  local src=$D/mfem-src$A
  if [ ! -f "$src/makefile" ]; then
    git clone --quiet --depth 1 --branch v4.10 \
      https://github.com/mfem/mfem.git "$src"
  fi
  local blas="-Wl,-rpath,$OPENBLAS/lib -L$OPENBLAS/lib -lopenblas"
  export HIPCC_COMPILE_FLAGS_APPEND="-isystem $ROCM/include"
  # AMD platform, even if scripts/env.sh chooses NVIDIA for TFEL
  export HIP_PATH=$ROCM ROCM_PATH=$ROCM HIP_PLATFORM=amd
  cd "$src"
  # `make config` rewrites config/config.mk, on which all the objects
  # depend: it is only run again for a new extraction of the sources
  [ -f config/config.mk ] || make config \
    PREFIX="$D/install/mfem$A" \
    MFEM_USE_MPI=YES MFEM_USE_HIP=YES MFEM_USE_METIS_5=YES \
    MFEM_USE_SUITESPARSE=YES MFEM_USE_ZLIB=YES \
    MFEM_SHARED=YES MFEM_STATIC=NO \
    HIP_CXX="$ROCM/bin/hipcc" HIP_DIR="$ROCM" HIP_ARCH="$HIP_ARCH" \
    MPICXX="$MPI/bin/mpic++" \
    HYPRE_DIR="$HYPRE_DIR" \
    HYPRE_LIB="-Wl,-rpath,$HYPRE_DIR/lib -L$HYPRE_DIR/lib -lHYPRE" \
    METIS_DIR="$METIS" \
    METIS_LIB="-Wl,-rpath,$METIS/lib -L$METIS/lib -lmetis" \
    SUITESPARSE_DIR="$SUITESPARSE" \
    SUITESPARSE_OPT="-I$SUITESPARSE/include -I$SUITESPARSE/include/suitesparse" \
    LAPACK_LIB="$blas"
  on_cpus make -j "$JOBS" lib
  make install
  echo "MFEM installed in $D/install/mfem$A"
}

# clean environment, without scripts/env.sh: the view of the spack
# environment mmm-gpu contains the MFEM compiled with CUDA, whose headers
# would come before those of the MFEM compiled with HIP (the include
# directory of the view is the one of MPI). Only the TFEL of the prototype
# and MPI are in the paths, no other installation of TFEL or MFEM.
build_environment() {
  tfel=$P/install/tfel-mmm-gpu
  unset CMAKE_PREFIX_PATH PKG_CONFIG_PATH HIP_PLATFORM HIP_COMPILER
  export PATH=$tfel/bin:$MPI/bin:/usr/local/bin:/usr/bin:/bin
  export LD_LIBRARY_PATH=$tfel/lib
  export HIPCC_COMPILE_FLAGS_APPEND="-isystem $ROCM/include"
  # With ROCm 10.1 (packages of TheRock, AMD clang 24.0.0git d6f6cb69), the
  # kernels of finite strain behaviours in 3D compiled for gfx1036 at -O2 or
  # -O3 never end (examples rve and ssna303_3d_hypre, defect of the code
  # generation); they are right at -O0 and -O1
  # (reproducers/hip-finite-strain-3d-hang). In the steps
  # mfem-mgis and examples, hipcc only compiles the libraries of behaviours
  # generated by MFront (the HIP code of mfem-mgis is compiled by the clang
  # of ROCm, MFEM by its own step): -O1, given last, applies to them, for
  # this architecture only. The libraries of an existing build are not
  # compiled again when this option changes (remove their gpu-build
  # directories).
  if [ "$HIP_ARCH" = gfx1036 ]; then
    HIPCC_COMPILE_FLAGS_APPEND="$HIPCC_COMPILE_FLAGS_APPEND -O1"
  fi
}

build_mfem_mgis() {
  build_environment
  local b=$P/build/mfem-mgis-hip$A
  cmake -S "$P/mfem-mgis" -B "$b" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
    -DMPI_C_COMPILER="$MPI/bin/mpicc" -DMPI_CXX_COMPILER="$MPI/bin/mpicxx" \
    -DCMAKE_HIP_FLAGS="-isystem $ROCM/include" \
    -DMFEM_DIR="$D/install/mfem$A" \
    -DMFrontGenericInterface_DIR="$P/install/mgis/share/mgis/cmake" \
    -DTFELHOME="$tfel" -DTFEL_DIR="$tfel/share/tfel/cmake" \
    -Denable-gpu-tests="$GPU_TESTS" \
    -DCMAKE_INSTALL_PREFIX="$P/install/mfem-mgis-hip$A"
  on_cpus cmake --build "$b" -j "$JOBS"
  # the library of behaviours of the tests (BehaviourTest) depends on the
  # target check, not on build-tests
  on_cpus cmake --build "$b" -j "$JOBS" \
    --target build-tests BehaviourTest
  run_tests "$b"
}

# runs the tests of the given build; on the GPU, through `sg render` if
# the session does not have the group render yet (sg, a setuid program,
# removes LD_LIBRARY_PATH from the environment: it is given again)
#
# On the GPU, the kernels and the copies are serialised: with ROCm 10.1
# (HIP 7.16, packages of TheRock) on the integrated GPU gfx1036, the writes
# of a finished kernel to device memory are lost when the GPU stays idle
# about 300 us before the next operation, for instance while the code of a
# shared library is loaded at the first launch of one of its kernels
# (reproducers/hip-idle-write-loss, ROCm/rocm-systems#13469): wrong
# residuals (imposed components not zeroed) and memory faults of hypre.
#
# The OpenMPI of the environment mmm-gpu is built with CUDA: without
# OMPI_MCA_accelerator=null, each process opens a context on the NVIDIA GPU
# of the machine, which may be used by other jobs.
run_tests() {
  local serialize="export OMPI_MCA_accelerator=null; "
  if [ "$GPU_TESTS" = ON ]; then
    serialize="${serialize}export AMD_SERIALIZE_KERNEL=3 AMD_SERIALIZE_COPY=3; "
  fi
  # timeout of 600 s per test: a blocked test does not hold the cores for 1500 s
  local launcher=
  if [ -n "$CPUS" ]; then
    launcher="taskset -c '$CPUS' "
  fi
  local cmd="${serialize}cd '$1' && ${launcher}ctest --timeout 600 --output-on-failure"
  if [ "$GPU_TESTS" = ON ] && ! id -nG | grep -qw render; then
    sg render -c "export LD_LIBRARY_PATH='$LD_LIBRARY_PATH'; $cmd"
  else
    bash -c "$cmd"
  fi
}

build_examples() {
  build_environment
  cmake --install "$P/build/mfem-mgis-hip$A"
  local b=$P/build/examples-hip$A
  cmake -S "$P/mfem-mgis-examples" -B "$b" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
    -DMPI_C_COMPILER="$MPI/bin/mpicc" -DMPI_CXX_COMPILER="$MPI/bin/mpicxx" \
    -DCMAKE_HIP_FLAGS="-isystem $ROCM/include" \
    -DMFEMMGIS_DIR="$P/install/mfem-mgis-hip$A/share/mfem-mgis/cmake" \
    -DMFEM_DIR="$D/install/mfem$A" \
    -DMFEM_MGIS_EXAMPLES_GPU_TESTS="$GPU_TESTS"
  on_cpus cmake --build "$b" -j "$JOBS"
  run_tests "$b"
}

case "${1:-all}" in
  hypre) build_hypre ;;
  mfem) build_mfem ;;
  mfem-mgis) build_mfem_mgis ;;
  examples) build_examples ;;
  all) build_hypre; build_mfem; build_mfem_mgis; build_examples ;;
  *) echo "usage: $0 [hypre|mfem|mfem-mgis|examples|all]" >&2; exit 1 ;;
esac
