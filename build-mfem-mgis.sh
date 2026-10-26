#!/bin/bash
# Builds, tests and installs mfem-mgis, then builds and tests
# mfem-mgis-examples, mm-opera-hpc and the cases of cases/ (repositories
# `mfem-mgis/`, `mfem-mgis-examples/` and `mm-opera-hpc/`, branch
# `gpu-proto`), against MFEM 4.10 of the spack environment of the machine
# (PROTO_SPACK_ENVIRONMENT, see scripts/env.sh: `mmm-gpu` on the
# development machine), the MGIS installed by build-mgis.sh and the TFEL
# installed by build-tfel.sh.
#
# The other builds do not go through this script: MFEM built with HIP for
# an AMD GPU (scripts/build-hip-amd.sh), and the builds that validate the
# option enable-kokkos (build/mfem-mgis-kokkos-off and
# build/mfem-mgis-kokkos, see docs/validation.md).
#
# usage: ./build-mfem-mgis.sh [mfem-mgis|examples|opera|cases ...]
# Without argument: the four, in this order. mfem-mgis is installed in
# install/mfem-mgis, where the others find it. `cases`: cases of the
# measurement campaign specific to the prototype (directory cases/, built
# in build/cases).
#
# An existing build directory is reused, and reconfigured with the options
# of this script at each call. To start from scratch, delete it first.
#
# AmgX is used unless the machine sets PROTO_AMGX=OFF. It is the one given
# by the variable AMGX_DIR (the profile of the development machine defines
# it for the environment mmm-gpu), otherwise the one of the spack
# environment. A warning is printed if this AmgX does not have the patch for
# CUDA 13 (classical multigrid, see deps/spack/README.md).
#
# The MPI library is the package PROTO_MPI_PACKAGE (openmpi by default) of
# the spack environment, possibly an external package.
#
# The CUDA compilers (nvcc and its host compiler) are those of MFEM,
# chosen by mfem-mgis (see cmake/modules/MFEMMGISDeviceCode.cmake): the
# ones that an old build directory keeps in its cache are removed.
set -e
P=$(cd "$(dirname "$0")" && pwd)
source "$P/scripts/env.sh"
MFEM=$(proto_spack -e "$PROTO_SPACK_ENVIRONMENT" location -i mfem)
MPI=$(proto_spack -e "$PROTO_SPACK_ENVIRONMENT" \
  location -i "${PROTO_MPI_PACKAGE:-openmpi}")
AMGX_OPTIONS=()
if [ "${PROTO_AMGX:-ON}" = ON ]; then
  # AmgX, used by mfem-mgis through its C interface: the one of the
  # environment, or the installation given by AMGX_DIR
  AMGX=${AMGX_DIR:-$(proto_spack -e "$PROTO_SPACK_ENVIRONMENT" \
    location -i amgx)}
  # the patch for CUDA 13 is recorded in the specification of the
  # installation
  amgx_patch=$P/deps/spack/repos/spack_repo/mmm_gpu/packages/amgx/amgx-cuda13-sequence-kernel.patch
  if ! /usr/bin/grep -qs "$(sha256sum "$amgx_patch" | cut -d' ' -f1)" \
    "$AMGX/.spack/spec.json"; then
    echo "build-mfem-mgis.sh: warning: the AmgX of $AMGX is not patched" \
      "for CUDA 13 (classical multigrid), see deps/spack/README.md" >&2
  fi
  AMGX_OPTIONS=(-Denable-amgx=ON -DAMGX_DIR="$AMGX")
fi
# mpirun of the machine profile first, if any (PROTO_MPI_LAUNCHERS)
export PATH=${PROTO_MPI_LAUNCHERS:+$PROTO_MPI_LAUNCHERS:}$MPI/bin:$PATH
JOBS=${JOBS:-${PROTO_JOBS:-8}}
# launcher of the parallel tests given by find_package(MPI) (examples,
# mm-opera-hpc): the mpiexec of the machine profile, if any, instead of the
# one of the MPI library
MPIEXEC_OPTIONS=()
if [ -n "${PROTO_MPI_LAUNCHERS:-}" ]; then
  MPIEXEC_OPTIONS=(-DMPIEXEC_EXECUTABLE="$PROTO_MPI_LAUNCHERS/mpiexec")
fi
# usage: build name source-directory install(yes|no) cmake-options...
build() {
  local name=$1 source=$2 install=$3
  shift 3
  local b=$P/build/$name
  cmake -S "$source" -B "$b" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$PROTO_CC" -DCMAKE_CXX_COMPILER="$PROTO_CXX" \
    -DMPI_C_COMPILER="$MPI/bin/mpicc" -DMPI_CXX_COMPILER="$MPI/bin/mpicxx" \
    "${MPIEXEC_OPTIONS[@]}" \
    -UCMAKE_CUDA_COMPILER -UCMAKE_CUDA_HOST_COMPILER \
    "$@"
  cmake --build "$b" -j "$JOBS"
  # the tests are built by `check`
  cmake --build "$b" -j "$JOBS" --target check
  if [ "$install" = yes ]; then
    cmake --install "$b"
  fi
}
for what in ${@:-mfem-mgis examples opera cases}; do
  case "$what" in
    mfem-mgis)
      build mfem-mgis "$P/mfem-mgis" yes \
        -DMFEM_DIR="$MFEM" \
        "${AMGX_OPTIONS[@]}" \
        -DMFrontGenericInterface_DIR="$P/install/mgis/share/mgis/cmake" \
        -DTFEL_DIR="$PROTO_TFEL_INSTALL/share/tfel/cmake" \
        -DCMAKE_INSTALL_PREFIX="$P/install/mfem-mgis"
      ;;
    examples)
      build mfem-mgis-examples "$P/mfem-mgis-examples" no \
        -DMFEMMGIS_DIR="$P/install/mfem-mgis/share/mfem-mgis/cmake" \
        -DMFEM_DIR="$MFEM"
      ;;
    opera)
      build mm-opera-hpc "$P/mm-opera-hpc" no \
        -DMFEMMGIS_DIR="$P/install/mfem-mgis/share/mfem-mgis/cmake" \
        -DMFEM_DIR="$MFEM"
      ;;
    cases)
      build cases "$P/cases" no \
        -DMFEMMGIS_DIR="$P/install/mfem-mgis/share/mfem-mgis/cmake" \
        -DMFEM_DIR="$MFEM"
      ;;
    *)
      echo "usage: $0 [mfem-mgis|examples|opera|cases ...]" >&2
      exit 1
      ;;
  esac
done
