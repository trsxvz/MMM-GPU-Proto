#!/bin/bash
# Compile, teste et installe mfem-mgis, puis compile et teste
# mfem-mgis-examples et mm-opera-hpc (dépôts `mfem-mgis/`,
# `mfem-mgis-examples/` et `mm-opera-hpc/`, branche `gpu-proto`), contre
# MFEM 4.10 de l'environnement spack `mmm-gpu`, le MGIS installé par
# build-mgis.sh et le TFEL installé par build-tfel.sh.
#
# usage : ./build-mfem-mgis.sh [mfem-mgis|examples|opera ...]
# Sans argument : les trois, dans cet ordre. mfem-mgis est installé dans
# install/mfem-mgis, où les deux autres le trouvent.
#
# Un dossier de construction existant est réutilisé. Pour repartir de zéro,
# le supprimer avant.
set -e
P=$(cd "$(dirname "$0")" && pwd)
source "$P/scripts/env.sh"
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
MFEM=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
  location -i mfem)
MPI=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
  location -i openmpi)
export PATH=$MPI/bin:$PATH
JOBS=${JOBS:-8}
# usage: build name source-directory install(yes|no) cmake-options...
build() {
  local name=$1 source=$2 install=$3
  shift 3
  local b=$P/build/$name
  if [ ! -f "$b/CMakeCache.txt" ]; then
    cmake -S "$source" -B "$b" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
      -DMPI_C_COMPILER="$MPI/bin/mpicc" -DMPI_CXX_COMPILER="$MPI/bin/mpicxx" \
      "$@"
  fi
  cmake --build "$b" -j "$JOBS"
  # the tests are built by `check`
  cmake --build "$b" -j "$JOBS" --target check
  if [ "$install" = yes ]; then
    cmake --install "$b"
  fi
}
for what in ${@:-mfem-mgis examples opera}; do
  case "$what" in
    mfem-mgis)
      build mfem-mgis "$P/mfem-mgis" yes \
        -DMFEM_DIR="$MFEM" \
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
    *)
      echo "usage: $0 [mfem-mgis|examples|opera ...]" >&2
      exit 1
      ;;
  esac
done
