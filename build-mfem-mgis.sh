#!/bin/bash
# Compile, teste et installe mfem-mgis, puis compile et teste
# mfem-mgis-examples, mm-opera-hpc et les cas de cases/ (dépôts `mfem-mgis/`,
# `mfem-mgis-examples/` et `mm-opera-hpc/`, branche `gpu-proto`), contre
# MFEM 4.10 de l'environnement spack `mmm-gpu`, le MGIS installé par
# build-mgis.sh et le TFEL installé par build-tfel.sh.
#
# Les autres constructions ne passent pas par ce script : MFEM compilé avec
# HIP pour un GPU AMD (scripts/build-hip-amd.sh), et les constructions qui
# valident l'option enable-kokkos (build/mfem-mgis-kokkos-off et
# build/mfem-mgis-kokkos, voir docs/validation.md).
#
# usage : ./build-mfem-mgis.sh [mfem-mgis|examples|opera|cas ...]
# Sans argument : les quatre, dans cet ordre. mfem-mgis est installé dans
# install/mfem-mgis, où les autres le trouvent. `cas` : cas de la campagne
# de mesures propres au prototype (dossier cases/).
#
# Un dossier de construction existant est réutilisé, et reconfiguré avec
# les options de ce script à chaque appel. Pour repartir de zéro, le
# supprimer avant.
#
# AmgX est celui donné par la variable AMGX_DIR (scripts/env.sh la définit
# pour l'environnement mmm-gpu), sinon celui de l'environnement spack. Un
# avertissement est affiché si cet AmgX n'a pas le patch pour CUDA 13
# (multigrille classique, voir deps/spack/README.md).
#
# Les compilateurs CUDA (nvcc et son compilateur hôte) sont ceux de MFEM,
# choisis par mfem-mgis (voir cmake/modules/MFEMMGISDeviceCode.cmake) :
# ceux qu'un ancien dossier garde dans son cache sont retirés.
set -e
P=$(cd "$(dirname "$0")" && pwd)
source "$P/scripts/env.sh"
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
MFEM=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
  location -i mfem)
MPI=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
  location -i openmpi)
# AmgX, used by mfem-mgis through its C interface: the one of the
# environment, or the installation given by AMGX_DIR
AMGX=${AMGX_DIR:-$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" \
  -e "$P/deps/spack/mmm-gpu" location -i amgx)}
# the patch for CUDA 13 is recorded in the specification of the installation
amgx_patch=$P/deps/spack/repos/spack_repo/mmm_gpu/packages/amgx/amgx-cuda13-sequence-kernel.patch
if ! /usr/bin/grep -qs "$(sha256sum "$amgx_patch" | cut -d' ' -f1)" \
  "$AMGX/.spack/spec.json"; then
  echo "build-mfem-mgis.sh: warning: the AmgX of $AMGX is not patched" \
    "for CUDA 13 (classical multigrid), see deps/spack/README.md" >&2
fi
export PATH=$MPI/bin:$PATH
JOBS=${JOBS:-8}
# usage: build name source-directory install(yes|no) cmake-options...
build() {
  local name=$1 source=$2 install=$3
  shift 3
  local b=$P/build/$name
  cmake -S "$source" -B "$b" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
    -DMPI_C_COMPILER="$MPI/bin/mpicc" -DMPI_CXX_COMPILER="$MPI/bin/mpicxx" \
    -UCMAKE_CUDA_COMPILER -UCMAKE_CUDA_HOST_COMPILER \
    "$@"
  cmake --build "$b" -j "$JOBS"
  # the tests are built by `check`
  cmake --build "$b" -j "$JOBS" --target check
  if [ "$install" = yes ]; then
    cmake --install "$b"
  fi
}
for what in ${@:-mfem-mgis examples opera cas}; do
  case "$what" in
    mfem-mgis)
      build mfem-mgis "$P/mfem-mgis" yes \
        -DMFEM_DIR="$MFEM" \
        -Denable-amgx=ON -DAMGX_DIR="$AMGX" \
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
    cas)
      build cas "$P/cas" no \
        -DMFEMMGIS_DIR="$P/install/mfem-mgis/share/mfem-mgis/cmake" \
        -DMFEM_DIR="$MFEM"
      ;;
    *)
      echo "usage: $0 [mfem-mgis|examples|opera|cas ...]" >&2
      exit 1
      ;;
  esac
done
