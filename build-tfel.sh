#!/bin/bash
# Builds and installs the TFEL fork of the prototype.
# PROTO_MACHINE (see scripts/env.sh) selects the machine profile
# (workstation by default): its file deps/machines/<machine>/tfel.sh gives
# the build and installation directories, the generic-parallel
# configurations and the options of CMake specific to the machine. On the
# development machine (deps/machines/workstation/tfel.sh), PROTO_ENVIRONMENT
# selects the spack environment (mmm-gpu by default, or tfelgpu_local).
# To change the configurations of an existing build directory:
# cmake -Dgeneric-parallel-configurations="..." -U "tfel-substitutions*" .
# (with -DRAJA_DIR=<installation of RAJA>/lib/cmake/raja when adding raja,
# -DArccore_ROOT=<installation of Arccore> when adding arcane).
# The tests run in the environment of scripts/env.sh.
set -e
P=$(cd "$(dirname "$0")" && pwd)
MACHINE_DIR=$P/deps/machines/${PROTO_MACHINE:-workstation}
if [ ! -f "$MACHINE_DIR/tfel.sh" ]; then
  echo "build-tfel.sh: no file $MACHINE_DIR/tfel.sh" >&2
  exit 1
fi
source "$MACHINE_DIR/tfel.sh"
mkdir -p "$BUILD"
cd "$BUILD"
if [ ! -f CMakeCache.txt ]; then
  cmake "$P/tfel" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$INSTALL" \
    -Denable-website=OFF -Denable-python=OFF -Denable-python-bindings=OFF \
    -Denable-fortran=OFF \
    -Dgeneric-parallel-configurations="$CONFIGURATIONS" \
    "${TFEL_OPTIONS[@]}"
fi
make -j"${PROTO_JOBS:-16}" ${@:-install}
