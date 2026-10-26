# Environment of the prototype, to be sourced: `source scripts/env.sh`.
#
# PROTO_MACHINE selects the machine profile, a directory of deps/machines
# (workstation by default, the development machine of the prototype). Its
# file machine.sh defines the compilers, the libraries and the paths of the
# machine (see deps/machines/README.md):
# - PROTO_CC, PROTO_CXX, PROTO_FC: C, C++ and Fortran compilers;
# - PROTO_SPACK: spack executable, PROTO_SPACK_ENVIRONMENT: directory of the
#   spack environment providing the libraries (MFEM, hypre, MPI...);
# - PROTO_SPACK_VIEW: view of this environment;
# - PROTO_CUDA: CUDA toolkit (nvcc);
# - PROTO_TFEL_BUILD and PROTO_TFEL_INSTALL: build and installation
#   directories of TFEL;
# - PATH, LD_LIBRARY_PATH, CMAKE_PREFIX_PATH and the variables of the
#   toolchains of the machine.
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PROTO=$P
export PROTO_MACHINE=${PROTO_MACHINE:-workstation}
export PROTO_MACHINE_DIR=$P/deps/machines/$PROTO_MACHINE
if [ ! -f "$PROTO_MACHINE_DIR/machine.sh" ]; then
  echo "env.sh: unknown machine '$PROTO_MACHINE'" \
    "(no file $PROTO_MACHINE_DIR/machine.sh)" >&2
  return 1 2> /dev/null || exit 1
fi
source "$PROTO_MACHINE_DIR/machine.sh" || {
  return 1 2> /dev/null || exit 1
}
# spack of the machine, isolated from the configuration of the user
# (~/.spack): the configuration comes from the environment and from the
# instance of spack only
proto_spack() {
  SPACK_DISABLE_LOCAL_CONFIG=1 "$PROTO_SPACK" "$@"
}
