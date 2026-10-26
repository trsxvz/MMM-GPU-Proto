# TFEL settings of Jean Zay (H100 partition), sourced by build-tfel.sh, which
# defines P (root of the prototype). See deps/machines/workstation/tfel.sh
# for the variables to define.
#
# generic-parallel configurations: cuda-nvcc and openmp-host-gcc. MFront
# compiles the generated code with `-arch=native` (nvcc) and
# `-march=native` (gcc): TFEL, MGIS and mfem-mgis are built on a node of the
# partition gpu_p6 (deps/machines/jean-zay-h100/jobs/build.slurm), whose
# GPU and processors are the ones of the computations. The tests of TFEL
# are built only if PROTO_TFEL_TESTING=ON.
source "$P/scripts/env.sh"
BUILD=$PROTO_TFEL_BUILD
INSTALL=$PROTO_TFEL_INSTALL
CONFIGURATIONS="cuda-nvcc;openmp-host-gcc"
TFEL_OPTIONS=(-DCMAKE_C_COMPILER="$PROTO_CC" -DCMAKE_CXX_COMPILER="$PROTO_CXX"
  -Denable-testing="${PROTO_TFEL_TESTING:-OFF}"
  -DCUDAToolkit_ROOT="$PROTO_CUDA")
