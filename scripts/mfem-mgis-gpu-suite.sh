#!/bin/bash
# Runs the test suite of mfem-mgis on the GPU (MFEM_MGIS_DEVICE=cuda), one
# test at a time, since simultaneous tests share the GPU and fail. The
# stiffness matrices computed by kernels are added in a fixed order
# (MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered), so that the results, and
# the outcome of the tests, are the same at each run. hypre reserves one
# GiB per process on the GPU (MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=1) so that
# the tests on four processes fit in its memory.
#
# Modes:
# - kernels (used when no mode is given): all the computations forced to
#   the kernels of the batch path;
# - default: the computations selected by each test.
#
# The tests whose name selects a device (cuda, gpu, debug) are excluded:
# they run with `ctest -R "cuda|gpu|GPU|AmgX|HypreDevicePoolSize|okkos"`.
# `DeviceForallTest-cpu` and `DeviceOptionTest-none`, which require the
# absence of device, are expected to fail. `TwoLayerCube-1-MatrixFreeExecution`,
# which checks that the option `--execution` gives no warning about the
# variables of placement, is excluded: this suite sets them.
#
# usage: scripts/mfem-mgis-gpu-suite.sh [default|kernels] [build directory]
set -u
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mode=${1:-kernels}
build=${2:-$P/build/mfem-mgis}
source "$P/scripts/env.sh" > /dev/null
unset OMP_PROC_BIND OMP_PLACES
settings=(MFEM_MGIS_DEVICE=cuda
  MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=1
  MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered)
case "$mode" in
  default) ;;
  kernels)
    settings+=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
      MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
      MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
      MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
    ;;
  *)
    echo "usage: $0 [default|kernels] [build directory]" >&2
    exit 1
    ;;
esac
cd "$build" || exit 1
env "${settings[@]}" ctest -j 1 \
  -E "cuda|gpu|GPU|debug|^TwoLayerCube-1-MatrixFreeExecution$" --output-on-failure
