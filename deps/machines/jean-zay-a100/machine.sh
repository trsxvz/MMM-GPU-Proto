# Machine profile: Jean Zay (IDRIS), partition gpu_p5 (nodes of 8 A100
# 80 GB and 2 AMD EPYC of 32 cores, RHEL 9.6), sourced by scripts/env.sh,
# which defines P (root of the prototype). Same organisation as the profile
# jean-zay-h100, with the modules of the A100 partition (arch/a100), which
# also run on the login nodes.
#
# The prototype is deployed by scripts/offline/deploy.sh in a directory
# PROTO_DEPLOY (by default the parent directory of the prototype) which
# holds the instance of spack (spack/) and its packages (spack-packages/),
# see docs/install.md. The libraries come from the spack environment
# deps/spack/mmm-light, with the compilers, CUDA, MPI and CMake of the
# modules of IDRIS as external packages (spack-config.yaml).
module purge > /dev/null 2>&1
# Open MPI 4.1.8 of IDRIS: compiled with CUDA 13.0, PSM2 (Omni-Path) and the
# PMIx of Slurm, the variant matching the loaded compiler (gcc 14.2.0)
module load arch/a100 gcc/14.2.0 cuda/13.0.3 openmpi/4.1.8-cuda \
  cmake/3.31.4 > /dev/null 2>&1 || {
  echo "jean-zay-a100: the modules of the machine can't be loaded" >&2
  return 1
}
export PROTO_DEPLOY=${PROTO_DEPLOY:-$(dirname "$P")}
# the Open MPI of IDRIS has no mpirun (only srun): mpirun of
# scripts/mpirun-srun, first in the PATH, for the tests and the measurements
export PROTO_MPI_LAUNCHERS=$P/scripts/mpirun-srun
export PROTO_CC=gcc PROTO_CXX=g++ PROTO_FC=gfortran
export PROTO_SPACK=$PROTO_DEPLOY/spack/bin/spack
# builtin repository of spack, read by deps/spack/mmm-light/spack.yaml
export PROTO_SPACK_PACKAGES=$PROTO_DEPLOY/spack-packages
export PROTO_SPACK_ENVIRONMENT=$P/deps/spack/mmm-light
export PROTO_SPACK_VIEW=$PROTO_SPACK_ENVIRONMENT/.spack-env/view
export PROTO_CUDA=/lustre/fshomisc/sys/common/nvidia/cuda/13.0
export PROTO_TFEL_BUILD=$P/build/tfel PROTO_TFEL_INSTALL=$P/install/tfel
# AmgX is not part of mmm-light
export PROTO_AMGX=OFF
# parallel jobs of the builds: the cores of the job on the node
export PROTO_JOBS=${PROTO_JOBS:-${SLURM_CPUS_ON_NODE:-8}}
unset SPACK_ENV SPACK_ENV_VIEW
export CUDA_HOME=$PROTO_CUDA
export CMAKE_PREFIX_PATH=$PROTO_TFEL_INSTALL:$PROTO_SPACK_VIEW
export PATH=$PROTO_MPI_LAUNCHERS:$PROTO_TFEL_INSTALL/bin:$PROTO_SPACK_VIEW/bin:$PATH
export LD_LIBRARY_PATH=$PROTO_TFEL_INSTALL/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_PROC_BIND=spread OMP_PLACES=threads
# options of CMake of MGIS (build-mgis.sh): tests of the batch libraries of
# the cuda-nvcc configuration
PROTO_MGIS_OPTIONS=(-Denable-generic-parallel-cuda-tests=ON)
