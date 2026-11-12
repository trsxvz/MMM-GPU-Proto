#!/bin/bash
# Some tests of a build of the prototype only, after build.slurm, for
# instance to run again a test which failed. To be submitted from the
# deployment directory ROOT:
#   cd $ROOT && sbatch -A <project> \
#     MMM-GPU-Proto/deps/machines/adastra-mi300/jobs/ctest.slurm \
#     BUILD REGEX [ctest options...]
# where BUILD is a directory of ROOT/MMM-GPU-Proto/build (mfem-mgis,
# mfem-mgis-examples, mm-opera-hpc, cases...) and REGEX selects the tests
# (ctest -R).
#SBATCH --job-name=mmm-ctest-mi300
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gpus-per-node=4
#SBATCH --time=00:30:00
#SBATCH --output=%x-%j.out
set -u
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY ftp_proxy FTP_PROXY
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
export PROTO_MACHINE=adastra-mi300
P=$PROTO_DEPLOY/MMM-GPU-Proto
source "$P/scripts/env.sh"
build=${1:?usage: ctest.slurm BUILD REGEX [ctest options...]}
regex=${2:?usage: ctest.slurm BUILD REGEX [ctest options...]}
shift 2
cd "$P/build/$build" && ctest -R "$regex" --output-on-failure "$@"
