#!/bin/bash
# Build of MFEM (deploy.sh mfem), then build and tests of TFEL, MGIS,
# mfem-mgis, the examples, mm-opera-hpc and the cases of the prototype on a
# node of the partition mi300 (MI300A): MFront compiles the generated code
# for the processors of the node (`-march=native`) and for the MI300A, and
# the tests use one APU. To be submitted from the deployment directory ROOT,
# after `deploy.sh setup`:
#   cd $ROOT && sbatch -A <project> \
#     MMM-GPU-Proto/deps/machines/adastra-mi300/jobs/build.slurm [what...]
# where `what` is given to build-mfem-mgis.sh (mfem-mgis, examples, opera,
# cases; all by default).
#SBATCH --job-name=mmm-build-mi300
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gpus-per-node=4
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
set -eu
# no access to the internet from the deployment (docs/install.md)
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY ftp_proxy FTP_PROXY
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
D=$PROTO_DEPLOY/MMM-GPU-Proto/scripts/offline/deploy.sh
"$D" mfem adastra-mi300
"$D" build adastra-mi300 "$@"
