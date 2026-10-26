#!/bin/bash
# Build and tests of TFEL, MGIS, mfem-mgis, the examples, mm-opera-hpc and
# the cases of the prototype on a node of the partition gpu_p6 (H100): MFront
# compiles the generated code for the GPU and the processors of the node
# (`-arch=native`, `-march=native`), and the tests use one GPU. To be
# submitted from the deployment directory ROOT, after install.slurm:
#   cd $ROOT && sbatch -A <project>@h100 \
#     MMM-GPU-Proto/deps/machines/jean-zay-h100/jobs/build.slurm [what...]
# where `what` is given to build-mfem-mgis.sh (mfem-mgis, examples, opera,
# cases; all by default).
#SBATCH --job-name=mmm-build
#SBATCH --partition=gpu_p6
#SBATCH --constraint=h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --hint=nomultithread
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
set -eu
export TMPDIR=${JOBSCRATCH:-${TMPDIR:-/tmp}}
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
"$PROTO_DEPLOY/MMM-GPU-Proto/scripts/offline/deploy.sh" build jean-zay-h100 "$@"
