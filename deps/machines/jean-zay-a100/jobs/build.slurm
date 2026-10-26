#!/bin/bash
# Build and tests of TFEL, MGIS, mfem-mgis, the examples, mm-opera-hpc and
# the cases of the prototype on a node of the partition gpu_p5 (A100): MFront
# compiles the generated code for the GPU and the processors of the node
# (`-arch=native`, `-march=native`), and the tests use one GPU. To be
# submitted from the deployment directory ROOT, after install.slurm:
#   cd $ROOT && sbatch -A <project>@a100 \
#     MMM-GPU-Proto/deps/machines/jean-zay-a100/jobs/build.slurm [what...]
# where `what` is given to build-mfem-mgis.sh (mfem-mgis, examples, opera,
# cases; all by default).
#SBATCH --job-name=mmm-build-a100
#SBATCH --partition=gpu_p5
#SBATCH --constraint=a100
#SBATCH --nodes=1
# 4 tasks: the parallel tests of mfem-mgis run `mpirun -n 2` (up to 4
# processes), which needs as many slots in the job
#SBATCH --ntasks=4
#SBATCH --gres=gpu:1
# 16 cores (32 logical processors): billed as 2 GPUs on gpu_p5, where 16
# logical processors are billed as one GPU
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
set -eu
export TMPDIR=${JOBSCRATCH:-${TMPDIR:-/tmp}}
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
"$PROTO_DEPLOY/MMM-GPU-Proto/scripts/offline/deploy.sh" build jean-zay-a100 "$@"
