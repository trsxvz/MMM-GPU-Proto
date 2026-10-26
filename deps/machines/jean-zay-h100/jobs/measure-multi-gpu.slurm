#!/bin/bash
# Strong scaling of the beam on one whole node of gpu_p6 (4 H100, 96 cores):
# cases/measure-multi-gpu.sh, after build.slurm. To be submitted from the
# deployment directory ROOT:
#   cd $ROOT && [refinement=6] sbatch -A <project>@h100 \
#     MMM-GPU-Proto/deps/machines/jean-zay-h100/jobs/measure-multi-gpu.slurm
# The results are written in ROOT/MMM-GPU-Proto/results/multi-gpu-beam-<date>.
#SBATCH --job-name=mmm-multi-gpu
#SBATCH --partition=gpu_p6
#SBATCH --constraint=h100
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --hint=nomultithread
#SBATCH --qos=qos_gpu_h100-dev
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
set -u
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
export PROTO_MACHINE=jean-zay-h100
"$PROTO_DEPLOY/MMM-GPU-Proto/cases/measure-multi-gpu.sh"
