#!/bin/bash
# Strong scaling of the beam on one whole node of gpu_p5 (8 A100, 128 cores):
# cases/measure-multi-gpu.sh on 1, 2, 4 and 8 GPUs and 32, 64 and 128
# processes on the CPU, after build.slurm. To be submitted from the
# deployment directory ROOT:
#   cd $ROOT && [refinement=6] sbatch -A <project>@a100 \
#     MMM-GPU-Proto/deps/machines/jean-zay-a100/jobs/measure-multi-gpu.slurm
# The results are written in ROOT/MMM-GPU-Proto/results/multi-gpu-beam-<date>.
#SBATCH --job-name=mmm-multi-gpu-a100
#SBATCH --partition=gpu_p5
#SBATCH --constraint=a100
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:8
#SBATCH --hint=nomultithread
#SBATCH --qos=qos_gpu_a100-dev
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
set -u
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
export PROTO_MACHINE=jean-zay-a100
export gpu_processes=${gpu_processes:-"1 2 4 8"}
export cpu_processes=${cpu_processes:-"32 64 128"}
"$PROTO_DEPLOY/MMM-GPU-Proto/cases/measure-multi-gpu.sh"
