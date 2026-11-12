#!/bin/bash
# Strong scaling on one whole node of the partition mi300 (4 MI300A, 96
# cores): cases/measure-multi-gpu.sh on 1, 2 and 4 GPUs and 24, 48 and 96
# processes on the CPU, after build.slurm. To be submitted from the
# deployment directory ROOT:
#   cd $ROOT && [multi_gpu_case=polycrystal] [refinement=n] \
#     sbatch -A <project> \
#     MMM-GPU-Proto/deps/machines/adastra-mi300/jobs/measure-multi-gpu.slurm
# The results are written in ROOT/MMM-GPU-Proto/results/multi-gpu-<case>-<date>.
#SBATCH --job-name=mmm-multi-gpu-mi300
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
set -u
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY ftp_proxy FTP_PROXY
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
export PROTO_MACHINE=adastra-mi300
export gpu_processes=${gpu_processes:-"1 2 4"}
export cpu_processes=${cpu_processes:-"24 48 96"}
"$PROTO_DEPLOY/MMM-GPU-Proto/cases/measure-multi-gpu.sh"
