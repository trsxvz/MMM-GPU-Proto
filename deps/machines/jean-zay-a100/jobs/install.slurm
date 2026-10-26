#!/bin/bash
# Installation of the spack environment mmm-light on Jean Zay for the A100
# partition, on the CPU partition (the compilers of the modules of the A100
# partition run there), from the mirror of the sources extracted by
# `deploy.sh install` (docs/install.md). To be submitted from the deployment
# directory ROOT (parent directory of the prototype):
#   cd $ROOT && sbatch -A <project>@cpu \
#     MMM-GPU-Proto/deps/machines/jean-zay-a100/jobs/install.slurm
# With PROTO_BUILDCACHE (for instance $WORK/mmm-buildcache-a100) exported
# before sbatch, the installed packages are also copied in this binary
# cache.
#SBATCH --job-name=mmm-install-a100
#SBATCH --partition=cpu_p1
#SBATCH --qos=qos_cpu-dev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --hint=nomultithread
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
set -eu
# builds of spack in the temporary space of the job
export TMPDIR=${JOBSCRATCH:-${TMPDIR:-/tmp}}
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
# MFEM compiled with CUDA: up to several GB of memory per file for ptxas
# (160 GB on the nodes of cpu_p1)
export PROTO_JOBS=${PROTO_JOBS:-16}
"$PROTO_DEPLOY/MMM-GPU-Proto/scripts/offline/deploy.sh" install jean-zay-a100
