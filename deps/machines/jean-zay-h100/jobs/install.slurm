#!/bin/bash
# Installation of the spack environment mmm-light on Jean Zay, on the
# partition compil_h100 (the compilation node of the H100 partition, no hours
# consumed), from the mirror of the sources extracted by `deploy.sh install`
# (docs/install.md). The compilers of the modules of the H100 partition
# (arch/h100) are built for its processors (Sapphire Rapids): they stop
# with "Illegal instruction" on the login nodes and the other partitions
# (Cascade Lake). To be submitted from the deployment directory ROOT (parent
# directory of the prototype):
#   cd $ROOT && sbatch -A <project>@h100 \
#     MMM-GPU-Proto/deps/machines/jean-zay-h100/jobs/install.slurm
# With PROTO_BUILDCACHE (for instance $WORK/mmm-buildcache) exported before
# sbatch, the installed packages are also copied in this binary cache.
#SBATCH --job-name=mmm-install
#SBATCH --partition=compil_h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=nomultithread
#SBATCH --time=10:00:00
#SBATCH --output=%x-%j.out
set -eu
# builds of spack in the temporary space of the job
export TMPDIR=${JOBSCRATCH:-${TMPDIR:-/tmp}}
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
# MFEM compiled with CUDA: up to several GB of memory per file for ptxas
export PROTO_JOBS=${PROTO_JOBS:-16}
"$PROTO_DEPLOY/MMM-GPU-Proto/scripts/offline/deploy.sh" install jean-zay-h100
