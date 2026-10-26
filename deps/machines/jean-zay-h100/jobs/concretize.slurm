#!/bin/bash
# Concretization of the spack environment mmm-light for the H100 partition
# of Jean Zay (`deploy.sh concretize`), on the partition compil_h100: spack
# checks that the compilers run, and the gcc of the modules of arch/h100
# does not run on the login nodes (Cascade Lake). To be submitted from the
# deployment directory ROOT:
#   cd $ROOT && sbatch -A <project>@h100 \
#     MMM-GPU-Proto/deps/machines/jean-zay-h100/jobs/concretize.slurm
# The file spack.lock to give to `bundle.sh mirror` is then
# ROOT/MMM-GPU-Proto/deps/spack/mmm-light/spack.lock.
#SBATCH --job-name=mmm-concretize
#SBATCH --partition=compil_h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
set -eu
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
"$PROTO_DEPLOY/MMM-GPU-Proto/scripts/offline/deploy.sh" concretize jean-zay-h100
