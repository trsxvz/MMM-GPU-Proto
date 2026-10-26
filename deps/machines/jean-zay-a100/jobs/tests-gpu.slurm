#!/bin/bash
# Test suite of mfem-mgis on one A100 (scripts/mfem-mgis-gpu-suite.sh, modes
# kernels and default), the tests whose name selects a device, and the tests
# of the examples, of mm-opera-hpc and of the cases, after build.slurm. To
# be submitted from the deployment directory ROOT:
#   cd $ROOT && sbatch -A <project>@a100 \
#     MMM-GPU-Proto/deps/machines/jean-zay-a100/jobs/tests-gpu.slurm
# The outputs of ctest are in ROOT/MMM-GPU-Proto/build/*/Testing.
#SBATCH --job-name=mmm-tests-gpu-a100
#SBATCH --partition=gpu_p5
#SBATCH --constraint=a100
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:1
# 16 cores (32 logical processors): billed as 2 GPUs on gpu_p5, where 16
# logical processors are billed as one GPU
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
set -u
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
export PROTO_MACHINE=jean-zay-a100
P=$PROTO_DEPLOY/MMM-GPU-Proto
source "$P/scripts/env.sh"
status=0
for mode in kernels default; do
  echo "== mfem-mgis on the GPU, mode $mode"
  "$P/scripts/mfem-mgis-gpu-suite.sh" "$mode" || status=1
done
echo "== mfem-mgis, tests selecting a device"
(cd "$P/build/mfem-mgis" &&
  ctest -j 1 -R "cuda|gpu|GPU|HypreDevicePoolSize" --output-on-failure) ||
  status=1
for b in mfem-mgis-examples mm-opera-hpc cases; do
  echo "== $b"
  (cd "$P/build/$b" && ctest -j 1 --output-on-failure) || status=1
done
exit $status
