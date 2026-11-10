#!/bin/bash
# Layout of the data at the integration points, point by point or component by
# component, on one H100 and 24 cores (a quarter of a node of gpu_p6):
# cases/measure-layout.sh for the beam and the polycrystal (GPU and CPU),
# then, if Nsight Compute is found (variable NCU or PATH), the hardware
# counters of the kernels of the beam on the GPU (counters=1), after
# build.slurm. To be submitted from the deployment
# directory ROOT:
#   cd $ROOT && sbatch -A <project>@h100 \
#     MMM-GPU-Proto/deps/machines/jean-zay-h100/jobs/measure-layout.slurm
# The results are written in ROOT/MMM-GPU-Proto/results/layout-<case>-<date>
# and layout-counters-beam-<date>.
#SBATCH --job-name=mmm-layout-h100
#SBATCH --partition=gpu_p6
#SBATCH --constraint=h100
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --hint=nomultithread
#SBATCH --qos=qos_gpu_h100-dev
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
set -u
export PROTO_DEPLOY=${PROTO_DEPLOY:-$SLURM_SUBMIT_DIR}
export PROTO_MACHINE=jean-zay-h100
export CPU_PROCESSES=${CPU_PROCESSES:-24}
P=$PROTO_DEPLOY/MMM-GPU-Proto
# modules of the machine, for ncu
source "$P/scripts/env.sh" > /dev/null
status=0
for c in beam polycrystal; do
  "$P/cases/measure-layout.sh" "$c" || status=1
done
if [ -n "${NCU:-$(command -v ncu)}" ]; then
  counters=1 "$P/cases/measure-layout.sh" beam \
    "$P/results/layout-counters-beam-$(date +%Y-%m-%d-%H%M)" || status=1
else
  echo "ncu not found: no hardware counters"
fi
exit $status
