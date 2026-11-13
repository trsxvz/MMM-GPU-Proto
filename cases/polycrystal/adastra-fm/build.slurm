#!/bin/bash
# Separate build (build-fix.sh) in FM_FIX, on a node of the MI300A partition
# (architecture of the GPUs detected by CMake, -march=native of the nodes):
# variants FM_BUILD_VARIANTS (default B).
#SBATCH --job-name=fm-build
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gpus-per-node=4
#SBATCH --time=01:30:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
log "start, build in ${FM_FIX:?}"
if FM_BUILD_VARIANTS=${FM_BUILD_VARIANTS:-B} bash "$FM_FIX/build-fix.sh" \
  > "$FM_FIX/build.log" 2>&1; then
  log "built: $(grep -E '^(A|B) built' "$FM_FIX/build.log" | tr '\n' ' ')"
else
  log "build failed (see $FM_FIX/build.log)"
  exit 1
fi
