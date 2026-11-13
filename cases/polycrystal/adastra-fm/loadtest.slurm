#!/bin/bash
# Load of the GPUs: RVEs of FM_LOAD_GRAINS grains (cleaned meshes, refinement
# 1, P2) on the 4 MI300A of one node, FM_NSTEPS steps of 1 s, best solver of
# the GPU. Memory of the node sampled every 5 s (free: the memory of the
# MI300A is shared by the CPU and the GPU) and of the GPUs (rocm-smi) in
# <results>/memory-<label>.txt.
#SBATCH --job-name=fm-loadtest
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --time=01:30:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
out=$R/loadtest-$SLURM_JOB_ID
mkdir -p "$out"
log "start, ${FM_LOAD_GRAINS:?} grains, scripts $P"
export refinement=1 nsteps=${FM_NSTEPS:-10} modes=gpu solver_label=best
export solver_gpu=${FM_SOLVER_GPU:?} run_timeout_gpu=2400
for n in $FM_LOAD_GRAINS; do
  export mesh=$W/meshes/clean/v$n/${n}crystals.msh
  export vectors=$MESHES/v$n/vectors_${n}crystals.txt label=g${n}c
  (while :; do
     echo "$(date +%s) $(free -m | awk '/^Mem/ {print "used_MB", $3, "total_MB", $2}')" \
       "$(rocm-smi --showmeminfo vram --csv 2> /dev/null |
          awk -F, 'NR > 1 {printf "gpu%s_used %s ", $1, $3}')"
     sleep 5
   done) > "$out/memory-g${n}c.txt" 2>&1 &
  sampler=$!
  "$P/cases/measure-polycrystal-grains.sh" "$out"
  kill "$sampler"
  log "$n grains: $(tail -1 "$out/series.txt"); max used" \
    "$(awk '{if ($3 > m) m = $3} END {print m}' "$out/memory-g${n}c.txt") MB"
done
log end
