#!/bin/bash
# Kernel trace (rocprofv3) of the polycrystal of FM_GRAINS grains (cleaned
# mesh, r1, P2) on the GPUs of one node, FM_NSTEPS steps, best solver of the
# GPU: fragmentation of the work in kernels (one per grain and per phase),
# analysed by fm-chain/kernels.py. A computation without the profiler first,
# for its overhead.
#SBATCH --job-name=fm-profile
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
n=${FM_GRAINS:-400}
out=$R/profile-g$n-$SLURM_JOB_ID
mkdir -p "$out"
log "start, $n grains, scripts $P"
export mesh=$W/meshes/clean/v$n/${n}crystals.msh
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
export refinement=1 nsteps=${FM_NSTEPS:-3} modes=gpu solver_gpu=${FM_SOLVER_GPU:?}
export run_timeout_gpu=2400
label=g${n}c solver_label=plain "$P/cases/measure-polycrystal-grains.sh" "$out"
log "plain: $(tail -1 "$out/series.txt")"
profiler="/opt/rocm-6.4.3/bin/rocprofv3 --kernel-trace --stats --output-format csv -d $out/prof --" \
  label=g${n}c solver_label=rocprof "$P/cases/measure-polycrystal-grains.sh" "$out"
log "rocprof: $(tail -1 "$out/series.txt")"
newton=$(tail -1 "$out/series.txt" | awk '{print $11}')
python3 -I "$SLURM_SUBMIT_DIR/fm-chain/kernels.py" "$out/prof" "$newton" > "$out/kernels.txt" 2>&1
log "analysis: $(head -4 "$out/kernels.txt" | tr '\n' ' ')"
