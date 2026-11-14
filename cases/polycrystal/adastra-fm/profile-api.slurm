#!/bin/bash
# Where the GPUs wait: polycrystal of FM_GRAINS grains (cleaned mesh, r1, P2,
# FM_NSTEPS steps of FM_DURATION s, variables executable and extra_options
# of measure-polycrystal-grains.sh) on the GPUs of one node, traced by
# rocprofv3 (kernels, HIP runtime, memory copies), analysed by
# fm-chain/api-gaps.py (idle time of the GPUs attributed to the calls of the
# HIP runtime or to the host) and fm-chain/kernels.py.
#SBATCH --job-name=fm-profile-api
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
out=$R/profile-api-g$n-$SLURM_JOB_ID
mkdir -p "$out"
log "start, $n grains, scripts $P"
export mesh=$W/meshes/clean/v$n/${n}crystals.msh
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
export refinement=1 nsteps=${FM_NSTEPS:-2} duration=${FM_DURATION:-2}
export modes=gpu solver_gpu=${FM_SOLVER_GPU:?} run_timeout_gpu=2400
profiler="/opt/rocm-6.4.3/bin/rocprofv3 --kernel-trace --hip-runtime-trace --memory-copy-trace --output-format csv -d $out/prof --" \
  label=g${n}c solver_label=api "$P/cases/measure-polycrystal-grains.sh" "$out"
log "traced: $(tail -1 "$out/series.txt")"
python3 -I "$SLURM_SUBMIT_DIR/fm-chain/api-gaps.py" "$out/prof" > "$out/api-gaps.txt" 2>&1
newton=$(tail -1 "$out/series.txt" | awk '{print $11}')
python3 -I "$SLURM_SUBMIT_DIR/fm-chain/kernels.py" "$out/prof" "$newton" > "$out/kernels.txt" 2>&1
log "analysis: $(grep -A4 '^all processes' "$out/api-gaps.txt" | tr '\n' ' ' | cut -c1-600)"
