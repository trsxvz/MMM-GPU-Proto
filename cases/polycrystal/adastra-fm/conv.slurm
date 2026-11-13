#!/bin/bash
# Convergence in mesh: RVE of 100 grains, P2, refinement FM_REFINEMENT, on
# the GPUs of the allocation only.
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --time=01:30:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
k=${FM_REFINEMENT:?}
log "start, refinement $k on $SLURM_JOB_NUM_NODES nodes, tag ${FM_TAG:-none}, ${FM_NSTEPS:-5} steps"
# FM_MESHES: directory of the meshes (default: the raw meshes of Merope and
# gmsh; $W/meshes/clean: short edges collapsed, see JOURNAL.md); the
# orientations are those of the raw meshes (same grains)
export mesh=${FM_MESHES:-$MESHES}/v100/100crystals.msh vectors=$MESHES/v100/vectors_100crystals.txt
export label=g100 refinement=$k modes=gpu nsteps=${FM_NSTEPS:-5} run_timeout=4800
# FM_TAG: suffix of the results directory (variants of the study, for
# instance 20 steps of 0.25 s: FM_TAG=-dt025 FM_NSTEPS=20 FM_DURATION=5)
out=$R/conv${FM_TAG:-}-r$k
export duration=${FM_DURATION:-$nsteps}
"$P/cases/measure-polycrystal-grains.sh" "$out"
energy "$out/energy.txt"
log "end: $(tail -1 "$out/series.txt" 2>/dev/null)"
