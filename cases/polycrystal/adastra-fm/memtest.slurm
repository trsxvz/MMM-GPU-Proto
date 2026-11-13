#!/bin/bash
# Reading of the mesh of FM_GRAINS grains (cleaned) by 96 processes on one
# node (as the CPU computations of the series): memory (MaxRSS of the step)
# and time of Mesh::Constructor; P1, no refinement, one step of 0.01 s.
#SBATCH --job-name=fm-memtest
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
n=${FM_GRAINS:-6400}
out=$R/memtest-g$n
log "start, $n grains"
export mesh=$W/meshes/clean/v$n/${n}crystals.msh vectors=$MESHES/v$n/vectors_${n}crystals.txt
export label=g${n}c refinement=0 order=1 nsteps=1 duration=0.01 modes=cpu run_timeout_cpu=2400
"$P/cases/measure-polycrystal-grains.sh" "$out"
sacct -j "$SLURM_JOB_ID" -P --format=JobID,MaxRSS,MaxRSSTask,AveRSS,Elapsed,State > "$out/memory.txt"
log "end: $(tail -1 "$out/series.txt"); $(grep -E '\.[0-9]+\|' "$out/memory.txt" | tail -1)"
