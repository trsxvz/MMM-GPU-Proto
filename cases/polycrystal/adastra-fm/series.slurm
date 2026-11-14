#!/bin/bash
# Series in number of grains: RVE of FM_GRAINS grains (cleaned mesh), P2,
# refinement FM_SERIES_REFINEMENT (1: retained from the convergence study,
# see JOURNAL.md), FM_NSTEPS steps of 1 s, on all the GPUs (FM_GPU_REPS
# times) then all the cores (FM_CPU_REPS times) of the nodes of the
# allocation, with the solvers FM_SOLVER_GPU and FM_SOLVER_CPU (the best of
# each side in the study of the solvers). Time of each step in
# <computation>.steps, energy of each computation in series.txt.
# FM_MESH_GPU, FM_MESH_CPU: other meshes for the GPU and the CPU, for
# instance meshes partitioned for their numbers of processes
# (meshes/partition/, read with --mesh_mode Restart and
# FM_SERIES_REFINEMENT=0).
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
n=${FM_GRAINS:?}
k=${FM_SERIES_REFINEMENT:-1}
out=$R/series-g$n-$SLURM_JOB_ID
mkdir -p "$out"
log "start, $n grains (cleaned), refinement $k, $SLURM_JOB_NUM_NODES nodes," \
  "${FM_NSTEPS:-20} steps, GPU x${FM_GPU_REPS:-1} ($FM_SOLVER_GPU)," \
  "CPU x${FM_CPU_REPS:-1} ($FM_SOLVER_CPU), scripts $P"
export mesh=$W/meshes/clean/v$n/${n}crystals.msh
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
export refinement=$k nsteps=${FM_NSTEPS:-20} solver_label=best
export solver_gpu=$FM_SOLVER_GPU solver_cpu=$FM_SOLVER_CPU
export run_timeout_gpu=${FM_TIMEOUT_GPU:-2400}
export run_timeout_cpu=${FM_TIMEOUT_CPU:-10800}
export cpu_environment=${FM_CPU_ENVIRONMENT:-}
for i in $(seq "${FM_GPU_REPS:-1}"); do
  mesh=${FM_MESH_GPU:-$mesh} label=g${n}c-rep$i modes=gpu \
    "$P/cases/measure-polycrystal-grains.sh" "$out"
  log "gpu rep $i: $(tail -1 "$out/series.txt")"
done
for i in $(seq "${FM_CPU_REPS:-1}"); do
  mesh=${FM_MESH_CPU:-$mesh} label=g${n}c-rep$i modes=cpu \
    "$P/cases/measure-polycrystal-grains.sh" "$out"
  log "cpu rep $i: $(tail -1 "$out/series.txt")"
done
energy "$out/energy.txt"
log end
