#!/bin/bash
# Points of a scalability study on the same number of nodes, one after the
# other in one allocation: one wait in the queue instead of one per point
# for the large numbers of nodes. The computations are the ones of
# series.slurm, with its variables.
# - FM_GPU_POINTS: "grains[:computations]" on the GPUs, with the parts made
#   and checked beforehand (scaling.sh nodes).
# - FM_CPU_POINTS: "grains[:keep]" on the cores. The parts of each one are
#   written and checked here unless they already are (partition.slurm, with
#   its variables), read by the computation, then removed unless the point
#   ends with `:keep`: only one set weighs on the quota of files. The nodes
#   wait during the partitioning (2 to 8 minutes).
# A point whose parts are missing is left out and the next ones run.
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
chain=$SLURM_SUBMIT_DIR/fm-chain
log "start, $SLURM_JOB_NUM_NODES nodes, GPU ${FM_GPU_POINTS:-none}, CPU ${FM_CPU_POINTS:-none}"
for point in ${FM_GPU_POINTS:-}; do
  IFS=: read -r n computations <<< "$point"
  d=$W/meshes/part/v$n-np$((4 * SLURM_JOB_NUM_NODES))
  if [ ! -e "$d/checked" ]; then
    log "GPU, $n grains: no checked parts in $d"
    continue
  fi
  FM_GRAINS=$n FM_MESH_GPU=$d/mesh. FM_GPU_REPS=${computations:-1} FM_CPU_REPS=0 \
    MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels bash "$chain/series.slurm"
done
for point in ${FM_CPU_POINTS:-}; do
  IFS=: read -r n keep <<< "$point"
  np=$((96 * SLURM_JOB_NUM_NODES))
  if FM_PARTITIONS=$n:$np bash "$chain/partition.slurm"; then
    FM_GRAINS=$n FM_MESH_CPU=$W/meshes/part/v$n-np$np/mesh. FM_GPU_REPS=0 \
      FM_CPU_REPS=1 bash "$chain/series.slurm"
  fi
  [ "$keep" = keep ] ||
    FM_PARTITIONS= FM_REMOVE_PARTITIONS=$n:$np bash "$chain/partition.slurm"
done
log end
