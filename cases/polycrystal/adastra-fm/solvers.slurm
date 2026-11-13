#!/bin/bash
# Choice of the solvers on each side, before the series: RVE of 100 grains
# (cleaned mesh), P2, refinement FM_REFINEMENT (1 by default), FM_NSTEPS
# steps of 1 s (10 by default), on all the GPUs (FM_MODE=gpu) or all the cores
# (FM_MODE=cpu) of one node, for the configurations FM_CONFIGS among (hypre
# only on Adastra; NewtonK0: tangent applied without assembly, CG
# preconditioned by one cycle of BoomerAMG built once on the elastic
# operator; Newton: PCG and BoomerAMG rebuilt at each iteration):
# - k0ew1: NewtonK0, Eisenstat-Walker 1 from 0.1 (best on Jean Zay);
# - k0ew1-001: NewtonK0, Eisenstat-Walker 1 from 0.01;
# - k0ew2: NewtonK0, Eisenstat-Walker 2 from 0.1;
# - newtonew1: Newton, Eisenstat-Walker 1 from 0.1;
# - k0ew1-nogtl: k0ew1 with MPICH_GPU_SUPPORT_ENABLED=0 (CPU only: the
#   support of the GPU buffers by cray-mpich is enabled by the machine
#   profile for all the computations).
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
k=${FM_REFINEMENT:-1}
mode=${FM_MODE:?FM_MODE (gpu or cpu) required}
out=$R/solvers-r$k
log "start, 100 grains (cleaned), refinement $k, $mode: ${FM_CONFIGS:?}"
export mesh=$W/meshes/clean/v100/100crystals.msh vectors=$MESHES/v100/vectors_100crystals.txt
export label=g100c refinement=$k nsteps=${FM_NSTEPS:-10} modes=$mode
ew1="--adaptive-linear-solver-tolerance EisenstatWalker1"
ew2="--adaptive-linear-solver-tolerance EisenstatWalker2"
declare -A configs=(
  [k0ew1]="--non-linear-solver NewtonK0 $ew1 --adaptive-linear-solver-initial-tolerance 0.1"
  [k0ew1-001]="--non-linear-solver NewtonK0 $ew1 --adaptive-linear-solver-initial-tolerance 0.01"
  [k0ew2]="--non-linear-solver NewtonK0 $ew2 --adaptive-linear-solver-initial-tolerance 0.1"
  [newtonew1]="--non-linear-solver Newton $ew1 --adaptive-linear-solver-initial-tolerance 0.1"
  [k0ew1-nogtl]="--non-linear-solver NewtonK0 $ew1 --adaptive-linear-solver-initial-tolerance 0.1"
)
for c in $FM_CONFIGS; do
  export solver_label=$c solver_gpu="${configs[$c]}" solver_cpu="${configs[$c]}"
  if [ "$c" = k0ew1-nogtl ]; then
    export cpu_environment=MPICH_GPU_SUPPORT_ENABLED=0
  else
    export cpu_environment=""
  fi
  # one results directory per job: series.txt is appended to
  "$P/cases/measure-polycrystal-grains.sh" "$out/$mode-$SLURM_JOB_ID"
  log "$mode $c: $(tail -1 "$out/$mode-$SLURM_JOB_ID/series.txt")"
done
log end
