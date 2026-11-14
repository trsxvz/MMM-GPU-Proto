#!/bin/bash
# Prediction by linear extrapolation on the polycrystal: RVE of FM_GRAINS
# grains (cleaned mesh, P2, refinement 1), FM_NSTEPS steps over 20 s, one
# material for all the grains, with the options FM_COMMON (default: the
# line search), without and with `--prediction LinearExtrapolation`, the
# two variants being run one after the other at each repetition, on the
# same nodes: FM_GPU_REPS times on all the GPUs, then FM_CPU_REPS times on
# all the cores. The solvers are FM_SOLVER (default: NewtonK0, first choice
# of Eisenstat and Walker from 0.01, on both sides). FM_VARIANTS (default
# "ls ls-pred") selects the variants. A line per computation in
# <results>/series.txt (see cases/measure-polycrystal-grains.sh).
#SBATCH --job-name=fm-pred-ab
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
n=${FM_GRAINS:?}
out=$R/pred-ab-g$n-$SLURM_JOB_ID
mkdir -p "$out"
solver=${FM_SOLVER:---non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01}
common=${FM_COMMON:---line-search Backtracking}
export executable=${FM_EXECUTABLE:?}
log "start, $n grains (cleaned), $SLURM_JOB_NUM_NODES nodes," \
  "${FM_NSTEPS:-20} steps, GPU x${FM_GPU_REPS:-2}, CPU x${FM_CPU_REPS:-0}," \
  "solvers $solver, options $common, executable $executable, scripts $P"
export mesh=$W/meshes/clean/v$n/${n}crystals.msh
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
export refinement=1 nsteps=${FM_NSTEPS:-20} duration=20 solver_label=best
export solver_gpu=$solver solver_cpu=$solver
export run_timeout_gpu=${FM_TIMEOUT_GPU:-2400}
export run_timeout_cpu=${FM_TIMEOUT_CPU:-10800}
# usage: run_variants mode repetition
run_variants() {
  local v options
  for v in ${FM_VARIANTS:-ls ls-pred}; do
    options="--single-material $common"
    [ "$v" = ls-pred ] && options="$options --prediction LinearExtrapolation"
    extra_options=$options label=g${n}c-$v-rep$2 modes=$1 \
      "$P/cases/measure-polycrystal-grains.sh" "$out"
    log "$1 $v rep $2: $(tail -1 "$out/series.txt")"
  done
}
for i in $(seq "${FM_GPU_REPS:-2}"); do
  run_variants gpu "$i"
done
for i in $(seq "${FM_CPU_REPS:-0}"); do
  run_variants cpu "$i"
done
energy "$out/energy.txt"
log end
