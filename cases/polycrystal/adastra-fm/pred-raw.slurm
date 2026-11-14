#!/bin/bash
# Prediction by linear extrapolation on a case whose behaviour fails to
# integrate at some estimates: RVE of 100 grains, raw mesh (P2, refinement
# 1), 5 steps of 1 s, one material for all the grains, on all the cores,
# NewtonK0 with the first choice of Eisenstat and Walker from 0.01 (the
# case i3-raw of ls-identity.slurm), with the options FM_OPTIONS (default:
# the line search and the extrapolation) and the executable FM_EXECUTABLE.
# Logged: the line of the computation (see
# cases/measure-polycrystal-grains.sh), the extrapolated estimates
# discarded (the resolution then starts again from the default prediction),
# the estimates moved after an integration failure (without line search)
# and the failed resolutions. FM_LABEL names the computation.
#SBATCH --job-name=fm-pred-raw
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
out=$R/pred-raw-$SLURM_JOB_ID
mkdir -p "$out"
export executable=${FM_EXECUTABLE:?}
options=${FM_OPTIONS---line-search Backtracking --prediction LinearExtrapolation}
label=${FM_LABEL:-raw}
log "start, 100 grains (raw mesh), 5 steps, options $options," \
  "executable $executable, scripts $P"
export mesh=$MESHES/v100/100crystals.msh
export vectors=$MESHES/v100/vectors_100crystals.txt
export refinement=1 nsteps=5 solver_label=best
solver="--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance"
solver="$solver EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01"
export solver_gpu=$solver solver_cpu=$solver
export run_timeout_gpu=2400 run_timeout_cpu=${FM_TIMEOUT_CPU:-3600}
extra_options="--single-material $options" label=$label modes=cpu \
  "$P/cases/measure-polycrystal-grains.sh" "$out" > "$out/$label.out" 2>&1
log "$label: $(tail -1 "$out/series.txt")"
t=$out/$label-r1-cpu-n1-best.txt
log "$label: $(grep -c 'starts again from the default prediction' "$t")" \
  "extrapolated estimates discarded," \
  "$(grep -c 'integration failure' "$t") estimates moved after an integration failure," \
  "$(grep -c 'the resolution failed' "$t") failed resolutions"
energy "$out/energy.txt"
log end
