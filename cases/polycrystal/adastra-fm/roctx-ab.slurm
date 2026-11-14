#!/bin/bash
# Cost of the ranges of roctx without profiler: the final computation of
# the polycrystal (FM_GRAINS grains, mesh partitioned beforehand, one
# material for all the grains, line search, prediction by linear
# extrapolation, FM_NSTEPS steps over 20 s) on all the GPUs of the
# allocation, with the executables FM_EXECUTABLES (`<name>=<path>` separated
# by spaces, for instance the executable of the measurements and the ones
# built without and with the option enable-roctx of mfem-mgis), run one
# after the other at each of the FM_GPU_REPS repetitions, on the same nodes.
# A line per computation in <results>/series.txt (see
# cases/measure-polycrystal-grains.sh).
#SBATCH --job-name=fm-roctx-ab
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
n=${FM_GRAINS:?}
np=$((SLURM_JOB_NUM_NODES * 4))
out=$R/roctx-ab-g$n-$SLURM_JOB_ID
mkdir -p "$out"
export mesh=${FM_MESH:-$W/meshes/part/v$n-np$np/mesh.}
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
export refinement=0 nsteps=${FM_NSTEPS:-20} duration=20 solver_label=best
export modes=gpu run_timeout_gpu=${FM_TIMEOUT_GPU:-2400}
export solver_gpu=${FM_SOLVER_GPU:---non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01}
export extra_options=${FM_EXTRA_OPTIONS:---single-material --line-search Backtracking --prediction LinearExtrapolation --mesh_mode Restart}
log "start, $n grains, $SLURM_JOB_NUM_NODES nodes, mesh $mesh, $nsteps steps," \
  "GPU x${FM_GPU_REPS:-2}, executables ${FM_EXECUTABLES:?}, options" \
  "$extra_options, solver $solver_gpu, scripts $P"
for i in $(seq "${FM_GPU_REPS:-2}"); do
  for e in $FM_EXECUTABLES; do
    executable=${e#*=} label=g${n}c-${e%%=*}-rep$i \
      "$P/cases/measure-polycrystal-grains.sh" "$out"
    log "${e%%=*} rep $i: $(tail -1 "$out/series.txt")"
  done
done
energy "$out/energy.txt"
log end
