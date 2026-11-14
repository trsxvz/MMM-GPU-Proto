#!/bin/bash
# Trace of the final computation, in the style of Nsight Systems: polycrystal
# of FM_GRAINS grains on all the GPUs of the allocation (one process per
# GPU), final configuration (one material for all the grains, line search,
# prediction by linear extrapolation, NewtonK0 with the first choice of
# Eisenstat and Walker from 0.01, mesh partitioned beforehand and read with
# --mesh_mode Restart), under rocprofv3 with the timed sections, the time
# steps and the resolutions declared as ranges of roctx (executable built
# with the option enable-roctx of mfem-mgis, fm-chain/build-roctx.sh).
#
# The number of nodes and the time limit are given at the submission:
#   sbatch -A genden15 --nodes=1 --time=01:30:00 \
#     --export=ALL,FM_RUN=<name>,FM_PROTO=<scripts>,FM_GRAINS=400 \
#     fm-chain/profile-trace.slurm
#
# Variables:
# - FM_GRAINS (required): number of grains;
# - FM_MESH: partitioned mesh (default
#   <scratch>/meshes/part/v<grains>-np<processes>/mesh.), FM_VECTORS:
#   orientations (default the ones of the series);
# - FM_NSTEPS, FM_DURATION: steps and simulated time in s (default 20 and
#   20);
# - FM_EXECUTABLE: executable (default
#   <scratch>/fix-roctx/bin/uniaxial-polycrystal-on);
# - FM_SOLVER_GPU, FM_EXTRA_OPTIONS: options of the solvers and of the
#   polycrystal (default: the final configuration);
# - FM_TRACE_RANKS: traced processes (default `node:0`, the processes of
#   the first node; `all`, `node:<i>,<j>` or ranks `0-3,128`, see
#   fm-chain/trace-ranks.sh): the trace of all the processes of a large
#   computation would be too large to be read;
# - FM_TRACE_OPTIONS: options of rocprofv3 (default: --runtime-trace
#   --marker-trace --output-format pftrace csv; see fm-chain/trace-ranks.sh
#   for the buffer of the trace of Perfetto);
# - FM_ANALYSIS: 1 to analyse the traces in the job, 0 to only write the
#   commands of the analysis in the journal of the chain, for the login
#   node (default: 1 up to 4 nodes, 0 beyond: the analysis takes a few
#   minutes on one core, during which a large allocation would wait);
# - FM_REFERENCE=1: the same computation without profiler first, on the
#   same nodes, to measure the overhead of the trace (default 0);
# - FM_TIMEOUT_GPU: time limit of a computation in s (default 2400).
#
# Results in <scratch>/fm/<FM_RUN>/profile-trace-g<grains>-<job>:
# - trace/rank<rank>-<host>/trace_results.pftrace: trace of each traced
#   process, to open in https://ui.perfetto.dev;
# - trace/rank<rank>-<host>/trace_*.csv: the same events as tables;
# - trace-zones.txt: time and occupancy of the GPU per zone
#   (fm-chain/trace-zones.py); api-gaps.txt, kernels.txt: idle time of the
#   GPUs and kernels (fm-chain/api-gaps.py, fm-chain/kernels.py);
# - series.txt: one line per computation
#   (cases/measure-polycrystal-grains.sh).
#SBATCH --job-name=fm-profile-trace
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
out=$R/profile-trace-g$n-$SLURM_JOB_ID
mkdir -p "$out"
export mesh=${FM_MESH:-$W/meshes/part/v$n-np$np/mesh.}
export vectors=${FM_VECTORS:-$MESHES/v$n/vectors_${n}crystals.txt}
export refinement=0 nsteps=${FM_NSTEPS:-20}
export duration=${FM_DURATION:-20}
export modes=gpu run_timeout_gpu=${FM_TIMEOUT_GPU:-2400}
export solver_gpu=${FM_SOLVER_GPU:---non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01}
export executable=${FM_EXECUTABLE:-$W/fix-roctx/bin/uniaxial-polycrystal-on}
export extra_options=${FM_EXTRA_OPTIONS:---single-material --line-search Backtracking --prediction LinearExtrapolation --mesh_mode Restart}
export FM_TRACE_DIR=$out/trace FM_TRACE_RANKS=${FM_TRACE_RANKS:-node:0}
log "start, $n grains, $SLURM_JOB_NUM_NODES nodes, $np processes, mesh $mesh," \
  "$nsteps steps over $duration s, executable $executable," \
  "options $extra_options, solver $solver_gpu, traced processes" \
  "$FM_TRACE_RANKS, scripts $P"
if { [ ! -f "${mesh}000000" ] && [ ! -f "$mesh" ]; } || [ ! -x "$executable" ]; then
  log "mesh $mesh or executable missing"
  exit 1
fi
if [ "${FM_REFERENCE:-0}" = 1 ]; then
  label=g${n}c solver_label=plain \
    "$P/cases/measure-polycrystal-grains.sh" "$out"
  log "without profiler: $(tail -1 "$out/series.txt")"
fi
profiler="$SLURM_SUBMIT_DIR/fm-chain/trace-ranks.sh" label=g${n}c \
  solver_label=trace "$P/cases/measure-polycrystal-grains.sh" "$out"
log "traced: $(tail -1 "$out/series.txt")"
log "traces: $(find "$FM_TRACE_DIR" -name '*.pftrace' | wc -l) processes," \
  "$(du -sh "$FM_TRACE_DIR" 2> /dev/null | cut -f1)," \
  "pftrace $(find "$FM_TRACE_DIR" -name '*.pftrace' -exec du -ch {} + 2> /dev/null | tail -1 | cut -f1)"
energy "$out/energy.txt"
zones="python3 -I fm-chain/trace-zones.py $FM_TRACE_DIR --steps > $out/trace-zones.txt"
analysis=${FM_ANALYSIS:-$((SLURM_JOB_NUM_NODES <= 4))}
if [ "$analysis" != 1 ]; then
  log "analysis to run on the login node, from $SLURM_SUBMIT_DIR: $zones"
  log end
  exit 0
fi
(cd "$SLURM_SUBMIT_DIR" && eval "$zones" 2>&1)
log "zones: $(grep -E '^mean of|^steady regime' "$out/trace-zones.txt" | tr '\n' ';')"
python3 -I "$SLURM_SUBMIT_DIR/fm-chain/api-gaps.py" "$FM_TRACE_DIR" \
  > "$out/api-gaps.txt" 2>&1
newton=$(tail -1 "$out/series.txt" | awk '{print $11}')
python3 -I "$SLURM_SUBMIT_DIR/fm-chain/kernels.py" "$FM_TRACE_DIR" "$newton" \
  > "$out/kernels.txt" 2>&1
log end
