#!/bin/bash
# Profile of the host of the final computation: polycrystal of FM_GRAINS
# grains on all the GPUs of the allocation (one process per GPU), final
# configuration (see fm-chain/profile-trace.slurm), the processes
# FM_TRACE_RANKS running under Linux perf, which samples their call stacks
# (fm-chain/trace-ranks.sh). It tells which functions of the host run while
# the GPU waits, for instance during the first assembly of the tangent
# operator, which the traces of the GPU do not show.
#
#   sbatch -A genden15 --nodes=1 --time=01:30:00 \
#     --export=ALL,FM_RUN=<name>,FM_PROTO=<scripts>,FM_GRAINS=600 \
#     fm-chain/profile-host.slurm
#
# Variables:
# - FM_GRAINS (required): number of grains;
# - FM_MESH, FM_VECTORS, FM_EXECUTABLE, FM_SOLVER_GPU, FM_EXTRA_OPTIONS,
#   FM_TIMEOUT_GPU: as in fm-chain/profile-trace.slurm (FM_EXECUTABLE is
#   required);
# - FM_NSTEPS, FM_DURATION: steps and simulated time in s (default 1 and
#   1: the setup and the first step);
# - FM_TRACE_RANKS: sampled processes (default `0`);
# - FM_TRACE_OPTIONS: options of `perf record` (see
#   fm-chain/trace-ranks.sh).
#
# Results in <scratch>/fm/<FM_RUN>/profile-host-g<grains>-<job>:
# - trace/rank<rank>-<host>/perf.data: samples of each sampled process
#   (250 MB for a computation of one minute), to read with `perf report`;
# - report-rank<rank>.txt, written by the command of `perf report` given
#   in the journal of the chain, to run on the login node: for each
#   process, share of the samples of each function with its callees and
#   without them, and number of its samples (above 0.3 % of the samples of
#   the process, all threads included);
# - series.txt: one line per computation
#   (cases/measure-polycrystal-grains.sh).
#SBATCH --job-name=fm-profile-host
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
out=$R/profile-host-g$n-$SLURM_JOB_ID
mkdir -p "$out"
export mesh=${FM_MESH:-$W/meshes/part/v$n-np$np/mesh.}
export vectors=${FM_VECTORS:-$MESHES/v$n/vectors_${n}crystals.txt}
export refinement=0 nsteps=${FM_NSTEPS:-1}
export duration=${FM_DURATION:-1}
export modes=gpu run_timeout_gpu=${FM_TIMEOUT_GPU:-2400}
export solver_gpu=${FM_SOLVER_GPU:---non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01}
export executable=${FM_EXECUTABLE:?}
export extra_options=${FM_EXTRA_OPTIONS:---single-material --line-search Backtracking --prediction LinearExtrapolation --mesh_mode Restart}
export FM_TRACE_DIR=$out/trace FM_TRACE_RANKS=${FM_TRACE_RANKS:-0}
export FM_TRACE_TOOL=perf
log "start, $n grains, $SLURM_JOB_NUM_NODES nodes, $np processes, mesh $mesh," \
  "$nsteps steps over $duration s, executable $executable," \
  "options $extra_options, solver $solver_gpu, sampled processes" \
  "$FM_TRACE_RANKS, scripts $P"
if [ ! -f "${mesh}000000" ] || [ ! -x "$executable" ]; then
  log "mesh $mesh or executable missing"
  exit 1
fi
profiler="$SLURM_SUBMIT_DIR/fm-chain/trace-ranks.sh" label=g${n}c \
  solver_label=perf "$P/cases/measure-polycrystal-grains.sh" "$out"
log "sampled: $(tail -1 "$out/series.txt")"
# the report is left to the login node, where it takes a few seconds: on a
# node, `perf report` did not end in 20 minutes (job 5996554)
for d in "$FM_TRACE_DIR"/rank*; do
  [ -f "$d/perf.data" ] || continue
  rank=${d##*/rank}
  log "$(du -h "$d/perf.data" | cut -f1) of samples; report, on the login node:" \
    "perf report -i $d/perf.data --stdio --children -n -g none" \
    "--percent-limit 0.3 --sort dso,symbol | grep -v '^\$'" \
    "> $out/report-rank${rank%%-*}.txt"
done
log end
