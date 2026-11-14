#!/bin/bash
# Launcher placed before the executable of a MPI computation (variable
# `profiler` of cases/measure-polycrystal-grains.sh): the selected
# processes run under rocprofv3, the profiler of ROCm, the others run
# directly, without any overhead. A trace of all the processes of a large
# computation would be too large to be read.
#
# usage: trace-ranks.sh <executable> [arguments...]
#
# Variables read:
# - FM_TRACE_DIR (required): directory of the traces, one subdirectory
#   rank<rank>-<host> per traced process;
# - FM_TRACE_RANKS: traced processes (default `node:0`):
#   - `all`;
#   - `node:<i>[,<j>...]`: the processes of the nodes of these indices in
#     the allocation (SLURM_NODEID), for instance the 4 processes of the
#     first node;
#   - `<a>[-<b>][,...]`: ranks and ranges of ranks (SLURM_PROCID), for
#     instance `0-3,128`;
# - FM_TRACE_TOOL: `rocprofv3` (default) or `rocprof-sys`
#   (rocprof-sys-sample, formerly Omnitrace: the same events of the GPU,
#   with the sampling of the call stacks of the host and the calls of MPI,
#   in one Perfetto trace per process);
# - FM_TRACE_OPTIONS: options of the profiler (default for rocprofv3:
#   `--runtime-trace --marker-trace --output-format pftrace csv`, that is
#   the calls of the HIP runtime, the ranges of roctx, the kernels and the
#   memory copies, for Perfetto and as tables; default for rocprof-sys:
#   `--trace --freq 50 --include mpip rocprofiler-sdk`). The trace of
#   Perfetto is limited by its buffer, 1 GB per process by default, beyond
#   which the events are discarded (20 time steps of the polycrystal take
#   0.3 GB; the tables are not limited). `--perfetto-buffer-size <KB>`
#   changes it, below 4 GB: with 4194304, Perfetto refuses the buffer
#   ("Invalid buffer sizes"), no trace is written and the program hangs
#   at its end (ROCm 6.4.3);
# - FM_ROCPROFV3: profiler (default: the rocprofv3 or the
#   rocprof-sys-sample of ROCM_PATH).
set -u
rank=${SLURM_PROCID:-0}
node=${SLURM_NODEID:-0}
selection=${FM_TRACE_RANKS:-node:0}
traced=0
case "$selection" in
  all) traced=1 ;;
  node:*)
    selection=${selection#node:}
    for i in ${selection//,/ }; do
      [ "$i" = "$node" ] && traced=1
    done
    ;;
  *)
    for i in ${selection//,/ }; do
      if [ "$rank" -ge "${i%-*}" ] && [ "$rank" -le "${i#*-}" ]; then
        traced=1
      fi
    done
    ;;
esac 2> /dev/null
if [ $traced = 0 ]; then
  exec "$@"
fi
directory=${FM_TRACE_DIR:?}/rank$(printf %05d "$rank")-$(hostname -s)
mkdir -p "$directory"
rocm=${ROCM_PATH:-/opt/rocm}
case "${FM_TRACE_TOOL:-rocprofv3}" in
  rocprofv3)
    options=${FM_TRACE_OPTIONS:---runtime-trace --marker-trace --output-format pftrace csv}
    # shellcheck disable=SC2086
    exec "${FM_ROCPROFV3:-$rocm/bin/rocprofv3}" $options \
      -d "$directory" -o trace -- "$@"
    ;;
  rocprof-sys)
    options=${FM_TRACE_OPTIONS:---trace --freq 50 --include mpip rocprofiler-sdk}
    # shellcheck disable=SC2086
    exec "${FM_ROCPROFV3:-$rocm/bin/rocprof-sys-sample}" $options \
      --output "$directory" -- "$@"
    ;;
  *)
    echo "trace-ranks.sh: unknown profiler ${FM_TRACE_TOOL}" >&2
    exit 1
    ;;
esac
