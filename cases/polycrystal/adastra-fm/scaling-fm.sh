#!/bin/bash
# Weak and strong scalability of the FM (15/11/2026), 64 nodes at most,
# 600 grains per node at the base load (4.94 M unknowns per GPU, 206,000
# per core). Run on the login node, in the directory holding fm-chain/:
#   fm-chain/scaling-fm.sh meshes
#   fm-chain/scaling-fm.sh computations
# - Weak: 600 N grains on N nodes, N = 1 to 64 (GPU x2, CPU x1; the CPU
#   point of 38,400 grains on 64 nodes is the job 5989775).
# - Strong: 600 grains from 1 to 64 nodes, 2400 from 4 and 9600 from 16,
#   doubling the nodes (GPU x1, CPU x1), the base point being the one of
#   the weak series.
# The parts of the CPU of the strong series (33,000 files) are removed
# after each computation, and so are the 3072 parts of 19,200 grains: the
# quota of files of the group has no room for a set of 6144 parts next to
# them. The points on 64 nodes come last: their jobs wait the longest, and
# nothing waits for the files of their parts.
set -eu
export FM_RUN=2026-11-15-scaling
scaling=fm-chain/scaling.sh
# usage: after_mesh grains...
# dependency on the mesh jobs of these numbers of grains still in the queue
after_mesh() {
  local names ids
  names=$(printf 'fm-scal-mesh-g%s,' "$@")
  ids=$(squeue -h -u "$USER" -n "${names%,}" -o %i | paste -sd:)
  echo "${ids:+afterok:$ids}"
}
case ${1:?meshes or computations} in
  meshes)
    $scaling mesh 600 1200 2400 4800 9600 19200
    ;;
  computations)
    # GPU: weak series, twice each
    for n in 600 1200 2400 4800 9600 19200; do
      FM_AFTER=$(after_mesh $n) FM_GPU_REPS=2 $scaling gpu $n $((n / 600))
    done
    FM_AFTER=afterok:$($scaling check 38400:256) FM_GPU_REPS=2 \
      $scaling gpu 38400 64
    # GPU: strong series
    FM_AFTER=$(after_mesh 600) $scaling gpu 600 2 4 8 16 32 64
    FM_AFTER=$(after_mesh 2400) $scaling gpu 2400 8 16 32 64
    FM_AFTER=$(after_mesh 9600) $scaling gpu 9600 32 64
    # CPU: weak series, parts kept
    for n in 600 1200 2400 4800 9600; do
      FM_AFTER=$(after_mesh $n) $scaling cpu $n:$((n / 600)):keep
    done
    # CPU: 19,200 grains of the weak series and strong series up to 32
    # nodes, then the points on 64 nodes
    large=$(FM_AFTER=$(after_mesh 19200) $scaling cpu 19200:32)
    strong=$(FM_AFTER=$(after_mesh 600 2400 9600) FM_TIME_CPU=01:30:00 \
      FM_TIMEOUT_CPU=4500 $scaling cpu 9600:32 2400:32 600:32 2400:16 \
      600:16 2400:8 600:8 600:4 600:2)
    FM_AFTER=afterany:$large,afterany:$strong FM_TIME_CPU=01:30:00 \
      FM_TIMEOUT_CPU=3000 $scaling cpu 600:64 2400:64 9600:64
    ;;
esac
