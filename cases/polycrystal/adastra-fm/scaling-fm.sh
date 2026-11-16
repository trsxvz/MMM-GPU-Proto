#!/bin/bash
# Weak and strong scalability of the FM (15/11/2026), 64 nodes at most,
# 600 grains per node at the base load (4.94 M unknowns per GPU, 206,000
# per core). Run on the login node, in the directory holding fm-chain/:
#   fm-chain/scaling-fm.sh meshes
#   fm-chain/scaling-fm.sh computations
# - Weak: 600 N grains on N nodes, N = 1 to 64.
# - Strong: 2400 grains from 4 to 64 nodes, 4800 from 8, 9600 from 16 and
#   19,200 from 32, doubling the nodes, each one from the fewest nodes the
#   memory of the GPUs allows (the point of the weak series).
# - GPU x3, CPU x1 per point, with the executable writing the iterations
#   of the linear solver.
# All the points of a number of nodes are in one job (two on 64 nodes): one
# wait in the queue each.
# Quota of files of the group (3000 files left at the start): each job
# reads the kept parts of its point of the weak series first and removes
# them, which leaves room for the parts of its points of the strong series,
# written in the job and removed after each computation. The 6144 parts of
# the strong series on 64 nodes need the room freed by all the other jobs.
# The parts of the weak series are written again at the end.
# (A first study, with an executable which did not write these iterations
# and one job per point, also had 600 grains from 1 to 64 nodes:
# fm/2026-11-15-scaling.)
set -eu
export FM_RUN=${FM_RUN:-2026-11-15-scaling-cg}
scaling=fm-chain/scaling.sh
# usage: after_mesh grains...
# dependency on the mesh jobs of these numbers of grains still in the queue
after_mesh() {
  local names ids
  names=$(printf 'fm-scal-mesh-g%s,' "$@")
  ids=$(squeue -h -u "$USER" -n "${names%,}" -o %i | paste -sd:)
  echo "${ids:+afterok:$ids}"
}
# usage: points suffix grains...
points() {
  local suffix=$1
  shift
  printf "%s$suffix " "$@"
}
case ${1:?meshes or computations} in
  meshes)
    $scaling mesh 600 1200 2400 4800 9600 19200
    ;;
  computations)
    # RVEs of the strong series on each number of nodes, the largest first
    declare -A strong=([8]="2400" [16]="4800 2400" [32]="9600 4800 2400"
      [64]="19200 9600 4800 2400")
    declare -A time=([8]=02:00:00 [16]=02:30:00 [32]=03:00:00)
    kept=
    others=
    for nodes in 1 2 4 8 16 32; do
      weak=$((600 * nodes))
      rves="$weak ${strong[$nodes]:-}"
      job=$(FM_AFTER=$(after_mesh $rves) FM_GPU_POINTS=$(points :3 $rves) \
        FM_CPU_POINTS=$rves FM_TIME_CPU=${time[$nodes]:-01:30:00} \
        FM_TIMEOUT_GPU=600 FM_TIMEOUT_CPU=3600 $scaling nodes $nodes)
      kept+="$weak:$((96 * nodes)) "
      others+=":$job"
    done
    # 64 nodes: the GPUs and the CPU of the weak series, whose parts exist
    checked=$($scaling check 38400:6144)
    FM_AFTER=afterok:$checked FM_GPU_POINTS=$(points :3 38400 ${strong[64]}) \
      FM_CPU_POINTS=38400:keep FM_TIME_CPU=02:00:00 FM_TIMEOUT_GPU=600 \
      FM_TIMEOUT_CPU=4800 $scaling nodes 64
    # 64 nodes: the CPU of the strong series, the first parts written before
    first=$(FM_AFTER=afterany$others $scaling parts 19200:6144)
    last=$(FM_AFTER=afterany:$first FM_CPU_POINTS=${strong[64]} \
      FM_TIME_CPU=02:00:00 FM_TIMEOUT_CPU=2400 $scaling nodes 64)
    FM_AFTER=afterany:$last $scaling parts $kept
    ;;
esac
