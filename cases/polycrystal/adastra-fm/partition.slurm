#!/bin/bash
# Partitioned meshes of the polycrystal (meshes/partition/partition-mesh,
# built first on the node from FM_PARTITION_SOURCES, default
# <scratch>/partition-tool): for each "grains:processes" of FM_PARTITIONS,
# the cleaned mesh of <grains> grains refined once and split in
# <processes> parts, in <scratch>/meshes/part/v<grains>-np<processes>/mesh.<rank>.
# With FM_TEST=1, then compares on the CPU (100 grains, 96 processes) the
# computation reading the partitioned mesh ("Restart", no refinement) with
# the one building it ("FromScratch", one refinement), and the reading of
# the mesh on the GPU (400 grains, 4 processes).
#SBATCH --job-name=fm-partition
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --time=01:30:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
S=${FM_PARTITION_SOURCES:-$W/partition-tool}
log "start, partitions ${FM_PARTITIONS:?}, tool in $S"
# environment of the prototype for the build and the partitioning only: it
# redefines variables of the chain (P)
partition() {
  set +u
  source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
  set -u
  if ! make -C "$S" partition-mesh > "$S/build.log" 2>&1; then
    log "build of partition-mesh failed (see $S/build.log)"
    return 1
  fi
  for p in $FM_PARTITIONS; do
    n=${p%%:*} np=${p##*:}
    d=$W/meshes/part/v$n-np$np
    # partitions already made are kept
    if grep -qs " parts, " "$d/partition.log"; then
      log "v$n np$np: already partitioned"
      continue
    fi
    mkdir -p "$d"
    start=$(date +%s)
    # FM_PARTITION_TASKS_PER_NODE: fewer processes per node for the large
    # meshes (each process reads the whole serial mesh, about 9 GB for
    # 12,800 grains)
    srun -n "$np" ${FM_PARTITION_TASKS_PER_NODE:+--ntasks-per-node=$FM_PARTITION_TASKS_PER_NODE} \
      --cpu-bind=cores "$S/partition-mesh" \
        "$W/meshes/clean/v$n/${n}crystals.msh" 1 "$d/mesh." > "$d/partition.log" 2>&1
    log "v$n np$np: status $?, $(( $(date +%s) - start )) s, $(tail -1 "$d/partition.log"), $(du -sh "$d" | cut -f1)"
  done
}
(partition) || exit 1
if [ "${FM_TEST:-0}" = 1 ]; then
  out=$R/partition-test-$SLURM_JOB_ID
  mkdir -p "$out"
  export nsteps=3 duration=3 solver_label=best
  export solver_gpu=${FM_SOLVER_GPU:?} solver_cpu=${FM_SOLVER_CPU:?}
  export run_timeout_gpu=2400 run_timeout_cpu=2400
  section() {
    awk -F'|' -v s="$2" '$3 ~ ("--> " s) {gsub(/ /, "", $7); print $7; exit}' "$1"
  }
  base_options=$extra_options
  for mode in cpu gpu; do
    if [ $mode = cpu ]; then n=100 np=96; else n=400 np=4; fi
    export vectors=$MESHES/v$n/vectors_${n}crystals.txt
    mesh=$W/meshes/clean/v$n/${n}crystals.msh refinement=1 label=g${n}c-scratch \
      modes=$mode "$P/cases/measure-polycrystal-grains.sh" "$out"
    log "$mode scratch: $(tail -1 "$out/series.txt"); Mesh::Constructor $(section "$out/g${n}c-scratch-r1-$mode-n1-best.txt" Mesh::Constructor) s"
    mesh=$W/meshes/part/v$n-np$np/mesh. refinement=0 label=g${n}c-restart \
      extra_options="$base_options --mesh_mode Restart" modes=$mode \
      "$P/cases/measure-polycrystal-grains.sh" "$out"
    log "$mode restart: $(tail -1 "$out/series.txt"); Mesh::Constructor $(section "$out/g${n}c-restart-r0-$mode-n1-best.txt" Mesh::Constructor) s"
    if cmp -s "$out/g${n}c-scratch-r1-$mode-n1-best.res" "$out/g${n}c-restart-r0-$mode-n1-best.res"; then
      log "$mode: .res identical"
    else
      log "$mode: .res differ: $(tail -1 "$out/g${n}c-scratch-r1-$mode-n1-best.res") / $(tail -1 "$out/g${n}c-restart-r0-$mode-n1-best.res")"
    fi
  done
fi
log end
