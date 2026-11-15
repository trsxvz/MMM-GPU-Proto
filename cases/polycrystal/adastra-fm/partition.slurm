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
# FM_PARTITION_WORKERS=<w>: the parts are written by <w> processes of
# meshes/partition/split-mesh instead of one process per part, for the
# numbers of parts whose processes cannot each read the serial mesh (6144
# parts of 38,400 grains). Each partition is then checked by split-mesh
# --check on all the cores, which writes the file `checked` in its
# directory when the parts are consistent: a partition is only kept by the
# next jobs with this file. FM_CHECK: partitions "grains:processes" already
# made, checked the same way (same fingerprint of the elements for all the
# partitions of a mesh). The job fails if a partitioning or a check fails. With FM_SPLIT_TEST=1, first compares split-mesh on
# FM_SPLIT_TEST_CASE (default 400:96, one node, files on its local disk) to
# the parts made by MFEM alone (split-mesh --mfem) and by partition-mesh,
# and the CPU computations reading the parts of the two tools
# (FM_SPLIT_TEST_STEPS time steps, default 2; 0: no computation).
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
log "start, partitions ${FM_PARTITIONS:-none}, tool in $S"
# environment of the prototype for the build and the partitioning only: it
# redefines variables of the chain (P)
# split-mesh against MFEM alone and against partition-mesh: same files
split_test() {
  local c=${FM_SPLIT_TEST_CASE:-400:96}
  local n=${c%%:*} np=${c##*:} t=/tmp/split-test-$SLURM_JOB_ID status=0
  local m=$W/meshes/clean/v$n/${n}crystals.msh
  mkdir -p "$t"/{split,few,mfem,old}
  srun -n "$np" --cpu-bind=cores "$S/split-mesh" "$m" 1 "$np" "$t/split/mesh." > "$t/split.log" 2>&1 || status=1
  srun -n 7 --cpu-bind=cores "$S/split-mesh" "$m" 1 "$np" "$t/few/mesh." > "$t/few.log" 2>&1 || status=1
  srun -n "$np" --cpu-bind=cores "$S/split-mesh" --mfem "$m" 1 "$t/mfem/mesh." > "$t/mfem.log" 2>&1 || status=1
  srun -n "$np" --cpu-bind=cores "$S/partition-mesh" "$m" 1 "$t/old/mesh." > "$t/old.log" 2>&1 || status=1
  log "split test v$n np$np: status $status, $(tail -1 "$t/split.log"); $(tail -2 "$t/split.log" | head -1)"
  for other in few mfem old; do
    local same=0 different=0
    for f in "$t"/split/mesh.*; do
      if cmp -s "$f" "$t/$other/$(basename "$f")"; then
        same=$((same + 1))
      else
        different=$((different + 1))
      fi
    done
    log "split test, $np processes of split-mesh against $other: $same files identical, $different different"
    # files different from the ones of partition-mesh: the computations
    # reading them are compared below
    if [ "$other" = old ]; then
      [ "$different" = 0 ] || touch "$t/compute"
      continue
    fi
    [ "$different" = 0 ] && [ "$same" = "$np" ] || status=1
  done
  for parts in split old; do
    srun -n "$SLURM_NTASKS" --cpu-bind=cores "$S/split-mesh" --check "$np" "$t/$parts/mesh." > "$t/check-$parts.log" 2>&1 || status=1
    log "split test, check of $parts: $(tail -3 "$t/check-$parts.log" | tr '\n' ' ')"
  done
  # a damaged part must be refused: a shared vertex replaced in one part
  mkdir -p "$t/bad"
  cp "$t"/split/mesh.* "$t/bad/"
  awk '/^shared_vertices [1-9]/ && !done {print; getline; print ($1 == 0 ? 1 : 0); done = 1; next} {print}' \
    "$t/split/mesh.000001" > "$t/bad/mesh.000001"
  if srun -n "$SLURM_NTASKS" --cpu-bind=cores "$S/split-mesh" --check "$np" "$t/bad/mesh." > "$t/check-bad.log" 2>&1; then
    log "split test: damaged part NOT detected"
    status=1
  else
    log "split test, damaged part: $(grep -h -e "^group" -e "check" "$t/check-bad.log" | tr '\n' ' ')"
  fi
  rm -rf "$t/bad" "$t/few" "$t/mfem"
  [ "$status" = 0 ] || log "split test FAILED (logs in $t on $(hostname))"
  return "$status"
}
# usage: check_partition grains processes
# reads back the parts of a partition on all the cores (split-mesh --check)
# and, if they are consistent, writes the file `checked` in its directory:
# only a checked partition is kept by the next jobs
check_partition() {
  local n=$1 np=$2 d=$W/meshes/part/v$1-np$2 start s
  start=$(date +%s)
  rm -f "$d/checked"
  srun -n "$SLURM_NTASKS" --cpu-bind=cores "$S/split-mesh" --check "$np" "$d/mesh." > "$d/check.log" 2>&1
  s=$?
  log "v$n np$np: check status $s, $(( $(date +%s) - start )) s, $(tail -3 "$d/check.log" | tr '\n' ' ')"
  [ "$s" = 0 ] && touch "$d/checked"
  return "$s"
}
partition() {
  set +u
  source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
  set -u
  if ! make -C "$S" partition-mesh split-mesh > "$S/build.log" 2>&1; then
    log "build of partition-mesh and split-mesh failed (see $S/build.log)"
    return 1
  fi
  if [ "${FM_SPLIT_TEST:-0}" = 1 ]; then
    split_test || return 1
  fi
  local status=0 s
  for p in ${FM_PARTITIONS:-}; do
    n=${p%%:*} np=${p##*:}
    d=$W/meshes/part/v$n-np$np
    # partitions already made and checked are kept
    if [ -e "$d/checked" ]; then
      log "v$n np$np: already partitioned and checked"
      continue
    fi
    mkdir -p "$d"
    start=$(date +%s)
    # FM_PARTITION_TASKS_PER_NODE: fewer processes per node for the large
    # meshes (each process reads the whole serial mesh, about 9 GB for
    # 12,800 grains)
    if [ -n "${FM_PARTITION_WORKERS:-}" ]; then
      srun -n "$FM_PARTITION_WORKERS" ${FM_PARTITION_TASKS_PER_NODE:+--ntasks-per-node=$FM_PARTITION_TASKS_PER_NODE} \
        --cpu-bind=cores "$S/split-mesh" \
          "$W/meshes/clean/v$n/${n}crystals.msh" 1 "$np" "$d/mesh." > "$d/partition.log" 2>&1
      s=$?
      log "v$n np$np: status $s, $(( $(date +%s) - start )) s, $FM_PARTITION_WORKERS processes of split-mesh, $(tail -2 "$d/partition.log" | tr '\n' ' '), $(du -sh "$d" | cut -f1)"
    else
      srun -n "$np" ${FM_PARTITION_TASKS_PER_NODE:+--ntasks-per-node=$FM_PARTITION_TASKS_PER_NODE} \
        --cpu-bind=cores "$S/partition-mesh" \
          "$W/meshes/clean/v$n/${n}crystals.msh" 1 "$d/mesh." > "$d/partition.log" 2>&1
      s=$?
      log "v$n np$np: status $s, $(( $(date +%s) - start )) s, $(tail -1 "$d/partition.log"), $(du -sh "$d" | cut -f1)"
    fi
    # a partition whose tool failed is not checked, and the job fails
    if [ "$s" != 0 ]; then
      status=1
      continue
    fi
    check_partition "$n" "$np" || status=1
  done
  # partitions already made, read back
  for p in ${FM_CHECK:-}; do
    check_partition "${p%%:*}" "${p##*:}" || status=1
  done
  [ "$status" = 0 ] || log "partitions: FAILURES"
  return "$status"
}
(partition) || exit 1
if [ "${FM_SPLIT_TEST:-0}" = 1 ]; then
  c=${FM_SPLIT_TEST_CASE:-400:96}
  n=${c%%:*} t=/tmp/split-test-$SLURM_JOB_ID out=$R/split-test-$SLURM_JOB_ID
  if [ -e "$t/compute" ] && [ "${FM_SPLIT_TEST_STEPS:-2}" != 0 ]; then
    # CPU computations reading the parts of partition-mesh and of split-mesh
    mkdir -p "$out"
    for parts in old split; do
      mesh=$t/$parts/mesh. vectors=$MESHES/v$n/vectors_${n}crystals.txt \
        refinement=0 label=g${n}c-$parts nsteps=${FM_SPLIT_TEST_STEPS:-2} \
        duration=${FM_SPLIT_TEST_STEPS:-2} solver_label=best \
        solver_cpu=${FM_SOLVER_CPU:?} run_timeout_cpu=2400 modes=cpu \
        extra_options="$extra_options --mesh_mode Restart" \
        "$P/cases/measure-polycrystal-grains.sh" "$out"
      log "split test, computation on the parts of $parts: $(tail -1 "$out/series.txt")"
    done
    log "split test, .res: $(tail -1 "$out/g${n}c-old-r0-cpu-n1-best.res") / $(tail -1 "$out/g${n}c-split-r0-cpu-n1-best.res")"
  elif [ ! -e "$t/compute" ]; then
    log "split test: files identical to the ones of partition-mesh, no computation needed"
  fi
  rm -rf "$t"
fi
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
