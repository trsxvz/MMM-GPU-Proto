#!/bin/bash
# One material for all the grains (option --single-material of the
# polycrystal: orientations given point by point, one batch of integration)
# against one material per grain, on one node. The executable is built first
# on the node by build-fix.sh in FM_FIX (default <scratch>/fix2: modified
# mfem-mgis and mm-opera-hpc in patched/, variant B only).
# - CPU, FM_CPU_GRAINS grains (default 100), 96 processes:
#   - the executable FM_EXECUTABLE_OLD (default: the one of the series, one
#     material per grain) and the new one with one material per grain give
#     identical outputs (rotation matrices kept, same values);
#   - one material for all the grains against one material per grain:
#     largest differences of the macroscopic values and of the averages per
#     grain (compare-outputs.py), the sums of the assembly being made in
#     another order;
# - GPU, FM_GPU_GRAINS grains (default 400), 4 MI300A: times of both.
# FM_NSTEPS steps of 1 s (default 3), best solver of each side.
#SBATCH --job-name=fm-single-ab
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
out=$R/single-ab-$SLURM_JOB_ID
mkdir -p "$out"
export FM_FIX=${FM_FIX:-$W/fix2}
log "start, build in $FM_FIX, scripts $P"
if ! FM_BUILD_VARIANTS=B bash "$FM_FIX/build-fix.sh" > "$out/build.log" 2>&1; then
  log "build failed (see $out/build.log)"
  exit 1
fi
new=$FM_FIX/bin/uniaxial-polycrystal-B
old=${FM_EXECUTABLE_OLD:-$W/bin/uniaxial-polycrystal-9767ea67}
log "built: $new"
export refinement=1 nsteps=${FM_NSTEPS:-3} solver_label=best
export solver_gpu=${FM_SOLVER_GPU:?} solver_cpu=${FM_SOLVER_CPU:?}
export run_timeout_gpu=2400 run_timeout_cpu=2400
compare() {
  # largest scaled differences of the macroscopic values (time, gradient,
  # stresses) and of the averages per grain (comment lines removed)
  local a=$out/$1 b=$out/$2 r
  r=$(python3 -I "$P/cases/compare-outputs.py" "$a.res" "$b.res" \
    --rtol 1e-6 --groups 0 1-3 4-6 2>&1 | tail -1)
  grep -v '^#' "$a/avgStressPolycrystal" > "$out/avg-$1.txt"
  grep -v '^#' "$b/avgStressPolycrystal" > "$out/avg-$2.txt"
  r="$r; per grain: $(python3 -I "$P/cases/compare-outputs.py" \
    "$out/avg-$1.txt" "$out/avg-$2.txt" --rtol 1e-6 2>&1 | tail -1)"
  echo "$r"
}
n=${FM_CPU_GRAINS:-100}
export mesh=$W/meshes/clean/v$n/${n}crystals.msh
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
executable=$old label=g${n}c-old modes=cpu \
  "$P/cases/measure-polycrystal-grains.sh" "$out"
log "cpu old: $(tail -1 "$out/series.txt")"
executable=$new label=g${n}c-pergrain modes=cpu \
  "$P/cases/measure-polycrystal-grains.sh" "$out"
log "cpu per grain: $(tail -1 "$out/series.txt")"
executable=$new extra_options=--single-material label=g${n}c-single modes=cpu \
  "$P/cases/measure-polycrystal-grains.sh" "$out"
log "cpu single: $(tail -1 "$out/series.txt")"
a=g${n}c-old-r1-cpu-n1-best
b=g${n}c-pergrain-r1-cpu-n1-best
c=g${n}c-single-r1-cpu-n1-best
if cmp -s "$out/$a.res" "$out/$b.res" &&
  diff -rq "$out/$a" "$out/$b" > "$out/diff-cpu.txt" 2>&1; then
  log "cpu: old and new executables, one material per grain: identical outputs"
else
  log "cpu: old and new executables differ (see $out/diff-cpu.txt)"
fi
log "cpu: single material against one per grain: $(compare "$b" "$c")"
n=${FM_GPU_GRAINS:-400}
export mesh=$W/meshes/clean/v$n/${n}crystals.msh
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
executable=$new label=g${n}c-pergrain modes=gpu \
  "$P/cases/measure-polycrystal-grains.sh" "$out"
log "gpu per grain: $(tail -1 "$out/series.txt")"
executable=$new extra_options=--single-material label=g${n}c-single modes=gpu \
  "$P/cases/measure-polycrystal-grains.sh" "$out"
log "gpu single: $(tail -1 "$out/series.txt")"
log end
