#!/bin/bash
# Mean stresses of the polycrystal computed with the weights of the points
# cached (mfem-mgis), against the current computation, on one node. The
# executables are built first on the node by build-fix.sh in FM_FIX
# (default <scratch>/fix, holding build-fix.sh and the modified files in
# patched/): A (current mfem-mgis) and B (modified one), both timing the
# macroscopic stresses computed after each resolution
# (computeMacroscopicCauchyStress, mm-opera-hpc).
# - CPU, FM_CPU_GRAINS grains (default 100), 96 processes: the computations
#   are deterministic, the outputs of A and B must be identical
#   (<computation>.res and the averages per grain);
# - GPU, FM_GPU_GRAINS grains (default 400), 4 MI300A: time of the mean
#   stresses and wall clock time of A and B.
# FM_NSTEPS steps of 1 s (default 3), best solver of each side.
#SBATCH --job-name=fm-fix-ab
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
out=$R/fix-ab-$SLURM_JOB_ID
mkdir -p "$out"
F=${FM_FIX:-$W/fix}
log "start, build in $F, scripts $P"
if ! bash "$F/build-fix.sh" > "$out/build.log" 2>&1; then
  log "build failed (see $out/build.log)"
  exit 1
fi
log "built: $(grep -E '^(A|B) built' "$out/build.log" | tr '\n' ' ')"
FM_EXECUTABLE_A=$F/bin/uniaxial-polycrystal-A
FM_EXECUTABLE_B=$F/bin/uniaxial-polycrystal-B
export refinement=1 nsteps=${FM_NSTEPS:-3} solver_label=best
export solver_gpu=${FM_SOLVER_GPU:?} solver_cpu=${FM_SOLVER_CPU:?}
export run_timeout_gpu=2400 run_timeout_cpu=2400
n=${FM_CPU_GRAINS:-100}
export mesh=$W/meshes/clean/v$n/${n}crystals.msh
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
for v in A B; do
  e=FM_EXECUTABLE_$v
  executable=${!e} label=g${n}c-$v modes=cpu \
    "$P/cases/measure-polycrystal-grains.sh" "$out"
  log "cpu $v: $(tail -1 "$out/series.txt")"
done
a=$out/g${n}c-A-r1-cpu-n1-best
b=$out/g${n}c-B-r1-cpu-n1-best
if cmp -s "$a.res" "$b.res" && diff -rq "$a" "$b" > "$out/diff-cpu.txt" 2>&1; then
  log "cpu: outputs of A and B identical"
else
  log "cpu: outputs of A and B differ (see $out/diff-cpu.txt)"
fi
n=${FM_GPU_GRAINS:-400}
export mesh=$W/meshes/clean/v$n/${n}crystals.msh
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
for v in A B; do
  e=FM_EXECUTABLE_$v
  executable=${!e} label=g${n}c-$v modes=gpu \
    "$P/cases/measure-polycrystal-grains.sh" "$out"
  log "gpu $v: $(tail -1 "$out/series.txt")"
done
log end
