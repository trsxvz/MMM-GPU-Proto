#!/bin/bash
# ParaView export of the polycrystal (mfem-mgis patches 0302 to 0305,
# mm-opera-hpc patches 0054 and 0055, JOURNAL.md, 2026-11-18): FM_GRAINS
# grains (default 600) on the GPUs of the allocation (FM_MODE=gpu, default)
# or on its cores (FM_MODE=cpu), mesh partitioned for the processes,
# FM_NSTEPS steps of 1 s (default 20), options of the measurements of the
# FM, element matrices added in a fixed order on the GPUs
# (MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered), displacement and von Mises
# stress exported at each step (paraview_export=yes,
# --enable-export-von-Mises-stress), three computations on the same nodes:
# - a: FM_EXECUTABLE_A (before the patches), two collections
#   (OutputFile-Uniaxial-polycrystal, PolycrystalParaviewOutput);
# - b: FM_EXECUTABLE_B, one collection with both fields
#   (PolycrystalParaviewResults);
# - b32: the same in single precision (--enable-paraview-single-precision).
# The computations on the cores are reproducible bit for bit, not the ones
# on the GPUs (hypre on the GPUs, see JOURNAL.md, 2026-11-18): on the
# GPUs, the fields of the three computations differ by the rounding of the
# computation, and FM_TOLERANCE (relative, default 1e-6) is used for the
# comparisons of the fields of b to a (split) and to the reference, the
# comparison in single precision then holding for the coordinates only.
# Checks (compare-vtu.py, the job fails otherwise):
# - the mesh and the fields of b are the ones of a, bit for bit (split);
# - the reals of b32 are the roundings to the nearest single precision
#   values of the ones of b, the other arrays being the same (float);
# - the macroscopic stresses (.res) and the averages per grain of the three
#   computations are the same, bit for bit.
# With FM_REFERENCE_EXPORT (directory of a computation writing the two
# collections, for instance with the atomic assembly, which is not
# reproducible), its fields are compared to the ones of b within the
# relative tolerance FM_TOLERANCE (default 1e-6), the largest differences
# being written. Reported in the journal: size and number of files of each
# export, times of the post-processings and of the writing of the files
# (sections ParaviewExportResults::Save, PolycrystalParaviewOutput::Save of
# B; NLEP::executePostProcessings of all). The exported files are then
# removed (quota of files), except the last cycle of b and of b32, kept for a
# check in ParaView, unless FM_KEEP_EXPORT=yes (all kept).
#SBATCH --job-name=fm-export-ab
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --time=01:30:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
A=${FM_EXECUTABLE_A:?}
B=${FM_EXECUTABLE_B:?}
n=${FM_GRAINS:-600}
mode=${FM_MODE:-gpu}
np=$((SLURM_JOB_NUM_NODES * 4))
[ "$mode" = cpu ] && np=$((SLURM_JOB_NUM_NODES * 96))
nsteps=${FM_NSTEPS:-20}
out=$R/export-ab-g$n-$mode-$SLURM_JOB_ID
mkdir -p "$out"
log "start, $n grains, $SLURM_JOB_NUM_NODES nodes, $mode, $nsteps steps," \
  "A $A, B $B, scripts $P"
export mesh=$W/meshes/part/v$n-np$np/mesh.
export vectors=$MESHES/v$n/vectors_${n}crystals.txt refinement=0
export nsteps duration=$nsteps solver_label=best modes=$mode
export run_timeout_gpu=${FM_TIMEOUT_GPU:-3000} paraview_export=yes
export run_timeout_cpu=${FM_TIMEOUT_CPU:-5400}
export solver_gpu="--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01"
export solver_cpu=$solver_gpu
options="--single-material --line-search Backtracking --prediction LinearExtrapolation --mesh_mode Restart --enable-export-von-Mises-stress"
export MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels
export MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered
failed=0
# usage: section name log: time of a section (maximum over the processes)
section() {
  grep -E -- "--> $1 " "$2" | head -1 | awk -F'|' '{gsub(/ /, "", $7); print $7}'
}
# usage: run label executable [options]
run() {
  local label=$1 e=$2 name d line
  shift 2
  name=$label-r0-$mode-n$SLURM_JOB_NUM_NODES-best
  executable=$e label=$label extra_options="$options $*" \
    "$P/cases/measure-polycrystal-grains.sh" "$out" > "$out/$label.out" 2>&1
  line=$(tail -1 "$out/series.txt")
  log "$label: $line"
  [ "$(awk '{print $5, $9}' <<< "$line")" = "0 $nsteps" ] || failed=1
  for d in "$out/$name"/*/; do
    [ -d "$d" ] || continue
    log "$label: $(basename "$d"): $(du -sb "$d" | cut -f1) bytes," \
      "$(find "$d" -type f | wc -l) files"
  done
  log "$label: total $(section NLEP::solve "$out/$name.txt") s of resolution," \
    "post-processings $(section NLEP::executePostProcessings "$out/$name.txt") s," \
    "ParaviewExportResults::Save $(section ParaviewExportResults::Save "$out/$name.txt") s," \
    "PolycrystalParaviewOutput::Save $(section PolycrystalParaviewOutput::Save "$out/$name.txt") s," \
    "L2Projection::Solve $(section L2Projection::Solve "$out/$name.txt") s"
}
run a "$A"
run b "$B"
run b32 "$B" --enable-paraview-single-precision
d=$out/%s-r0-$mode-n$SLURM_JOB_NUM_NODES-best
# relative tolerance of the comparisons of the fields: none on the cores
tolerance=${FM_TOLERANCE:-1e-6}
[ "$mode" = cpu ] && tolerance=""
# printf patterns of the directories of the computations
da=$(printf "$d" a) db=$(printf "$d" b) db32=$(printf "$d" b32)
for x in b b32; do
  y=$(printf "$d" $x)
  if cmp -s "$da.res" "$y.res" &&
    cmp -s "$da/avgStressPolycrystal" "$y/avgStressPolycrystal"; then
    log "$x: same .res and averages per grain as a"
  elif [ -n "$tolerance" ]; then
    log "$x: .res or averages per grain differ from a (not reproducible" \
      "on the GPUs), axial stress $(tail -1 "$y.res" | awk '{print $7}')" \
      "against $(tail -1 "$da.res" | awk '{print $7}')"
  else
    log "$x: .res or averages per grain differ from a"
    failed=1
  fi
done
compare() {
  local label=$1
  shift
  if python3 -I "$SLURM_SUBMIT_DIR/fm-chain/compare-vtu.py" "$@" \
    > "$out/compare-$label.txt" 2>&1; then
    log "$label: $(tail -1 "$out/compare-$label.txt")"
  else
    log "$label: DIFFERENCES ($(tail -1 "$out/compare-$label.txt"))"
    failed=1
  fi
}
compare split-a-b split "$da" "$db" $tolerance
if [ -z "$tolerance" ]; then
  compare float-b-b32 float "$db/PolycrystalParaviewResults" \
    "$db32/PolycrystalParaviewResults"
fi
if [ -n "${FM_REFERENCE_EXPORT:-}" ]; then
  if python3 -I "$SLURM_SUBMIT_DIR/fm-chain/compare-vtu.py" split \
    "$FM_REFERENCE_EXPORT" "$db" "${FM_TOLERANCE:-1e-6}" \
    > "$out/compare-reference-b.txt" 2>&1; then
    log "reference-b: $(grep -E 'largest|compared' "$out/compare-reference-b.txt" | tr '\n' ';')"
  else
    log "reference-b: DIFFERENCES beyond ${FM_TOLERANCE:-1e-6}" \
      "($(grep -E 'largest|compared' "$out/compare-reference-b.txt" | tr '\n' ';'))"
  fi
fi
energy "$out/energy.txt"
# the exported files are removed (quota of files), except the last cycle of
# b and b32
if [ "${FM_KEEP_EXPORT:-no}" != yes ]; then
  last=$(printf "Cycle%06d" $((nsteps - 1)))
  for x in b b32; do
    c=$(printf "$d" $x)/PolycrystalParaviewResults
    mkdir -p "$out/kept-$x"
    mv "$c/$last" "$c"/*.pvd "$out/kept-$x/" 2> /dev/null
  done
  rm -rf "$da/OutputFile-Uniaxial-polycrystal" "$da/PolycrystalParaviewOutput" \
    "$db/PolycrystalParaviewResults" "$db32/PolycrystalParaviewResults"
fi
if [ $failed = 0 ]; then
  touch "$out/passed"
  log "end, same fields, roundings in single precision, same outputs"
else
  log "end, FAILURES"
  exit 1
fi
