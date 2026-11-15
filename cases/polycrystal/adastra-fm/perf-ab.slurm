#!/bin/bash
# Performance of the polycrystal on the GPUs (JOURNAL.md, 2026-11-15): the
# final computation of the polycrystal (FM_GRAINS grains, mesh partitioned
# beforehand, one material for all the grains, line search, prediction by
# linear extrapolation, FM_NSTEPS steps over FM_DURATION s, 20 and 20 by
# default) on all the GPUs of the allocation, in the variants FM_VARIANTS,
# run one after the other at each of the FM_GPU_REPS repetitions, on the
# same nodes:
# - ref: the executable of the measurements of the FM (FM_REFERENCE,
#   default <scratch>/bin/uniaxial-polycrystal-pred), before the changes;
# - deployed: the executable of FM_FIX built with the MFEM of the deployed
#   prototype (without Umpire): work vectors, global matrix assembled
#   without product of hypre;
# - umpire: the one built with MFEM with Umpire;
# - umpire-kernels: the same, with the mean stresses computed by kernels
#   (MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels): the final
#   configuration;
# - amg-<name>: the final configuration with the options of BoomerAMG of the
#   elastic operator FM_BOOMERAMG_<name> (as --hypre-boomeramg-parameters,
#   with `:` for `=` and `+` for `,`, which sbatch --export does not
#   accept in a value).
# With FM_CPU_REPS, the variants ref and umpire are also run on all the
# cores (the assembly of the global matrix is the same code on the host),
# to be done in a job of its own with a few steps (FM_GPU_REPS=0,
# FM_NSTEPS=2, FM_DURATION=2): 20 steps of 400 grains take about 70 min on
# 96 cores.
# A line per computation in <results>/series.txt (see
# cases/measure-polycrystal-grains.sh), the time of the main sections of
# each computation in <results>/sections.txt (a section timed at several
# places of the profiling tree has one entry per place).
# The job fails, so that a job submitted with --dependency=afterok does not
# start, if a computation fails or does not end its steps, or if the
# computations differ; otherwise the file <results>/passed is written
# (FM_REQUIRE of the next jobs, see common.sh):
# - on the GPUs, the axial stress of each computation is compared to the one
#   of the first computation with the relative tolerance FM_TOLERANCE
#   (default 1e-5: the element matrices are added in an order which changes
#   from one run to the next, and the axial stress of 400 grains varied by
#   1.4e-6 between the runs of one executable, job 5982471);
# - on the cores, the main output file of each computation (.res) is the
#   one of the first computation, bit for bit.
#SBATCH --job-name=fm-perf-ab
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
F=${FM_FIX:?}
n=${FM_GRAINS:?}
np=$((SLURM_JOB_NUM_NODES * 4))
out=$R/perf-ab-g$n-$SLURM_JOB_ID
mkdir -p "$out"
reference=${FM_REFERENCE:-$W/bin/uniaxial-polycrystal-pred}
variants=${FM_VARIANTS:-ref deployed umpire umpire-kernels}
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
export refinement=0 nsteps=${FM_NSTEPS:-20} duration=${FM_DURATION:-20}
export solver_label=best
export run_timeout_gpu=${FM_TIMEOUT_GPU:-2400}
export run_timeout_cpu=${FM_TIMEOUT_CPU:-5400}
solver=${FM_SOLVER:---non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01}
export solver_gpu=$solver solver_cpu=$solver
options=${FM_EXTRA_OPTIONS:---single-material --line-search Backtracking --prediction LinearExtrapolation --mesh_mode Restart}
log "start, $n grains, $SLURM_JOB_NUM_NODES nodes, $nsteps steps," \
  "GPU x${FM_GPU_REPS:-2}, CPU x${FM_CPU_REPS:-0}, variants $variants," \
  "reference $reference, executables of $F, options $options," \
  "solver $solver, scripts $P"
failed=0
# usage: run_variant variant mode repetition
run_variant() {
  local v=$1 mode=$2 e extra=$options kernels= parameters name line
  case $v in
    ref) e=$reference ;;
    deployed) e=$F/bin/uniaxial-polycrystal-deployed ;;
    umpire) e=$F/bin/uniaxial-polycrystal-umpire ;;
    umpire-kernels)
      e=$F/bin/uniaxial-polycrystal-umpire
      kernels=Kernels
      ;;
    amg-*)
      e=$F/bin/uniaxial-polycrystal-umpire
      kernels=Kernels
      name=FM_BOOMERAMG_${v#amg-}
      parameters=${!name:?}
      parameters=${parameters//:/=}
      extra="$extra --hypre-boomeramg-parameters ${parameters//+/,}"
      ;;
    *)
      log "unknown variant $v"
      failed=1
      return
      ;;
  esac
  if [ "$mode" = gpu ]; then
    export mesh=${FM_MESH:-$W/meshes/part/v$n-np$np/mesh.}
  else
    export mesh=${FM_MESH_CPU:-$W/meshes/part/v$n-np$((SLURM_JOB_NUM_NODES * 96))/mesh.}
  fi
  name=g${n}c-$v-rep$3
  ${kernels:+env MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=$kernels} \
    env executable="$e" extra_options="$extra" label="$name" modes="$mode" \
    "$P/cases/measure-polycrystal-grains.sh" "$out" > "$out/$name-$mode.out" 2>&1
  line=$(tail -1 "$out/series.txt")
  log "$mode $v rep $3: $line"
  # status 0 and all the steps
  [ "$(awk '{print $5, $9}' <<< "$line")" = "0 $nsteps" ] || failed=1
  # time of the main sections (maximum over the processes) and number of
  # levels of the multigrid of the elastic operator
  python3 -I - "$out/$name-r0-$mode-n$SLURM_JOB_NUM_NODES-best.txt" "$name $mode" \
    >> "$out/sections.txt" << 'PY'
import re, sys
sections = ("NLEP::solve", "FPS::integrate", "NS::computeResidual",
            "NS::getJacobian", "MFEM::Mult(r,c)", "PRED::residual",
            "computeMacroscopicCauchyStress", "common::post_processing_step")
values, levels = {}, "-"
for line in open(sys.argv[1], errors="replace"):
    f = [x.strip() for x in line.split("|")]
    for s in sections:
        if "--> " + s in f and len(f) > f.index("--> " + s) + 4:
            i = f.index("--> " + s)
            values.setdefault(s, []).append(
                f"{f[i + 1]} calls, min {f[i + 2]}, mean {f[i + 3]}, max {f[i + 4]}")
    m = re.match(r"\s*Num levels = (\d+)", line)
    if m and levels == "-":
        levels = m.group(1)
print(sys.argv[2] + ": levels " + levels + "; "
      + "; ".join(f"{s}: {' + '.join(values.get(s, ['-']))}" for s in sections))
PY
}
for i in $(seq "${FM_GPU_REPS:-2}"); do
  for v in $variants; do
    run_variant "$v" gpu "$i"
  done
done
for i in $(seq "${FM_CPU_REPS:-0}"); do
  for v in ${FM_CPU_VARIANTS:-ref umpire}; do
    run_variant "$v" cpu "$i"
  done
done
energy "$out/energy.txt"
# axial stresses of the GPU computations (column 8 of series.txt) compared
# to the first one
gpu=$(awk '$2 == "gpu"' "$out/series.txt" | wc -l)
if ! awk -v tol="${FM_TOLERANCE:-1e-5}" '$2 == "gpu" {
    if (ref == "") ref = $8
    d = ($8 - ref) / ref; if (d < 0) d = -d
    if (d > tol) { print $1, $8, "differs from", ref; bad = 1 } }
  END { exit bad }' "$out/series.txt" > "$out/stresses.txt" 2>&1; then
  log "axial stresses: $(tr '\n' ';' < "$out/stresses.txt")"
  failed=1
fi
# main output files of the CPU computations compared to the first one
cpu=0 first=""
for f in "$out"/*-cpu-n"$SLURM_JOB_NUM_NODES"-best.res; do
  [ -e "$f" ] || continue
  cpu=$((cpu + 1))
  first=${first:-$f}
  if ! cmp -s "$first" "$f"; then
    log "$(basename "$f") differs from $(basename "$first")"
    failed=1
  fi
done
if [ $failed = 0 ]; then
  touch "$out/passed"
  log "end, all the computations ended: $gpu on the GPUs, same axial" \
    "stresses within ${FM_TOLERANCE:-1e-5}; $cpu on the cores, same output" \
    "files bit for bit"
else
  log "end, FAILURES (the next jobs of the chain do not start)"
  exit 1
fi
