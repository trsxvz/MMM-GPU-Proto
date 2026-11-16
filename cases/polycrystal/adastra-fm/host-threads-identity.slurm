#!/bin/bash
# Setup of the host shared between the threads of OpenMP (option
# enable-openmp of mfem-mgis, JOURNAL.md, 2026-11-16): the outputs of a
# computation do not depend on the number of threads. The computation runs
# on the cores (deterministic): 100 grains, cleaned mesh without refinement,
# one material for all the grains, options of the measurements, FM_NSTEPS
# steps of 1 s (default 1), with the executable FM_EXECUTABLE and 4
# processes per node having 24 cores each, as the processes of a GPU
# computation, for each number of threads per process of FM_THREADS
# (default "1 7 24"; OMP_NUM_THREADS, the processes keeping their cores).
# Compared to the outputs of the first number: <computation>.res, the
# averages per grain and the lines of the iterations of the logs.
# The job fails if a computation fails or if outputs differ; otherwise the
# file <results>/passed is written.
#SBATCH --job-name=fm-host-threads-identity
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
e=${FM_EXECUTABLE:?}
out=$R/host-threads-identity-$SLURM_JOB_ID
mkdir -p "$out"
export mesh=$W/meshes/clean/v100/100crystals.msh
export vectors=$MESHES/v100/vectors_100crystals.txt
export refinement=0 nsteps=${FM_NSTEPS:-1} duration=${FM_NSTEPS:-1}
export modes=cpu cpu_threads=24 solver_label=best
export run_timeout_cpu=${FM_TIMEOUT_CPU:-2400}
export solver_cpu=${FM_SOLVER:---non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01}
export extra_options=${FM_EXTRA_OPTIONS:---single-material --line-search Backtracking --prediction LinearExtrapolation}
log "start, 100 grains, $nsteps steps, 4 processes of 24 cores," \
  "threads ${FM_THREADS:-1 7 24}, executable $e, scripts $P"
lines() {
  grep -E 'iteration|Solving time step|Solution at time|resolution failed|dividing|integration fail' "$1.txt"
}
failed=0
first=""
for t in ${FM_THREADS:-1 7 24}; do
  executable=$e label=threads$t cpu_environment="OMP_NUM_THREADS=$t" \
    "$P/cases/measure-polycrystal-grains.sh" "$out" > "$out/threads$t.out" 2>&1
  line=$(tail -1 "$out/series.txt")
  log "threads$t: $line"
  [ "$(awk '{print $5, $9}' <<< "$line")" = "0 $nsteps" ] || failed=1
  c=$out/threads$t-r0-cpu-n1-best
  if [ -z "$first" ]; then
    first=$c
  elif [ -s "$first.res" ] && cmp "$first.res" "$c.res" > "$out/diff-threads$t.txt" 2>&1 &&
    diff -rq "$first" "$c" >> "$out/diff-threads$t.txt" 2>&1 &&
    diff <(lines "$first") <(lines "$c") >> "$out/diff-threads$t.txt" 2>&1; then
    log "threads$t: outputs identical to the ones of $(basename "$first")" \
      "($(wc -c < "$c.res") bytes of .res, $(lines "$c" | wc -l) lines of iterations)"
  else
    log "threads$t: outputs differ from the ones of $(basename "$first")" \
      "(see $out/diff-threads$t.txt)"
    failed=1
  fi
done
energy "$out/energy.txt"
if [ $failed = 0 ]; then
  touch "$out/passed"
  log "end, same outputs bit for bit for all the numbers of threads"
else
  log "end, FAILURES (the next jobs of the chain do not start)"
  exit 1
fi
