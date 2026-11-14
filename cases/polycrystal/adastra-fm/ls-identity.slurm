#!/bin/bash
# Options disabled: the executable built with the line search of mfem-mgis
# (FM_EXECUTABLE_B) gives the outputs of the previous one (FM_EXECUTABLE_A,
# default: the one of <scratch>/fix3, mfem-mgis 0e949cd441 and mm-opera-hpc
# 79b0577), bit for bit, on the CPU (deterministic), 100 grains, 96
# processes, one material for all the grains:
# - i1-ew1, i1-ew2: cleaned mesh, 3 steps of 1 s, NewtonK0 with the first
#   and the second choice of Eisenstat and Walker (solvers of the series);
# - i2-newton: same case, Newton solver (HyprePCG and BoomerAMG) with the
#   second choice;
# - i3-raw: raw mesh, 5 steps, NewtonK0 with the first choice (integration
#   failures: estimate moved towards the previous one, tolerance after a
#   reduction).
# Compared: <computation>.res, the averages per grain, and the lines of the
# iterations of the logs.
#SBATCH --job-name=fm-ls-identity
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
out=$R/ls-identity-$SLURM_JOB_ID
mkdir -p "$out"
A=${FM_EXECUTABLE_A:-$W/fix3/bin/uniaxial-polycrystal-B}
B=${FM_EXECUTABLE_B:?}
log "start, A $A, B $B, scripts $P"
export refinement=1 solver_label=best
export run_timeout_gpu=2400 run_timeout_cpu=2400
export vectors=$MESHES/v100/vectors_100crystals.txt
k0="--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance"
ew1="$k0 EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01"
ew2="$k0 EisenstatWalker2 --adaptive-linear-solver-initial-tolerance 0.1"
newton="--non-linear-solver Newton --adaptive-linear-solver-tolerance"
newton="$newton EisenstatWalker2 --adaptive-linear-solver-initial-tolerance 0.1"
lines() {
  grep -E 'Fixed point iteration|Newton iteration|Fixed Point iteration|Solving time step|Solution at time|resolution failed|dividing|integration fail' "$1.txt"
}
# usage: run_case name mesh steps solver
run_case() {
  local a b v e
  export mesh=$2 nsteps=$3 solver_cpu=$4 solver_gpu=$4
  for v in A B; do
    e=$A
    [ $v = B ] && e=$B
    executable=$e extra_options=--single-material label=$1-$v modes=cpu \
      "$P/cases/measure-polycrystal-grains.sh" "$out" > "$out/$1-$v.out" 2>&1
    log "$1 $v: $(tail -1 "$out/series.txt")"
  done
  a=$out/$1-A-r1-cpu-n1-best
  b=$out/$1-B-r1-cpu-n1-best
  if [ -s "$a.res" ] && cmp "$a.res" "$b.res" > "$out/diff-$1.txt" 2>&1 &&
    diff -rq "$a" "$b" >> "$out/diff-$1.txt" 2>&1 &&
    diff <(lines "$a") <(lines "$b") >> "$out/diff-$1.txt" 2>&1; then
    log "$1: identical outputs ($(wc -c < "$a.res") bytes of .res," \
      "$(lines "$a" | wc -l) lines of iterations," \
      "$(grep -c 'integration failure' "$a.txt") estimates moved after an integration failure)"
  else
    log "$1: outputs differ (see $out/diff-$1.txt)"
  fi
}
clean=$W/meshes/clean/v100/100crystals.msh
for c in ${FM_CASES:-i1-ew1 i1-ew2 i2-newton i3-raw}; do
  case $c in
    i1-ew1) run_case $c "$clean" 3 "$ew1" ;;
    i1-ew2) run_case $c "$clean" 3 "$ew2" ;;
    i2-newton) run_case $c "$clean" 3 "$newton" ;;
    i3-raw) run_case $c "$MESHES/v100/100crystals.msh" 5 "$ew1" ;;
  esac
done
log end
