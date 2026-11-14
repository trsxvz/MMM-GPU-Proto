#!/bin/bash
# Prediction by linear extrapolation of mfem-mgis
# (LINEAR_EXTRAPOLATION_PREDICTION) and option --prediction of the
# polycrystal: separate build (build-fix.sh in FM_FIX, variant B), unit
# tests of mfem-mgis (FM_TESTS: the prediction and the line search,
# FM_OTHER_TESTS: other tests using the non linear solvers), then the
# polycrystal of the tests of mm-opera-hpc (5 grains, 4 time steps of
# 0.05 s, behaviour library of the deployed prototype), with and without
# the extrapolation, for NewtonK0 with the first choice of Eisenstat and
# Walker and the line search, and for the Newton solver: axial stresses at
# each time, lines of Newton iterations, ratios ||r||/||r_0|| at the
# iteration 0 of the resolutions (below 1: the estimate is better than the
# default one).
#SBATCH --job-name=fm-pred-tests
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
F=${FM_FIX:?}
out=$R/pred-tests-$SLURM_JOB_ID
mkdir -p "$out"
log "start, build in $F"
if [ "${FM_BUILD:-1}" = 1 ]; then
  if ! FM_BUILD_VARIANTS=B bash "$F/build-fix.sh" > "$F/build.log" 2>&1; then
    log "build failed (see $F/build.log)"
    cp "$F/build.log" "$out/build.log"
    exit 1
  fi
  cp "$F/build.log" "$out/build.log"
  log "built: $(grep -E '^(A|B) built' "$F/build.log" | tr '\n' ' ')"
fi
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
set -u
# launchers of the parallel tests (mpirun) of the machine profile
export PATH=${PROTO_MPI_LAUNCHERS:+$PROTO_MPI_LAUNCHERS:}$PATH
tests=${FM_TESTS:-LinearExtrapolationPredictionTest LineSearchTest NewtonLineSearchTest AdaptiveLinearSolverToleranceTest FixedPointSolverTest}
others=${FM_OTHER_TESTS-UniaxialTensileTest PredictionTest PredictionTest2 ImposedPressureTest ImposedPressureTest2 TwoLayerCube StationaryNonLinearHeatTransferTest}
cd "$F/build/mfem-mgis" || exit 1
# shellcheck disable=SC2086
if cmake --build . --target $tests BehaviourTest -j 32 > "$out/build-tests.log" 2>&1; then
  ctest -R "^(${tests// /|})" --output-on-failure > "$out/ctest.log" 2>&1
  log "tests: $(grep -E 'tests passed|tests failed' "$out/ctest.log")" \
    "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout)' "$out/ctest.log" | tr '\n' ' ')"
else
  log "build of the tests failed (see $out/build-tests.log)"
fi
# shellcheck disable=SC2086
if [ -n "$others" ] &&
  cmake --build . --target $others BehaviourTest -j 32 > "$out/build-others.log" 2>&1; then
  ctest -R "^(${others// /|})" --output-on-failure > "$out/ctest-others.log" 2>&1
  log "other tests: $(grep -E 'tests passed|tests failed' "$out/ctest-others.log")" \
    "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout|Not Run)' "$out/ctest-others.log" | tr '\n' ' ')"
elif [ -n "$others" ]; then
  log "build of the other tests failed (see $out/build-others.log)"
fi
# polycrystal of the tests (polycrystal/CMakeLists.txt of mm-opera-hpc)
cd "$F/build/mm-opera-hpc-B/polycrystal" || exit 1
exe=$F/bin/uniaxial-polycrystal-B
lib=$D/MMM-GPU-Proto/build/mm-opera-hpc/polycrystal/BehavioursPolycrystal/cpu/libBehaviour.so
check=$F/src/mm-opera-hpc/cmake/run-and-check.py
k0="--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --line-search Backtracking"
"$exe" --disable-post-processing --nstep 1 --duration 0.2 -l "$lib" \
  --prediction Unknown > "$out/polycrystal-prediction-unknown.txt" 2>&1
log "polycrystal, unknown prediction: status $?," \
  "$(grep -c 'error: invalid prediction' "$out/polycrystal-prediction-unknown.txt") error message"
# usage: ratios log: ratios ||r||/||r_0|| at the iteration 0
ratios() {
  sed -n 's/.*iteration  0 *: ||r|| = [^,]*, ||r||\/||r_0|| = //p' "$1" | tr '\n' ' '
}
for t in "k0-default:$k0" "k0-extrapolation:$k0 --prediction LinearExtrapolation" \
  "newton-default:" "newton-extrapolation:--prediction LinearExtrapolation"; do
  name=${t%%:*}
  opts=${t#*:}
  t0=$(date +%s)
  # shellcheck disable=SC2086
  python3 -I "$check" --file "$out/polycrystal-$name.res" --tolerance 0.01 \
    --check 0.2 7 1.7994e7 -- "$exe" --disable-post-processing --nstep 4 \
    --duration 0.2 -l "$lib" $opts \
    --macroscopic-stress-output-file "$out/polycrystal-$name.res" \
    > "$out/polycrystal-$name.txt" 2>&1
  log "polycrystal, $name: status $?, $(($(date +%s) - t0)) s," \
    "Szz $(awk 'NF >= 7 {printf "%s ", $7}' "$out/polycrystal-$name.res" 2> /dev/null)," \
    "$(grep -c '||r||/||r_0||' "$out/polycrystal-$name.txt") lines of iterations," \
    "$(grep -c 'the resolution failed' "$out/polycrystal-$name.txt") failed resolutions," \
    "ratios at the iteration 0: $(ratios "$out/polycrystal-$name.txt")"
done
for s in k0 newton; do
  log "polycrystal, $s: relative differences of Szz with and without extrapolation:" \
    "$(paste "$out/polycrystal-$s-default.res" "$out/polycrystal-$s-extrapolation.res" |
      awk 'NF >= 14 {d = ($14 - $7) / $7; printf "%.3g ", d}')"
done
log end
