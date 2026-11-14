#!/bin/bash
# Line search of the non linear solvers of mfem-mgis: separate build
# (build-fix.sh in FM_FIX, variant B), unit tests of mfem-mgis on the line
# search (FM_TESTS), other tests of mfem-mgis using the non linear solvers
# (FM_OTHER_TESTS), then tests of the polycrystal of mm-opera-hpc run by
# hand with the behaviour library of the deployed prototype (the libraries
# of the separate build are not built): commands and checks of
# polycrystal/CMakeLists.txt.
#SBATCH --job-name=fm-ls-tests
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
out=$R/ls-tests-$SLURM_JOB_ID
mkdir -p "$out"
log "start, build in $F"
if ! FM_BUILD_VARIANTS=B bash "$F/build-fix.sh" > "$F/build.log" 2>&1; then
  log "build failed (see $F/build.log)"
  cp "$F/build.log" "$out/build.log"
  exit 1
fi
cp "$F/build.log" "$out/build.log"
log "built: $(grep -E '^(A|B) built' "$F/build.log" | tr '\n' ' ')"
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
set -u
# launchers of the parallel tests (mpirun) of the machine profile
export PATH=${PROTO_MPI_LAUNCHERS:+$PROTO_MPI_LAUNCHERS:}$PATH
tests=${FM_TESTS:-LineSearchTest NewtonLineSearchTest AdaptiveLinearSolverToleranceTest FixedPointSolverTest}
# other tests of the non linear solvers and of the predictions, unchanged
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
# tests of the polycrystal (polycrystal/CMakeLists.txt of mm-opera-hpc)
cd "$F/build/mm-opera-hpc-B/polycrystal" || exit 1
exe=$F/bin/uniaxial-polycrystal-B
lib=$D/MMM-GPU-Proto/build/mm-opera-hpc/polycrystal/BehavioursPolycrystal/cpu/libBehaviour.so
check=$F/src/mm-opera-hpc/cmake/run-and-check.py
t0=$(date +%s)
"$exe" --disable-post-processing --nstep 1 --duration 0.2 -l "$lib" \
  --non-linear-solver NewtonK0 \
  --adaptive-linear-solver-tolerance EisenstatWalker1 \
  --non-linear-solver-maximum-iterations 1 \
  --macroscopic-stress-output-file "$out/polycrystal-iterations.res" \
  > "$out/polycrystal-iterations.txt" 2>&1
s=$?
log "polycrystal, 1 iteration at most: status $s, $(($(date +%s) - t0)) s," \
  "$(grep -c 'dividing time step by 2' "$out/polycrystal-iterations.txt") divisions of the time step," \
  "$(grep -c 'Maximum number of substeps reached' "$out/polycrystal-iterations.txt") failure message," \
  "$(grep -c 'Fixed point iteration  2' "$out/polycrystal-iterations.txt") second iterations," \
  "$(grep -c 'the resolution failed' "$out/polycrystal-iterations.txt") failed resolutions"
"$exe" --disable-post-processing --nstep 1 -l "$lib" \
  --non-linear-solver-maximum-iterations 0 \
  > "$out/polycrystal-iterations-0.txt" 2>&1
log "polycrystal, 0 iteration at most: status $? (refusal expected)"
"$exe" --disable-post-processing --nstep 1 --duration 0.2 -l "$lib" \
  --line-search Unknown > "$out/polycrystal-line-search-unknown.txt" 2>&1
log "polycrystal, unknown line search: status $?," \
  "$(grep -c 'error: invalid line search' "$out/polycrystal-line-search-unknown.txt") error message"
for t in "default:" "line-search:--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --line-search Backtracking" \
  "line-search-newton:--line-search Backtracking"; do
  name=${t%%:*}
  opts=${t#*:}
  t0=$(date +%s)
  # shellcheck disable=SC2086
  python3 -I "$check" --file "$out/polycrystal-$name.res" --tolerance 0.01 \
    --check 0.2 7 1.7994e7 -- "$exe" --disable-post-processing --nstep 1 \
    --duration 0.2 -l "$lib" $opts \
    --macroscopic-stress-output-file "$out/polycrystal-$name.res" \
    > "$out/polycrystal-$name.txt" 2>&1
  log "polycrystal, $name: status $?, $(($(date +%s) - t0)) s," \
    "Szz $(awk 'NF >= 7 {s = $7} END {print s}' "$out/polycrystal-$name.res" 2> /dev/null)," \
    "$(grep -c '||r||/||r_0||' "$out/polycrystal-$name.txt") Newton iterations," \
    "$(grep -c 'line search' "$out/polycrystal-$name.txt") line search messages"
done
log end
