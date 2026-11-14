#!/bin/bash
# Tests added after the review of the patches of 13-14/11/2026: weights of
# the integration points kept by the problems
# (IntegrationPointsWeightsTest), failure of the line search of NewtonK0
# (FixedPointSolverTest), material axes modified in place and invalidation
# of the kept rotation matrices (BatchIntegrationTest), option
# --single-material of the polycrystal compared to one material per grain.
# Separate build (build-fix.sh in FM_FIX, variant B, run here unless
# FM_BUILD=0), then:
# - the unit tests of mfem-mgis FM_TESTS and other tests using the
#   materials, the rotations and the non linear solvers (FM_OTHER_TESTS);
# - the tests of the polycrystal of mm-opera-hpc by ctest, with the
#   behaviours built here (FM_POLYCRYSTAL_CTEST=0 to skip them);
# - the comparison of the tests of the polycrystal run by hand with the
#   behaviour library of the deployed prototype, on one process and on two
#   (cmake/run-and-compare.py: largest relative differences of the main
#   output file and of the averages per grain).
#SBATCH --job-name=fm-review-tests
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
out=$R/review-tests-$SLURM_JOB_ID
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
tests=${FM_TESTS:-IntegrationPointsWeightsTest BatchIntegrationTest FixedPointSolverTest MaterialArraysTest}
others=${FM_OTHER_TESTS-BehaviourIntegratorsSelectionTest NonLinearEvolutionProblemPostProcessingsTest LinearExtrapolationPredictionTest LineSearchTest NewtonLineSearchTest AdaptiveLinearSolverToleranceTest UniaxialTensileTest PredictionTest PredictionTest2 ImposedPressureTest ImposedPressureTest2 TwoLayerCube StationaryNonLinearHeatTransferTest GradientsKernelsTest HeatTransferKernelsTest MicromorphicKernelsTest PostProcessingsKernelsTest}
cd "$F/build/mfem-mgis" || exit 1
targets=$(cmake --build . --target help | sed -n 's/^\.\.\. \(.*\)$/\1/p')
# the behaviours of the tests (BehaviourTest and the generic-parallel
# behaviours of each backend) are not built by the executables of the tests
behaviours="BehaviourTest $(grep -E '^GenericParallelBehaviour-[a-z]*$' <<< "$targets" | tr '\n' ' ')"
# usage: run_tests name tests: builds the executables which exist in this
# version of mfem-mgis and runs their tests
run_tests() {
  local selected="" t
  for t in $2; do
    if grep -qx "$t" <<< "$targets"; then
      selected="$selected $t"
    else
      log "no executable $t in this build"
    fi
  done
  selected=${selected# }
  # shellcheck disable=SC2086
  if ! cmake --build . --target $selected $behaviours -j 32 > "$out/build-$1.log" 2>&1; then
    log "$1: build failed (see $out/build-$1.log)"
    return
  fi
  ctest -R "^(${selected// /|})" -E "cuda|hip|gpu|GPU" --output-on-failure \
    > "$out/ctest-$1.log" 2>&1
  log "$1: $(grep -E 'tests passed|tests failed' "$out/ctest-$1.log")" \
    "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout|Not Run)' "$out/ctest-$1.log" | tr '\n' ' ')"
}
run_tests tests "$tests"
[ -n "$others" ] && run_tests others "$others"
# tests of the polycrystal by ctest (polycrystal/CMakeLists.txt of
# mm-opera-hpc), with the behaviours of this build
src=$F/src/mm-opera-hpc
exe=$F/bin/uniaxial-polycrystal-B
if [ "${FM_POLYCRYSTAL_CTEST:-1}" = 1 ]; then
  cd "$F/build/mm-opera-hpc-B" || exit 1
  if cmake --build . --target BehavioursPolycrystal -j 16 > "$out/build-behaviours.log" 2>&1; then
    (cd polycrystal && ctest -R '^UniaxialPolycrystal' -E "gpu" --output-on-failure) \
      > "$out/ctest-polycrystal.log" 2>&1
    log "polycrystal, ctest: $(grep -E 'tests passed|tests failed|No tests were found' "$out/ctest-polycrystal.log")" \
      "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout|Not Run)' "$out/ctest-polycrystal.log" | tr '\n' ' ')"
    log "polycrystal, ctest, comparisons: $(grep -E "maximum relative difference" "$out/ctest-polycrystal.log" | sed 's/^[0-9]*: //' | tr '\n' ';')"
  else
    log "polycrystal: build of the behaviours failed (see $out/build-behaviours.log)"
  fi
fi
# the same comparisons by hand, with the behaviour library of the deployed
# prototype
lib=$D/MMM-GPU-Proto/build/mm-opera-hpc/polycrystal/BehavioursPolycrystal/cpu/libBehaviour.so
for n in 1 2; do
  launcher=()
  [ $n = 2 ] && launcher=(mpirun -n 2)
  t0=$(date +%s)
  python3 -I "$src/cmake/run-and-compare.py" --directory "$out/comparison-n$n" \
    --compare polycrystal.res 1,2-4,5-7 "${FM_RES_TOLERANCE:-1e-9}" \
    --compare avgStressPolycrystal 1,2- 2e-5 --option=--single-material -- \
    "${launcher[@]}" "$exe" --disable-paraview-export --nstep 2 --duration 0.4 \
    --mesh "$src/polycrystal/mesh/5crystals.msh" \
    --vect "$src/polycrystal/mesh/vectors_5crystals.txt" \
    --macroscopic-stress-output-file polycrystal.res -l "$lib" \
    > "$out/comparison-n$n.txt" 2>&1
  log "polycrystal, single material against one per grain, $n process(es): status $?," \
    "$(($(date +%s) - t0)) s, $(grep -c 'Nombre de grains : *5' "$out/comparison-n$n.txt") line of the number of grains," \
    "$(grep -E "maximum relative difference|run-and-compare.py: " "$out/comparison-n$n.txt" | tr '\n' ';')"
done
log end
