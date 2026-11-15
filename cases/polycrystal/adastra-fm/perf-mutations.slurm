#!/bin/bash
# Mutations of the code of the performance of the polycrystal on the GPUs
# (JOURNAL.md, 2026-11-15): each one must make its tests fail.
# - global matrix assembled from the local matrices
#   (GlobalStiffnessMatrixTest, on the host and on a MI300A): values kept by
#   hypre not freed, transposed blocks kept by MFEM not freed, received
#   values ignored, sent values of another entry;
# - mean thermodynamic forces computed by kernels
#   (IntegrationPointsWeightsTest, on the host and on a MI300A): weights
#   ignored, stride of the points ignored, computation by kernels never
#   selected, measure doubled, thermodynamic forces of an orthotropic
#   behaviour not rotated;
# - work vectors of the non linear solvers and of the residual
#   (FixedPointSolverTest, NewtonLineSearchTest,
#   LinearExtrapolationPredictionTest, with the kernels forced): the same
#   vector for all the uses, inner forces not added to the residual on the
#   true degrees of freedom;
# - option --hypre-boomeramg-parameters of the polycrystal
#   (UniaxialPolycrystalTest-BoomerAMGParameters): options not given.
# Separate build on the local disk of the node (build-perf.sh in FM_FIX,
# variant umpire: the MFEM with Umpire of FM_FIX must exist), in another
# directory than the one of the executables of the measurements, which are
# not replaced (FM_FIX of build-perf.sh is a directory of the local disk
# holding a link to the sources and to the installation of MFEM).
# FM_MUTATIONS selects the mutations (default: all).
# A mutation is detected if its tests fail in at least one mode. The job
# fails, so that a job submitted with --dependency=afterok does not start,
# if the tests of the original sources fail, if a mutation is not applied
# or not built, or if a mutation is not detected; otherwise the file
# <results>/perf-mutations-<job>/passed is written.
#SBATCH --job-name=fm-perf-mutations
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
out=$R/perf-mutations-$SLURM_JOB_ID
mkdir -p "$out"
B=${TMPDIR:-/tmp}/fix-perf-mutations-$SLURM_JOB_ID
# the executables built here are mutated: they are kept on the local disk
L=$B/fix
mkdir -p "$L"
ln -s "$F/src.tgz" "$F/install" "$L/"
log "start, sources $F/src.tgz ($(sha256sum "$F/src.tgz" | cut -c1-16)), build in $B"
build() {
  JOBS=96 FM_FIX=$L FM_BUILD_DIR=$B FM_BUILD_VARIANTS=umpire \
    bash "$F/build-perf.sh" > "$out/build.log" 2>&1
}
if ! build; then
  log "build failed (see $out/build.log)"
  rm -rf "$B"
  exit 1
fi
load_build_environment
src=$B/src/mfem-mgis
opera=$B/src/mm-opera-hpc
global=$src/src/Batch/GlobalStiffnessMatrix.cxx
global_kernels=$src/src/Batch/GlobalStiffnessMatrixKernels.cxx
points_kernels=$src/src/Batch/IntegrationPointsFunctionsKernels.cxx
base=$src/src/NonLinearEvolutionProblemImplementationBase.cxx
problem=$src/src/NonLinearEvolutionProblemImplementation.cxx
solver=$src/src/NonLinearSolverBase.cxx
polycrystal=$opera/polycrystal/UniaxialPolycrystal.cxx
files=("$global" "$global_kernels" "$points_kernels" "$base" "$problem"
  "$solver" "$polycrystal")
mkdir -p "$out/original"
for f in "${files[@]}"; do
  cp "$f" "$out/original/$(basename "$f")"
done
restore() {
  local f
  for f in "${files[@]}"; do
    cmp -s "$out/original/$(basename "$f")" "$f" ||
      cp "$out/original/$(basename "$f")" "$f"
  done
}
# shellcheck disable=SC2046
behaviours=$(behaviour_targets "$B/build/mfem-mgis-umpire" | tr '\n' ' ')
# outcome of the tests of the current sources: failed (in at least one
# mode), passed, or broken (not built); names of the mutations which are
# not detected, of the original tests which fail
outcome=passed
undetected=""
# usage: run_tests name executables regex modes: builds the library and the
# executables, then runs the tests selected by the regular expression in
# each mode of test_settings (common.sh)
run_tests() {
  local mode settings summary
  cd "$B/build/mfem-mgis-umpire" || exit 1
  # shellcheck disable=SC2086
  if ! cmake --build . --target $2 $behaviours -j 48 > "$out/build-$1.log" 2>&1; then
    log "$1: build failed (see $out/build-$1.log)"
    outcome=broken
    return
  fi
  for mode in $4; do
    mapfile -t settings < <(test_settings "$mode")
    env "${settings[@]}" ctest -R "$3" -E "cuda|hip|gpu|GPU|debug" \
      --output-on-failure > "$out/ctest-$1-$mode.log" 2>&1
    summary=$(grep -E 'tests passed|tests failed|No tests were found' "$out/ctest-$1-$mode.log")
    log "$1, $mode: ${summary:-no summary}," \
      "$(grep -c 'check failed' "$out/ctest-$1-$mode.log") failed checks;" \
      "first: $(grep -m 2 'check failed' "$out/ctest-$1-$mode.log" | tr '\n' ' ' | cut -c1-330)"
    if ! grep -q '^100% tests passed' <<< "$summary"; then
      [ "$outcome" = broken ] || outcome=failed
    fi
  done
}
# usage: run_polycrystal name: builds the polycrystal and runs the tests of
# the options of the linear solver of the elastic operator
run_polycrystal() {
  local summary
  cd "$B/build/mm-opera-hpc-umpire" || exit 1
  if ! cmake --build . --target uniaxial-polycrystal BehavioursPolycrystal -j 16 > "$out/build-$1.log" 2>&1; then
    log "$1: build failed (see $out/build-$1.log)"
    outcome=broken
    return
  fi
  (cd polycrystal && ctest -R '^UniaxialPolycrystalTest-(BoomerAMG|AmgX)Parameters' --output-on-failure) \
    > "$out/ctest-$1.log" 2>&1
  summary=$(grep -E 'tests passed|tests failed|No tests were found' "$out/ctest-$1.log")
  log "$1: ${summary:-no summary}" \
    "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout|Not Run)' "$out/ctest-$1.log" | tr '\n' ' ')"
  grep -q '^100% tests passed' <<< "$summary" || outcome=failed
}
# usage: verdict name expected: compares the outcome of the tests which
# just ran to the expected one (passed for the original sources, failed for
# a mutation)
verdict() {
  if [ "$outcome" != "$2" ]; then
    log "$1: tests $outcome, $2 expected"
    undetected="$undetected $1"
  fi
  outcome=passed
}
# usage: mutate name file pattern replacement [file pattern replacement...]:
# each pattern must be found once in its file
mutate() {
  local name=$1
  shift
  restore
  if ! python3 -I - "$@" << 'PY'
import sys
args = sys.argv[1:]
for i in range(0, len(args), 3):
    p, a, b = args[i:i + 3]
    s = open(p).read()
    if s.count(a) != 1:
        sys.exit(f"pattern found {s.count(a)} times in {p}")
    open(p, "w").write(s.replace(a, b))
PY
  then
    log "$name: mutation not applied"
    undetected="$undetected $name"
    return 1
  fi
}
selected() {
  case " ${FM_MUTATIONS:-all} " in
    *" all "* | *" $1 "*) return 0 ;;
  esac
  return 1
}
global_tests() {
  run_tests "$1" GlobalStiffnessMatrixTest '^GlobalStiffnessMatrixTest-np[124]$' "host gpu-kernels"
}
means_tests() {
  run_tests "$1" IntegrationPointsWeightsTest '^IntegrationPointsWeightsTest' "host gpu-kernels"
}
solvers_tests() {
  run_tests "$1" "FixedPointSolverTest NewtonLineSearchTest LinearExtrapolationPredictionTest" \
    '^(FixedPointSolverTest|NewtonLineSearchTest|LinearExtrapolationPredictionTest)' "host host-kernels"
}
# usage: mutation name tests file pattern replacement [file pattern
# replacement...]: applies a mutation, runs its tests, which must fail
mutation() {
  local name=$1 tests=$2
  shift 2
  selected "${name%%-*}" || return 0
  mutate "$name" "$@" || return 0
  "$tests" "$name"
  verdict "$name" failed
}
if selected original; then
  global_tests original-global
  verdict original-global passed
  means_tests original-means
  verdict original-means passed
  solvers_tests original-solvers
  verdict original-solvers passed
  run_polycrystal original-polycrystal
  verdict original-polycrystal passed
fi
# global matrix
mutation g1-values-kept-by-hypre-not-freed global_tests "$global" \
  '      if (use_device) {
        freeValuesKeptByHypre(diag);
        freeValuesKeptByHypre(offd);
      }' ''
mutation g2-received-values-ignored global_tests "$global_kernels" \
  '                                       : a.received_values[c - a.local_entries];' \
  '                                       : real{0};'
mutation g3-sent-values-of-another-entry global_tests "$global_kernels" \
  '      s[i] = a[entries[i]];' '      s[i] = a[entries[0]];'
mutation g4-transposed-blocks-not-freed global_tests "$global" \
  '      this->matrix->ResetTranspose();' ''
# mean thermodynamic forces computed by kernels
mutation k1-weights-ignored means_tests "$points_kernels" \
  '          s += weights[i] * v[i * stride];' \
  '          s += v[i * stride];'
mutation k2-stride-of-the-points-ignored means_tests "$points_kernels" \
  '          s += weights[i] * v[i * stride];' \
  '          s += weights[i] * v[i];'
mutation k3-kernels-never-selected means_tests "$base" \
  '      if (!internals::usePostProcessingsKernels(may_throw)) {
        return false;
      }
      const auto ow = this->getIntegrationPointsWeights(ctx, m, b);' \
  '      if (true || (!internals::usePostProcessingsKernels(may_throw))) {
        return false;
      }
      const auto ow = this->getIntegrationPointsWeights(ctx, m, b);'
mutation k4-measure-doubled means_tests "$base" \
  '      measure = workspace.measure;' \
  '      measure = 2 * workspace.measure;'
mutation k5-thermodynamic-forces-not-rotated means_tests "$base" \
  '      auto forces = mgis::FieldView<const real>{};
      if (om->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {' \
  '      auto forces = mgis::FieldView<const real>{};
      if (false && (om->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC)) {'
# work vectors
mutation v1-same-work-vector solvers_tests "$solver" \
  '    auto &w = this->work_vectors[static_cast<std::size_t>(v)];' \
  '    auto &w = this->work_vectors[static_cast<std::size_t>(v) % 1];'
mutation v2-inner-forces-not-added solvers_tests "$problem" \
  '      P->MultTranspose(f, fa);
      r += fa;' \
  '      P->MultTranspose(f, fa);'
# options of the linear solver of the elastic operator of the polycrystal
mutation p1-options-of-boomeramg-not-given run_polycrystal "$polycrystal" \
  '    if (boomeramg.empty()) {
      return std::string{s};
    }' \
  '    if (true || boomeramg.empty()) {
      return std::string{s};
    }'
mutation p2-options-of-amgx-accepted-with-boomeramg run_polycrystal "$polycrystal" \
  '  if ((s != "AmgX") && ((!amgx.empty()) || p.amgx_single_precision)) {' \
  '  if (false && (s != "AmgX")) {'
restore
rm -rf "$B"
if [ -z "$undetected" ]; then
  touch "$out/passed"
  log "end, all the mutations detected"
else
  log "end, FAILURES:$undetected"
  exit 1
fi
