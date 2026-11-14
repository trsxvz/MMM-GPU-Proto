#!/bin/bash
# Mutations of the code covered by the tests added after the review of the
# patches of 13-14/11/2026: each one must make its test fail.
# - weights of the integration points kept by the problems
#   (IntegrationPointsWeightsTest): weights doubled, weights of another
#   element, weights computed again at each call, weights ignored by the
#   integrals;
# - rotation matrices kept by the materials (BatchIntegrationTest,
#   MaterialArraysTest): invalidation without effect, transposed matrices
#   not invalidated, matrices computed at each call;
# - line search of NewtonK0 (FixedPointSolverTest): resolution going on
#   after a failed search, integration failure of a trial ignored;
# - option --single-material of the polycrystal (comparison to one material
#   per grain, cmake/run-and-compare.py, with the behaviour library of the
#   deployed prototype): all the points in the first grain, averages of the
#   next grain.
# Separate build in FM_FIX (build-fix.sh, variant B), which must not be the
# one of an executable used by measurements: the library built last is a
# mutated one. FM_MUTATIONS selects the mutations (default: all).
#SBATCH --job-name=fm-review-mutations
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
out=$R/review-mutations-$SLURM_JOB_ID
mkdir -p "$out"
log "start, build in $F"
if ! FM_BUILD_VARIANTS=B bash "$F/build-fix.sh" > "$F/build.log" 2>&1; then
  log "build failed (see $F/build.log)"
  exit 1
fi
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
set -u
export PATH=${PROTO_MPI_LAUNCHERS:+$PROTO_MPI_LAUNCHERS:}$PATH
src=$F/src/mfem-mgis
opera=$F/src/mm-opera-hpc
base=$src/src/NonLinearEvolutionProblemImplementationBase.cxx
means=$src/include/MFEMMGIS/NonLinearEvolutionProblemImplementation.ixx
material=$src/src/Material.cxx
fixed_point=$src/src/FixedPointSolver.cxx
polycrystal=$opera/polycrystal/UniaxialPolycrystal.cxx
files=("$base" "$means" "$material" "$fixed_point" "$polycrystal")
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
cd "$F/build/mfem-mgis" || exit 1
targets=$(cmake --build . --target help | sed -n 's/^\.\.\. \(.*\)$/\1/p')
behaviours="BehaviourTest $(grep -E '^GenericParallelBehaviour-[a-z]*$' <<< "$targets" | tr '\n' ' ')"
# usage: run_tests name executables regex: builds the library and the
# executables, then runs the tests selected by the regular expression
run_tests() {
  cd "$F/build/mfem-mgis" || exit 1
  # shellcheck disable=SC2086
  if ! cmake --build . --target $2 $behaviours -j 32 > "$out/build-$1.log" 2>&1; then
    log "$1: build failed (see $out/build-$1.log)"
    return
  fi
  ctest -R "$3" -E "cuda|hip|gpu|GPU|debug" --output-on-failure \
    > "$out/ctest-$1.log" 2>&1
  log "$1: $(grep -E 'tests passed|tests failed|No tests were found' "$out/ctest-$1.log")," \
    "$(grep -c 'check failed' "$out/ctest-$1.log") failed checks;" \
    "first: $(grep -m 2 'check failed' "$out/ctest-$1.log" | tr '\n' ' ' | cut -c1-330)"
}
# usage: run_polycrystal name: builds the polycrystal and compares one
# material for all the grains to one material per grain
run_polycrystal() {
  cd "$F/build/mm-opera-hpc-B" || exit 1
  if ! cmake --build . --target uniaxial-polycrystal -j 8 > "$out/build-$1.log" 2>&1; then
    log "$1: build failed (see $out/build-$1.log)"
    return
  fi
  local lib=$D/MMM-GPU-Proto/build/mm-opera-hpc/polycrystal/BehavioursPolycrystal/cpu/libBehaviour.so
  python3 -I "$opera/cmake/run-and-compare.py" --directory "$out/$1" \
    --compare polycrystal.res 1,2-4,5-7 1e-9 \
    --compare avgStressPolycrystal 1,2- 2e-5 --option=--single-material -- \
    "$F/build/mm-opera-hpc-B/polycrystal/uniaxial-polycrystal" \
    --disable-paraview-export --nstep 2 --duration 0.4 \
    --mesh "$opera/polycrystal/mesh/5crystals.msh" \
    --vect "$opera/polycrystal/mesh/vectors_5crystals.txt" \
    --macroscopic-stress-output-file polycrystal.res -l "$lib" \
    > "$out/$1.txt" 2>&1
  log "$1: status $?, $(grep -c 'MISMATCH' "$out/$1.txt") mismatches;" \
    "$(grep -E "maximum relative difference|run-and-compare.py: " "$out/$1.txt" | tr '\n' ';' | cut -c1-600)"
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
    return 1
  fi
}
selected() {
  case " ${FM_MUTATIONS:-all} " in
    *" all "* | *" $1 "*) return 0 ;;
  esac
  return 1
}
weights_tests() {
  run_tests "$1" IntegrationPointsWeightsTest '^IntegrationPointsWeightsTest'
}
rotation_tests() {
  run_tests "$1" "BatchIntegrationTest MaterialArraysTest" \
    '^(BatchIntegrationTest|MaterialArraysTest)$'
}
line_search_tests() {
  run_tests "$1" FixedPointSolverTest '^FixedPointSolverTest(-kernels-cpu)?$'
}
if selected original; then
  weights_tests original-weights
  rotation_tests original-rotations
  line_search_tests original-line-search
  run_polycrystal original-polycrystal
fi
# weights of the integration points
selected w1 && mutate w1-weights-doubled "$base" \
  '              obi->getIntegrationPointWeight(tr, ip);' \
  '              2 * obi->getIntegrationPointWeight(tr, ip);' &&
  weights_tests w1-weights-doubled
selected w2 && mutate w2-weights-of-another-element "$base" \
  '        auto& tr = *(fes.GetElementTransformation(ve[i]));' \
  '        auto& tr = *(fes.GetElementTransformation(ve[elements.Size() - 1 - i]));' &&
  weights_tests w2-weights-of-another-element
selected w3 && mutate w3-weights-computed-at-each-call "$base" \
  '    if (const auto p = this->integration_points_weights.find(key);
        p != this->integration_points_weights.end()) {' \
  '    if (const auto p = this->integration_points_weights.find(key);
        false && (p != this->integration_points_weights.end())) {' \
  "$base" \
  '    return &(this->integration_points_weights.emplace(key, std::move(w))' \
  '    return &(this->integration_points_weights.insert_or_assign(key, std::move(w))' &&
  weights_tests w3-weights-computed-at-each-call
selected w4 && mutate w4-weights-ignored-by-the-integrals "$means" \
  '            s[k] += wo * thforces(o, k);' \
  '            s[k] += thforces(o, k);' &&
  weights_tests w4-weights-ignored-by-the-integrals
# rotation matrices kept by the materials
selected r1 && mutate r1-invalidation-without-effect "$material" \
  '  void Material::invalidateRotationMatrices() noexcept {
    this->kept_rotation_matrices_computed = {false, false};' \
  '  void Material::invalidateRotationMatrices() noexcept {' &&
  rotation_tests r1-invalidation-without-effect
selected r2 && mutate r2-transposed-matrices-not-invalidated "$material" \
  '  void Material::invalidateRotationMatrices() noexcept {
    this->kept_rotation_matrices_computed = {false, false};' \
  '  void Material::invalidateRotationMatrices() noexcept {
    this->kept_rotation_matrices_computed[0] = false;' &&
  rotation_tests r2-transposed-matrices-not-invalidated
selected r3 && mutate r3-matrices-computed-at-each-call "$material" \
  '    if (!this->kept_rotation_matrices_computed[i]) {' \
  '    if (true || (!this->kept_rotation_matrices_computed[i])) {' &&
  rotation_tests r3-matrices-computed-at-each-call
# line search of NewtonK0
selected l1 && mutate l1-resolution-going-on-after-a-failed-search "$fixed_point" \
  '          add_iteration_information(norm, true, &ls_information);
          break;
        }
        if (ls.reductions != 0) {' \
  '          add_iteration_information(norm, true, &ls_information);
        }
        if (ls.reductions != 0) {' &&
  line_search_tests l1-resolution-going-on-after-a-failed-search
selected l2 && mutate l2-integration-failure-of-a-trial-ignored "$fixed_point" \
  '              if (!this->integrate(u, it_type)) {
                return {};
              }
              return update_residual(u);' \
  '              static_cast<void>(this->integrate(u, it_type));
              return update_residual(u);' &&
  line_search_tests l2-integration-failure-of-a-trial-ignored
restore
cd "$F/build/mfem-mgis" || exit 1
cmake --build . -j 32 > "$out/build-restored.log" 2>&1
cmake --install . > /dev/null 2>&1
# one material for all the grains of the polycrystal (with the restored
# library)
selected s1 && mutate s1-all-the-points-in-the-first-grain "$polycrystal" \
  '        grains.of_points[o] = g;' \
  '        grains.of_points[o] = 0 * g;' &&
  run_polycrystal s1-all-the-points-in-the-first-grain
selected s2 && mutate s2-averages-of-the-next-grain "$polycrystal" \
  '            const auto g = static_cast<std::size_t>(grains.of_points[o]);' \
  '            const auto g = static_cast<std::size_t>((grains.of_points[o] + 1) % grains.number);' &&
  run_polycrystal s2-averages-of-the-next-grain
restore
log end
