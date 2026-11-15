#!/bin/bash
# Mutations of the code of the performance of the polycrystal on the GPUs
# (JOURNAL.md, 2026-11-15): each one must make its tests fail.
# - global matrix assembled from the local matrices
#   (GlobalStiffnessMatrixTest, on the host and on a MI300A): values kept by
#   hypre not freed, received values ignored, sent values of another entry,
#   prolongation left on the host;
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
scripts=$P
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
set -u
P=$scripts
export PATH=${PROTO_MPI_LAUNCHERS:+$PROTO_MPI_LAUNCHERS:}$PATH
unset OMP_PROC_BIND OMP_PLACES
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
cd "$B/build/mfem-mgis-umpire" || exit 1
targets=$(cmake --build . --target help | sed -n 's/^\.\.\. \(.*\)$/\1/p')
behaviours="BehaviourTest $(grep -E '^GenericParallelBehaviour-[a-z]*$' <<< "$targets" | tr '\n' ' ')"
device=${PROTO_GPU_DEVICE:-cuda}
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
# usage: run_tests name executables regex modes: builds the library and the
# executables, then runs the tests selected by the regular expression in
# each mode: host, kernels (kernels forced, on the host) or gpu (kernels
# forced, on a MI300A)
run_tests() {
  local mode settings
  cd "$B/build/mfem-mgis-umpire" || exit 1
  # shellcheck disable=SC2086
  if ! cmake --build . --target $2 $behaviours -j 48 > "$out/build-$1.log" 2>&1; then
    log "$1: build failed (see $out/build-$1.log)"
    return
  fi
  for mode in $4; do
    case $mode in
      host) settings=(OMP_NUM_THREADS=1) ;;
      kernels) settings=(OMP_NUM_THREADS=1 "${kernels[@]}") ;;
      gpu) settings=(MFEM_MGIS_DEVICE="$device"
        MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=1
        MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered "${kernels[@]}") ;;
    esac
    env "${settings[@]}" ctest -R "$3" -E "cuda|hip|gpu|GPU|debug" \
      --output-on-failure > "$out/ctest-$1-$mode.log" 2>&1
    log "$1, $mode: $(grep -E 'tests passed|tests failed|No tests were found' "$out/ctest-$1-$mode.log")," \
      "$(grep -c 'check failed' "$out/ctest-$1-$mode.log") failed checks;" \
      "first: $(grep -m 2 'check failed' "$out/ctest-$1-$mode.log" | tr '\n' ' ' | cut -c1-330)"
  done
}
# usage: run_polycrystal name: builds the polycrystal and runs the test of
# the options of BoomerAMG
run_polycrystal() {
  cd "$B/build/mm-opera-hpc-umpire" || exit 1
  if ! cmake --build . --target uniaxial-polycrystal BehavioursPolycrystal -j 16 > "$out/build-$1.log" 2>&1; then
    log "$1: build failed (see $out/build-$1.log)"
    return
  fi
  (cd polycrystal && ctest -R '^UniaxialPolycrystalTest-BoomerAMGParameters' --output-on-failure) \
    > "$out/ctest-$1.log" 2>&1
  log "$1: $(grep -E 'tests passed|tests failed|No tests were found' "$out/ctest-$1.log")" \
    "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout|Not Run)' "$out/ctest-$1.log" | tr '\n' ' ')"
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
global_tests() {
  run_tests "$1" GlobalStiffnessMatrixTest '^GlobalStiffnessMatrixTest-np[124]$' "host gpu"
}
means_tests() {
  run_tests "$1" IntegrationPointsWeightsTest '^IntegrationPointsWeightsTest' "host gpu"
}
solvers_tests() {
  run_tests "$1" "FixedPointSolverTest NewtonLineSearchTest LinearExtrapolationPredictionTest" \
    '^(FixedPointSolverTest|NewtonLineSearchTest|LinearExtrapolationPredictionTest)' "host kernels"
}
if selected original; then
  global_tests original-global
  means_tests original-means
  solvers_tests original-solvers
  run_polycrystal original-polycrystal
fi
# global matrix
selected g1 && mutate g1-values-kept-by-hypre-not-freed "$global" \
  '      if (use_device) {
        freeValuesKeptByHypre(diag);
        freeValuesKeptByHypre(offd);
      }' '' &&
  global_tests g1-values-kept-by-hypre-not-freed
selected g2 && mutate g2-received-values-ignored "$global_kernels" \
  '                                       : a.received_values[c - a.local_entries];' \
  '                                       : real{0};' &&
  global_tests g2-received-values-ignored
selected g3 && mutate g3-sent-values-of-another-entry "$global_kernels" \
  '      s[i] = a[entries[i]];' '      s[i] = a[entries[0]];' &&
  global_tests g3-sent-values-of-another-entry
selected g4 && mutate g4-prolongation-left-on-the-host "$global" \
  '    P.HypreRead();
    return T;' '    return T;' &&
  global_tests g4-prolongation-left-on-the-host
# mean thermodynamic forces computed by kernels
selected k1 && mutate k1-weights-ignored "$points_kernels" \
  '          s += weights[i] * v[i * stride];' \
  '          s += v[i * stride];' &&
  means_tests k1-weights-ignored
selected k2 && mutate k2-stride-of-the-points-ignored "$points_kernels" \
  '          s += weights[i] * v[i * stride];' \
  '          s += weights[i] * v[i];' &&
  means_tests k2-stride-of-the-points-ignored
selected k3 && mutate k3-kernels-never-selected "$base" \
  '      if (!internals::usePostProcessingsKernels(may_throw)) {
        return false;
      }
      const auto ow = this->getIntegrationPointsWeights(ctx, m, b);' \
  '      if (true || (!internals::usePostProcessingsKernels(may_throw))) {
        return false;
      }
      const auto ow = this->getIntegrationPointsWeights(ctx, m, b);' &&
  means_tests k3-kernels-never-selected
selected k4 && mutate k4-measure-doubled "$base" \
  '      measure = workspace.measure;' \
  '      measure = 2 * workspace.measure;' &&
  means_tests k4-measure-doubled
selected k5 && mutate k5-thermodynamic-forces-not-rotated "$base" \
  '      auto forces = mgis::FieldView<const real>{};
      if (om->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {' \
  '      auto forces = mgis::FieldView<const real>{};
      if (false && (om->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC)) {' &&
  means_tests k5-thermodynamic-forces-not-rotated
# work vectors
selected v1 && mutate v1-same-work-vector "$solver" \
  '    auto &w = this->work_vectors[static_cast<std::size_t>(v)];' \
  '    auto &w = this->work_vectors[static_cast<std::size_t>(v) % 1];' &&
  solvers_tests v1-same-work-vector
selected v2 && mutate v2-inner-forces-not-added "$problem" \
  '      P->MultTranspose(f, fa);
      r += fa;' \
  '      P->MultTranspose(f, fa);' &&
  solvers_tests v2-inner-forces-not-added
# options of BoomerAMG of the polycrystal
selected p1 && mutate p1-options-of-boomeramg-not-given "$polycrystal" \
  '    if (boomeramg.empty()) {
      return std::string{s};
    }' \
  '    if (true || boomeramg.empty()) {
      return std::string{s};
    }' &&
  run_polycrystal p1-options-of-boomeramg-not-given
rm -rf "$B"
log end
