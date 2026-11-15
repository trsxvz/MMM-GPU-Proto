#!/bin/bash
# Performance of the polycrystal on the GPUs (JOURNAL.md, 2026-11-15):
# separate build of MFEM with Umpire, MGIS, mfem-mgis and the polycrystal
# (build-perf.sh in FM_FIX, variants FM_BUILD_VARIANTS, default "umpire
# deployed") on the local disk of the node, then, in the same job (the build
# trees are lost at its end):
# - the tests FM_TESTS of mfem-mgis (executables, also the prefixes of the
#   names of their tests; default: the tests of the changed components and
#   the ones using the non linear solvers), for the variant umpire: on the
#   host, then on a MI300A with the settings of
#   scripts/mfem-mgis-gpu-suite.sh (modes kernels and default), then the
#   tests whose name selects a device;
# - the tests FM_DEPLOYED_TESTS on a MI300A for the variant deployed (MFEM
#   without Umpire), mode kernels;
# - the tests of the polycrystal by ctest (host), variant umpire;
# - the polycrystal of the tests (5 grains, 2 steps, final options) on the 4
#   GPUs with the variant umpire: each process seeing its GPU; all the
#   processes seeing the 4 GPUs (mfem-mgis then creates the allocators of
#   Umpire of MFEM on the GPU of each process, resources DEVICE::<id>: the
#   memory of each GPU is read with rocm-smi during the computation); and
#   with the mean stresses computed by kernels
#   (MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels).
# The job fails if the build or a test fails, so that a job submitted with
# --dependency=afterok does not start; when everything passes, the file
# <results>/perf-tests-<job>/passed is written (FM_REQUIRE of the next
# jobs, see common.sh).
#SBATCH --job-name=fm-perf-tests
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
variants=${FM_BUILD_VARIANTS:-umpire deployed}
out=$R/perf-tests-$SLURM_JOB_ID
mkdir -p "$out"
export FM_BUILD_DIR=${FM_BUILD_DIR:-${TMPDIR:-/tmp}/fix-perf-$SLURM_JOB_ID}
B=$FM_BUILD_DIR
failed=0
log "start, sources $F/src.tgz ($(sha256sum "$F/src.tgz" | cut -c1-16))," \
  "build in $B, variants $variants"
if ! JOBS=96 FM_BUILD_VARIANTS=$variants bash "$F/build-perf.sh" > "$out/build.log" 2>&1; then
  log "build failed (see $out/build.log)"
  rm -rf "$B"
  exit 1
fi
log "built: $(grep -E '^[a-zA-Z 0-9.]* built|: MFEM ' "$out/build.log" | tr '\n' ';')"
# env.sh redefines P, the deployed prototype: the scripts of the
# computations remain the ones of FM_PROTO
scripts=$P
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
set -u
P=$scripts
# launchers of the parallel tests (mpirun) of the machine profile
export PATH=${PROTO_MPI_LAUNCHERS:+$PROTO_MPI_LAUNCHERS:}$PATH
tests=${FM_TESTS:-GlobalStiffnessMatrixTest IntegrationPointsWeightsTest PostProcessingsKernelsTest MaterialArraysTest LineSearchTest NewtonLineSearchTest AdaptiveLinearSolverToleranceTest FixedPointSolverTest LinearExtrapolationPredictionTest UniaxialTensileTest PredictionTest PredictionTest2 ImposedPressureTest ImposedPressureTest2 TwoLayerCube StationaryNonLinearHeatTransferTest BatchIntegrationTest MicromorphicKernelsTest HeatTransferKernelsTest SubSteppingMessageTest}
deployed_tests=${FM_DEPLOYED_TESTS:-GlobalStiffnessMatrixTest IntegrationPointsWeightsTest FixedPointSolverTest NewtonLineSearchTest LinearExtrapolationPredictionTest}
device=${PROTO_GPU_DEVICE:-cuda}
# usage: run_ctest variant mode tests
# runs the tests of a variant in a mode: host, kernels, default (settings of
# scripts/mfem-mgis-gpu-suite.sh) or device (tests selecting a device)
run_ctest() {
  local v=$1 mode=$2 regex="^(${3// /|})" settings=() select=() t0 summary
  (
    unset OMP_PROC_BIND OMP_PLACES
    settings=(MFEM_MGIS_DEVICE="$device" MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=1
      MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered)
    select=(-j 1 -R "$regex"
      -E "cuda|hip|gpu|GPU|debug|^TwoLayerCube-1-MatrixFreeExecution$")
    case $mode in
      host)
        settings=(OMP_NUM_THREADS=1)
        select=(-j 8 -R "$regex" -E "cuda|hip|gpu|GPU")
        ;;
      kernels)
        settings+=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
          MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
          MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
          MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
        ;;
      default) ;;
      device)
        settings=(OMP_NUM_THREADS=1)
        select=(-j 1 -R "${regex}.*(hip|gpu|GPU)")
        ;;
    esac
    cd "$B/build/mfem-mgis-$v" || exit 1
    env "${settings[@]}" ctest "${select[@]}" --output-on-failure \
      > "$out/ctest-$v-$mode.log" 2>&1
  )
  summary=$(grep -E 'tests passed|tests failed|No tests were found' "$out/ctest-$v-$mode.log")
  log "$v, mode $mode: ${summary:-no summary}" \
    "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout|Not Run)' "$out/ctest-$v-$mode.log" | tr '\n' ' ')"
  grep -q '^100% tests passed' <<< "$summary" || failed=1
}
for v in $variants; do
  cd "$B/build/mfem-mgis-$v" || exit 1
  targets=$(cmake --build . --target help | sed -n 's/^\.\.\. \(.*\)$/\1/p')
  # the behaviours of the tests are not built by their executables
  behaviours="BehaviourTest $(grep -E '^GenericParallelBehaviour-[a-z]*$' <<< "$targets" | tr '\n' ' ')"
  selected=$tests
  [ "$v" = umpire ] || selected=$deployed_tests
  # shellcheck disable=SC2086
  if ! cmake --build . --target $selected $behaviours -j 48 > "$out/build-tests-$v.log" 2>&1; then
    log "$v: build of the tests failed (see $out/build-tests-$v.log)"
    failed=1
    continue
  fi
  if [ "$v" = umpire ]; then
    for mode in host kernels default device; do
      run_ctest "$v" "$mode" "$selected"
    done
  else
    run_ctest "$v" kernels "$selected"
  fi
done
# tests of the polycrystal by ctest (polycrystal/CMakeLists.txt of
# mm-opera-hpc), with the behaviours of this build
if [ -d "$B/build/mm-opera-hpc-umpire" ]; then
  cd "$B/build/mm-opera-hpc-umpire" || exit 1
  if cmake --build . --target BehavioursPolycrystal -j 16 > "$out/build-behaviours.log" 2>&1; then
    (cd polycrystal && ctest -R '^UniaxialPolycrystal' -E "gpu" --output-on-failure) \
      > "$out/ctest-polycrystal.log" 2>&1
    summary=$(grep -E 'tests passed|tests failed|No tests were found' "$out/ctest-polycrystal.log")
    log "polycrystal, ctest: ${summary:-no summary}" \
      "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout|Not Run)' "$out/ctest-polycrystal.log" | tr '\n' ' ')"
    grep -q '^100% tests passed' <<< "$summary" || failed=1
  else
    log "polycrystal: build of the behaviours failed (see $out/build-behaviours.log)"
    failed=1
  fi
fi
# polycrystal of the tests on the 4 GPUs, final options
# usage: g5 label
g5() {
  local src=$B/src/mm-opera-hpc/polycrystal/mesh line
  mesh=$src/5crystals.msh vectors=$src/vectors_5crystals.txt refinement=0 \
    nsteps=2 duration=0.4 modes=gpu label=$1 solver_label=perf \
    solver_gpu="--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01" \
    executable=$F/bin/uniaxial-polycrystal-umpire \
    extra_options="--single-material --line-search Backtracking --prediction LinearExtrapolation" \
    "$P/cases/measure-polycrystal-grains.sh" "$out" > "$out/$1.out" 2>&1
  line=$(tail -1 "$out/series.txt")
  log "polycrystal, 5 grains, $1: $line"
  # status 0 and the two steps
  [ "$(awk '{print $5, $9}' <<< "$line")" = "0 2" ] || failed=1
}
if [ -x "$F/bin/uniaxial-polycrystal-umpire" ]; then
  g5 g5-gpu-per-rank
  MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels g5 g5-kernels
  # all the GPUs seen by all the processes: memory used on each GPU, read
  # every second during the computation
  (while :; do
    rocm-smi --showmeminfo gtt --csv 2> /dev/null | tr '\n' ' '
    echo
    sleep 1
  done) > "$out/g5-all-gpus-memory.txt" &
  sampler=$!
  PROTO_GPU_PER_RANK='' g5 g5-all-gpus
  kill "$sampler" 2> /dev/null
  log "all the GPUs seen by all the processes, largest memory used per GPU (bytes):" \
    "$(python3 -I - "$out/g5-all-gpus-memory.txt" << 'PY'
import re, sys
largest = {}
for line in open(sys.argv[1]):
    for card, total, used in re.findall(r"(card\d+),(\d+),(\d+)", line):
        largest[card] = max(largest.get(card, 0), int(used))
print(" ".join(f"{c}={u}" for c, u in sorted(largest.items())))
PY
)"
fi
rm -rf "$B"
if [ $failed = 0 ]; then
  touch "$out/passed"
  log "end, all the tests passed"
else
  log "end, FAILURES (the next jobs of the chain do not start)"
  exit 1
fi
