#!/bin/bash
# Tests of mfem-mgis using the non linear solvers (line search, adaptive
# tolerance, fixed point, predictions, tests of the Newton solver), on the
# GPU of a node of the MI300A partition, with the mfem-mgis of the separate
# build of FM_FIX (build-fix.sh, variant B, run here unless FM_BUILD=0).
# The settings are the ones of scripts/mfem-mgis-gpu-suite.sh, one test at
# a time:
# - kernels: device of the machine profile (hip), ordered assembly, all the
#   computations forced to the kernels of the batch path;
# - default: same device, the computations selected by each test;
# - device: the tests whose name selects a device, without these settings.
# FM_TESTS: executables of the tests (also the prefixes of their names),
# FM_MODES: modes run (default "kernels default device").
#SBATCH --job-name=fm-ls-gpu-tests
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
out=$R/ls-gpu-tests-$SLURM_JOB_ID
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
unset OMP_PROC_BIND OMP_PLACES
tests=${FM_TESTS:-LineSearchTest NewtonLineSearchTest AdaptiveLinearSolverToleranceTest FixedPointSolverTest LinearExtrapolationPredictionTest UniaxialTensileTest PredictionTest PredictionTest2 ImposedPressureTest ImposedPressureTest2 TwoLayerCube StationaryNonLinearHeatTransferTest BatchIntegrationTest MicromorphicKernelsTest HeatTransferKernelsTest SubSteppingMessageTest}
cd "$F/build/mfem-mgis" || exit 1
# only the executables defined by this version of mfem-mgis
targets=$(cmake --build . --target help | sed -n 's/^\.\.\. \(.*\)$/\1/p')
selected=""
for t in $tests; do
  if grep -qx "$t" <<< "$targets"; then
    selected="$selected $t"
  else
    log "no executable $t in this build"
  fi
done
# the behaviours of the tests (BehaviourTest and the generic-parallel
# behaviours of each backend) are not built by the executables of the tests
behaviours="BehaviourTest $(grep -E '^GenericParallelBehaviour-[a-z]*$' <<< "$targets" | tr '\n' ' ')"
# shellcheck disable=SC2086
if ! cmake --build . --target $selected $behaviours -j 32 > "$out/build-tests.log" 2>&1; then
  log "build of the tests failed (see $out/build-tests.log)"
  exit 1
fi
selected=${selected# }
regex="^(${selected// /|})"
device=${PROTO_GPU_DEVICE:-cuda}
for mode in ${FM_MODES:-kernels default device}; do
  settings=(MFEM_MGIS_DEVICE="$device" MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=1
    MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered)
  select=(-R "$regex"
    -E "cuda|hip|gpu|GPU|debug|^TwoLayerCube-1-MatrixFreeExecution$")
  case $mode in
    kernels)
      settings+=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
        MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
        MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
        MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
      ;;
    default) ;;
    device)
      settings=()
      select=(-R "${regex}.*(hip|gpu|GPU)")
      ;;
  esac
  t0=$(date +%s)
  env "${settings[@]}" ctest -j 1 "${select[@]}" --output-on-failure \
    > "$out/ctest-$mode.log" 2>&1
  log "GPU ($device), mode $mode, $(($(date +%s) - t0)) s:" \
    "$(grep -E 'tests passed|tests failed|No tests were found' "$out/ctest-$mode.log")" \
    "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout|Not Run)' "$out/ctest-$mode.log" | tr '\n' ' ')"
done
log end
