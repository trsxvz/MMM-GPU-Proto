#!/bin/bash
# Tests of mfem-mgis computing the mean values of the thermodynamic forces
# (BehaviourIntegratorsSelectionTest, MaterialArraysTest,
# BatchIntegrationTest), with the modified mfem-mgis built by build-fix.sh in
# FM_FIX (default <scratch>/fix), on a node of the MI300A partition.
#SBATCH --job-name=fm-fix-tests
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
F=${FM_FIX:-$W/fix}
out=$R/fix-tests-$SLURM_JOB_ID
mkdir -p "$out"
log "start, build of $F"
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
set -u
# launchers of the parallel tests (mpirun) of the machine profile
export PATH=${PROTO_MPI_LAUNCHERS:+$PROTO_MPI_LAUNCHERS:}$PATH
tests="BehaviourIntegratorsSelectionTest MaterialArraysTest BatchIntegrationTest"
cd "$F/build/mfem-mgis" || exit 1
# the behaviours of the tests (BehaviourTest and the generic-parallel
# behaviours of each backend) are built by the target check, not by the
# executables of the tests
behaviours="BehaviourTest $(cmake --build . --target help |
  sed -n 's/^\.\.\. \(GenericParallelBehaviour-[a-z]*\)$/\1/p' | tr '\n' ' ')"
log "behaviours of the tests: $behaviours"
# shellcheck disable=SC2086
if ! cmake --build . --target $tests $behaviours -j 32 > "$out/build.log" 2>&1; then
  log "build of the tests failed (see $out/build.log)"
  exit 1
fi
ctest -R "^(${tests// /|})" --output-on-failure > "$out/ctest.log" 2>&1
log "tests: $(grep -E 'tests passed|tests failed' "$out/ctest.log")" \
  "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout)' "$out/ctest.log" | tr '\n' ' ')"
log end
