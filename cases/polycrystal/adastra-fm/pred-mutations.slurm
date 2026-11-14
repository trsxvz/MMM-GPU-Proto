#!/bin/bash
# Mutations of the prediction by linear extrapolation of mfem-mgis and of
# its fallback on the default prediction (integration failure at the
# extrapolated estimate): each one must make
# LinearExtrapolationPredictionTest fail in at least one of its variants
# (sequential, 2 processes, integration by batches and kernels on the
# host). Separate build in FM_FIX, not used by the measurements (the library
# built last is a mutated one). FM_FALLBACK_ONLY=1: the mutations of the
# fallback only.
#SBATCH --job-name=fm-pred-mutations
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
out=$R/pred-mutations-$SLURM_JOB_ID
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
cd "$F/build/mfem-mgis" || exit 1
base=$F/src/mfem-mgis/src/NonLinearEvolutionProblemImplementationBase.cxx
impl=$F/src/mfem-mgis/src/NonLinearEvolutionProblemImplementation.cxx
cp "$base" "$out/original-base.cxx"
cp "$impl" "$out/original-impl.cxx"
t=LinearExtrapolationPredictionTest
mesh=$F/src/mfem-mgis/tests/cube.mesh
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
# usage: run name
run() {
  local s p k
  if ! cmake --build . --target $t BehaviourTest -j 32 > "$out/build-$1.log" 2>&1; then
    log "$1: build failed"
    return
  fi
  lib=$(find . -name 'libBehaviourTest.so' | head -1)
  ./tests/$t --mesh "$mesh" --library "$lib" > "$out/$1-sequential.txt" 2>&1
  s=$?
  mpirun -n 2 ./tests/$t --parallel --mesh "$mesh" --library "$lib" > "$out/$1-parallel.txt" 2>&1
  p=$?
  env "${kernels[@]}" ./tests/$t --device cpu --mesh "$mesh" --library "$lib" > "$out/$1-kernels.txt" 2>&1
  k=$?
  log "$1: sequential status $s, $(grep -c 'check failed' "$out/$1-sequential.txt") failed checks," \
    "parallel status $p, $(grep -c 'check failed' "$out/$1-parallel.txt") failed checks," \
    "kernels status $k, $(grep -c 'check failed' "$out/$1-kernels.txt") failed checks;" \
    "warnings $(grep -c 'starts again from the default prediction' "$out/$1-sequential.txt");" \
    "first: $(cat "$out/$1-sequential.txt" "$out/$1-parallel.txt" "$out/$1-kernels.txt" | grep -m 2 'check failed' | tr '\n' ' ' | cut -c1-330)"
}
run original
# usage: mutate name file pattern replacement
mutate() {
  cp "$out/original-base.cxx" "$base"
  cp "$out/original-impl.cxx" "$impl"
  if ! python3 -I - "$2" "$3" "$4" << 'PY'
import sys
p, a, b = sys.argv[1:4]
s = open(p).read()
if s.count(a) != 1:
    sys.exit(f"pattern found {s.count(a)} times")
open(p, "w").write(s.replace(a, b))
PY
  then
    log "$1: mutation not applied"
    return
  fi
  run "$1"
}
if [ "${FM_FALLBACK_ONLY:-0}" != 1 ]; then
  mutate m1-no-extrapolation "$base" \
    'if (first_resolution && this->previous_time_increment.has_value()) {' \
    'if (false && first_resolution && this->previous_time_increment.has_value()) {'
  mutate m2-each-resolution "$base" \
    'if (first_resolution && this->previous_time_increment.has_value()) {' \
    'if ((first_resolution || true) && this->previous_time_increment.has_value()) {'
  mutate m3-reference-of-the-prediction "$base" \
    'if ((*onorm > real{0}) &&' \
    'if (false && (*onorm > real{0}) &&'
  mutate m4-ratio-ignored "$base" \
    'add(this->u0, dt / (*(this->previous_time_increment)),' \
    'add(this->u0, real{1} + 0 * dt / (*(this->previous_time_increment)),'
  mutate m5-dirichlet-not-imposed "$base" \
    '    for (const auto& bc : this->dirichlet_boundary_conditions) {
      bc->updateImposedValues(this->u1, t + dt);
    }
    this->extrapolated_estimate_state =
        ExtrapolatedEstimateState::NOT_INTEGRATED;
  }  // end of extrapolateUnknowns' \
    '    static_cast<void>(t);
    this->extrapolated_estimate_state =
        ExtrapolatedEstimateState::NOT_INTEGRATED;
  }  // end of extrapolateUnknowns'
  mutate m6-revert-keeps-the-flag "$base" \
    '    this->hasStiffnessOperatorsBeenComputed = false;
    this->first_resolution_of_time_step = true;' \
    '    this->hasStiffnessOperatorsBeenComputed = false;'
  mutate m7-increment-not-released "$base" \
    '      this->previous_increment.Destroy();
      this->previous_time_increment.reset();' \
    '      static_cast<void>(this);'
  mutate m8-update-keeps-the-flag "$base" \
    '    this->u0 = this->u1;
    this->first_resolution_of_time_step = true;' \
    '    this->u0 = this->u1;'
fi
# fallback on the default prediction
mutate f1-no-fallback "$base" \
  '        if ((this->extrapolated_estimate_state ==
             ExtrapolatedEstimateState::INTEGRATION_FAILURE) &&' \
  '        if (false && (this->extrapolated_estimate_state ==
             ExtrapolatedEstimateState::INTEGRATION_FAILURE) &&'
mutate f2-estimate-not-restored "$base" \
  '          this->restoreDefaultPrediction(t, dt);
          this->solver->Mult(this->u0, this->u1);' \
  '          this->solver->Mult(this->u0, this->u1);'
mutate f3-imposed-values-not-restored "$base" \
  '    this->u1 = this->u0;
    for (const auto& bc : this->dirichlet_boundary_conditions) {
      bc->updateImposedValues(this->u1, t + dt);
    }
  }  // end of restoreDefaultPrediction' \
  '    this->u1 = this->u0;
    static_cast<void>(t);
    static_cast<void>(dt);
  }  // end of restoreDefaultPrediction'
mutate f4-prediction-operators-counted "$base" \
  '    if ((it == IntegrationType::PREDICTION_TANGENT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_SECANT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_ELASTIC_OPERATOR)) {
      return;
    }
    this->extrapolated_estimate_state =' \
  '    static_cast<void>(it);
    this->extrapolated_estimate_state ='
mutate f5-fallback-on-any-failure "$base" \
  '        if ((this->extrapolated_estimate_state ==
             ExtrapolatedEstimateState::INTEGRATION_FAILURE) &&' \
  '        if ((this->extrapolated_estimate_state !=
             ExtrapolatedEstimateState::NONE) &&'
mutate f6-no-warning "$base" \
  '          if (isMPIMainProcess(this->fe_discretization)) {
            warning(ctx.log(),' \
  '          if (false && isMPIMainProcess(this->fe_discretization)) {
            warning(ctx.log(),'
mutate f7-failure-not-reported-point-by-point "$impl" \
  '        if (!this->mgis_integrator->integrate(e, tr, ue, it)) {
          this->mgis_integrator->setTimeIncrement(dt);
          this->reportBehavioursIntegration(it, false);' \
  '        if (!this->mgis_integrator->integrate(e, tr, ue, it)) {
          this->mgis_integrator->setTimeIncrement(dt);'
mutate f8-failure-not-reported-parallel "$impl" \
  '    this->reportBehavioursIntegration(it, success);' \
  '    this->reportBehavioursIntegration(it, true);'
mutate f9-failure-not-reported-batches "$impl" \
  '        this->mgis_integrator->setTimeIncrement(dt);
        this->reportBehavioursIntegration(it, false);
        return false;
      }
    } else {' \
  '        this->mgis_integrator->setTimeIncrement(dt);
        return false;
      }
    } else {'
mutate f10-estimate-kept-after-success "$base" \
  '    this->extrapolated_estimate_state =
        success ? ExtrapolatedEstimateState::INTEGRATED
                : ExtrapolatedEstimateState::INTEGRATION_FAILURE;' \
  '    if (!success) {
      this->extrapolated_estimate_state =
          ExtrapolatedEstimateState::INTEGRATION_FAILURE;
    }'
cp "$out/original-base.cxx" "$base"
cp "$out/original-impl.cxx" "$impl"
log end
