#!/bin/bash
# Profiling sections declared as ranges of roctx to the profilers of ROCm
# (`mgis::setProfilingMarkers`, option enable-roctx of mfem-mgis): separate
# build of MGIS, mfem-mgis and the polycrystal (build-roctx.sh in FM_FIX,
# variants FM_BUILD_VARIANTS) on the local disk of the node, then, in the
# same job (the build trees are lost at its end):
# - the test of the profiling of MGIS (ProfilingTest);
# - for each variant, the test ProfilingMarkersTest of mfem-mgis on the host
#   and on a GPU (markers set with the option and a GPU, and only then)
#   and, for the variants with roctx, under rocprofv3 (the ranges of the
#   sections and of the steps must be in the trace of the markers);
# - the polycrystal of the tests of mm-opera-hpc (5 grains, 2 steps, 4
#   GPUs) with two processes out of four under rocprofv3
#   (fm-chain/trace-ranks.sh), then the analysis of fm-chain/trace-zones.py;
#   with FM_ROCPROF_SYS=1, the same computation with one process under
#   rocprof-sys-sample;
# - the tests of the polycrystal by ctest for the variant off
#   (FM_POLYCRYSTAL_CTEST=0 to skip them), and the other tests of mfem-mgis
#   FM_OTHER_TESTS (none by default).
#SBATCH --job-name=fm-roctx-tests
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
variants=${FM_BUILD_VARIANTS:-off on}
out=$R/roctx-tests-$SLURM_JOB_ID
mkdir -p "$out"
export FM_BUILD_DIR=${FM_BUILD_DIR:-${TMPDIR:-/tmp}/fix-roctx-$SLURM_JOB_ID}
B=$FM_BUILD_DIR
log "start, sources $F/src.tgz ($(sha256sum "$F/src.tgz" | cut -c1-16))," \
  "build in $B ($(df -h --output=avail "$(dirname "$B")" | tail -1 | tr -d ' ') available), variants $variants"
if ! FM_BUILD_VARIANTS=$variants bash "$F/build-roctx.sh" > "$F/build.log" 2>&1; then
  log "build failed (see $F/build.log)"
  cp "$F/build.log" "$out/build.log"
  rm -rf "$B"
  exit 1
fi
cp "$F/build.log" "$out/build.log"
log "built: $(grep -E '^[a-z]* built|mfem-mgis roctx|Found ROCTX' "$F/build.log" | tr '\n' ';')"
# env.sh redefines P, the deployed prototype: the scripts of the
# computations remain the ones of FM_PROTO
scripts=$P
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
set -u
P=$scripts
# launchers of the parallel tests (mpirun) of the machine profile
export PATH=${PROTO_MPI_LAUNCHERS:+$PROTO_MPI_LAUNCHERS:}$PATH
# the header of the roctx of roctracer has the functions used by mfem-mgis
echo '#include <roctracer/roctx.h>
int main() { roctxRangePushA("a"); return roctxRangePop(); }' > "$out/legacy.cxx"
if g++ -fsyntax-only -I"$ROCM_PATH/include" "$out/legacy.cxx" > "$out/legacy.log" 2>&1; then
  log "header roctracer/roctx.h: roctxRangePushA and roctxRangePop declared"
else
  log "header roctracer/roctx.h: compilation failed (see $out/legacy.log)"
fi
(cd "$B/build/mgis" && ctest -R '^ProfilingTest$' --output-on-failure) > "$out/ctest-mgis.log" 2>&1
log "mgis: $(grep -E '% tests passed|No tests were found' "$out/ctest-mgis.log")"
for v in $variants; do
  cd "$B/build/mfem-mgis-$v" || exit 1
  if ! cmake --build . --target ProfilingMarkersTest -j 8 > "$out/build-tests-$v.log" 2>&1; then
    log "$v: build of ProfilingMarkersTest failed (see $out/build-tests-$v.log)"
    continue
  fi
  ctest -R '^ProfilingMarkersTest' --output-on-failure -V > "$out/ctest-$v.log" 2>&1
  log "$v: $(grep -E '% tests passed' "$out/ctest-$v.log")," \
    "$(grep -E '^ *[0-9]+/[0-9]+ Test' "$out/ctest-$v.log" | sed 's/ \.\.*//; s/^ *//' | tr '\n' ';')" \
    "$(grep -E "in the trace of the markers" "$out/ctest-$v.log" | sed 's/^[0-9]*: -- //' | tr '\n' ';')"
  log "$v: libraries of roctx of the polycrystal:" \
    "$(ldd "$F/bin/uniaxial-polycrystal-$v" | grep -E 'roctx|rocprofiler|MFrontGenericInterface' | awk '{print $3}' | tr '\n' ' ')"
done
# polycrystal of the tests, 2 steps, 4 GPUs, 2 processes out of 4 traced
v=${FM_TRACE_VARIANT:-on}
if [ -x "$F/bin/uniaxial-polycrystal-$v" ]; then
  src=$B/src/mm-opera-hpc/polycrystal/mesh
  export FM_TRACE_DIR=$out/trace FM_TRACE_RANKS=${FM_TRACE_RANKS:-0,2}
  mesh=$src/5crystals.msh vectors=$src/vectors_5crystals.txt refinement=0 \
    nsteps=2 duration=0.4 modes=gpu label=g5 solver_label=trace \
    solver_gpu="--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01" \
    executable=$F/bin/uniaxial-polycrystal-$v \
    extra_options="--single-material --line-search Backtracking --prediction LinearExtrapolation" \
    profiler="$SLURM_SUBMIT_DIR/fm-chain/trace-ranks.sh" \
    "$P/cases/measure-polycrystal-grains.sh" "$out" > "$out/g5.out" 2>&1
  log "polycrystal, 5 grains, traced ranks $FM_TRACE_RANKS: $(tail -1 "$out/series.txt")"
  log "traces: $(cd "$FM_TRACE_DIR" 2> /dev/null && du -sh -- * | tr '\t\n' ' ;')"
  python3 -I "$SLURM_SUBMIT_DIR/fm-chain/trace-zones.py" "$FM_TRACE_DIR" \
    > "$out/trace-zones.txt" 2>&1
  log "zones: status $?, $(grep -E '^mean of' "$out/trace-zones.txt")"
  # the same computation under rocprof-sys-sample (sampling of the host,
  # MPI), one process traced
  if [ "${FM_ROCPROF_SYS:-0}" = 1 ]; then
    FM_TRACE_DIR=$out/trace-sys FM_TRACE_RANKS=0 FM_TRACE_TOOL=rocprof-sys \
      mesh=$src/5crystals.msh vectors=$src/vectors_5crystals.txt refinement=0 \
      nsteps=2 duration=0.4 modes=gpu label=g5 solver_label=sys \
      solver_gpu="--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01" \
      executable=$F/bin/uniaxial-polycrystal-$v \
      extra_options="--single-material --line-search Backtracking --prediction LinearExtrapolation" \
      profiler="$SLURM_SUBMIT_DIR/fm-chain/trace-ranks.sh" \
      "$P/cases/measure-polycrystal-grains.sh" "$out" > "$out/g5-sys.out" 2>&1
    log "polycrystal, 5 grains, rocprof-sys: $(tail -1 "$out/series.txt")"
    log "rocprof-sys: $(cd "$out/trace-sys" 2> /dev/null && find . -type f -printf '%p %s;' | head -c 900)"
  fi
fi
# tests of the polycrystal by ctest (polycrystal/CMakeLists.txt of
# mm-opera-hpc), with the behaviours of this build
if [ "${FM_POLYCRYSTAL_CTEST:-1}" = 1 ] && [ -d "$B/build/mm-opera-hpc-off" ]; then
  cd "$B/build/mm-opera-hpc-off" || exit 1
  if cmake --build . --target BehavioursPolycrystal -j 16 > "$out/build-behaviours.log" 2>&1; then
    (cd polycrystal && ctest -R '^UniaxialPolycrystal' -E "gpu" --output-on-failure) \
      > "$out/ctest-polycrystal.log" 2>&1
    log "polycrystal, ctest: $(grep -E 'tests passed|tests failed|No tests were found' "$out/ctest-polycrystal.log")" \
      "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout|Not Run)' "$out/ctest-polycrystal.log" | tr '\n' ' ')"
  else
    log "polycrystal: build of the behaviours failed (see $out/build-behaviours.log)"
  fi
fi
# other tests of mfem-mgis, on the host, for the variant off
others=${FM_OTHER_TESTS-}
if [ -n "$others" ]; then
  cd "$B/build/mfem-mgis-off" || exit 1
  targets=$(cmake --build . --target help | sed -n 's/^\.\.\. \(.*\)$/\1/p')
  behaviours="BehaviourTest $(grep -E '^GenericParallelBehaviour-[a-z]*$' <<< "$targets" | tr '\n' ' ')"
  # shellcheck disable=SC2086
  if cmake --build . --target $others $behaviours -j 32 > "$out/build-others.log" 2>&1; then
    ctest -R "^(${others// /|})" -E "cuda|hip|gpu|GPU" --output-on-failure \
      > "$out/ctest-others.log" 2>&1
    log "other tests: $(grep -E 'tests passed|tests failed' "$out/ctest-others.log")" \
      "$(grep -E '^ *[0-9]+ - .*(Failed|Exception|Timeout|Not Run)' "$out/ctest-others.log" | tr '\n' ' ')"
  else
    log "other tests: build failed (see $out/build-others.log)"
  fi
fi
rm -rf "$B"
log end
