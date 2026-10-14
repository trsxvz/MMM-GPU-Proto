#!/bin/bash
# Pressure imposed on the end of the beam (work package L7.6). With the
# inner forces computed by kernels, the pressure, independent of the
# unknowns, is assembled once for a unit pressure and added to the residual
# on the GPU, without integrator in the non linear form: the tangent
# operator is assembled by kernels without going through the assembly of
# MFEM, and can be applied without assembly. Comparison with an executable
# built with the previous version of mfem-mgis (first argument), where the
# pressure is a boundary integrator evaluated on the host at each residual
# and where the element matrices computed by kernels are assembled by MFEM
# on the host. Beam r4, 10 steps, NewtonK0 (best method of the work
# package L5), all the kernels, on the GPU and on 8 CPU processes, with and
# without pressure.
# Before each computation on the GPU, the script waits at most 10 minutes
# for another process to free the GPU, then reports it in the summary.
#
# The variable AFTER gives the executable of the new version (by default
# the one of the build directory of the case).
#
# usage: cases/measure-pressure-beam.sh executable-before [results directory]
work_package=pressure
before=$(realpath "$1")
shift
source "$(dirname "$0")/measure-common.sh" beam "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
matrix_free=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
linear_tolerance=1e-4
pressure=${PRESSURE:-1e8}
base="-nls fixed-point -fpm newton-k0 -k0 -pk0 -alst EisenstatWalker1 -alsti 0.1 -els amg -amgo AggressiveCoarseningLevels=1"
echo "pressure: $pressure Pa, executable before: $before" | tee -a "$out/summary.txt"

# waits at most 10 minutes for the GPU to be free of other processes
wait_for_gpu() {
  local i
  for ((i = 0; i != 120; ++i)); do
    if [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; then
      return 0
    fi
    sleep 5
  done
  echo "   GPU shared with another process" | tee -a "$out/summary.txt"
}

for version in before after; do
  if [ "$version" = before ]; then
    executable=$before
  else
    executable=${AFTER:-./beam}
  fi
  device=cuda
  options="$base -hdps 0.5"
  wait_for_gpu
  run "beam-gpu-without-pressure-$version" 1 1 "$cuda" "$batch" "${kernels[@]}"
  options="$base -hdps 0.5 -pa $pressure"
  wait_for_gpu
  run "beam-gpu-pressure-$version" 1 1 "$cuda" "$batch" "${kernels[@]}"
  if [ "$version" = after ]; then
    wait_for_gpu
    run "beam-gpu-pressure-matrix-free-$version" 1 1 "$cuda" "$batch" \
      "${matrix_free[@]}"
  fi
  device=""
  options="$base"
  run "beam-np8-without-pressure-$version" 1 8 "$generic" "$batch" \
    "${kernels[@]}"
  options="$base -pa $pressure"
  run "beam-np8-pressure-$version" 1 8 "$generic" "$batch" "${kernels[@]}"
done
python3 "$P/cases/summary-profiling.py" "$out"/beam-*.txt > "$out/profiling.txt" 2>&1
echo "end: $(date '+%F %T')" | tee -a "$out/summary.txt"
