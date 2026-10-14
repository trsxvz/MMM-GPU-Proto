#!/bin/bash
# Overhead of each placement of the computations given by the option
# --execution of mfem-mgis (task 9a): beam r4 (111,843 unknowns), 10 steps,
# Newton, PCG of hypre and BoomerAMG, everything on the CPU
# (`CPU_PROCESSES` processes, 8 by default, on a quiet machine), everything
# on the GPU (one process), and the eight placements of the finite
# elements, of the behaviours and of the linear solvers in mixed mode (one
# process and the GPU; the beam has no post-processings). The behaviour is
# the one of the case, for the CPU (generic interface) or for the GPU
# (generic-parallel interface, cuda backend), depending on the placement of
# the behaviours. The copies of the arrays of the materials between the
# host and the GPU are reported below the profiling table.
#
# The variable `cases_build` gives the build of the cases, built against a
# mfem-mgis which has the option --execution (by default build/cases-9a).
#
# usage: cases/measure-execution-beam.sh [results directory]
work_package=execution
cases_build=${cases_build:-$(cd "$(dirname "$0")/.." && pwd)/build/cases-9a}
source "$(dirname "$0")/measure-common.sh" beam "$@"
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
linear_tolerance=1e-4
summary="^dofs|^time|host to device|device to host"
device=""
# everything on the CPU
options="--execution cpu"
np=${CPU_PROCESSES:-8}
run "execution-cpu-${np}p" 1 "$np" "$generic"
# everything on the GPU
options="--execution gpu"
run execution-gpu 1 1 "$cuda"
# mixed placements
for fe in cpu gpu; do
  for bh in cpu gpu; do
    for lsd in cpu gpu; do
      library=$generic
      if [ "$bh" = gpu ]; then
        library=$cuda
      fi
      options="--execution mixed --finite-elements $fe --behaviours $bh --linear-solver-device $lsd"
      run "execution-mixed-fe-$fe-behaviours-$bh-solver-$lsd" 1 1 "$library"
    done
  done
done
echo "results in $out"
