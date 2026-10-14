#!/bin/bash
# Measurements of the work package L2 (behaviours integrated by batches on
# the CPU): integration point by point (library `generic`) against
# integration by batches (libraries `generic` and `generic-parallel`,
# backends `openmp` and `stlpar`), on 1 and 8 processes, then with several
# threads per process. To be run on a quiet machine, after
# build-mfem-mgis.sh.
#
# usage: cases/measure-L2.sh [beam|polycrystal] [results directory]
# (see cases/measure-common.sh)
work_package=L2
source "$(dirname "$0")/measure-common.sh" "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
for np in 8 1; do
  run "$case_name-np$np-point-by-point" 1 "$np" "$generic"
  run "$case_name-np$np-batches-generic" 1 "$np" "$generic" "$batch"
  run "$case_name-np$np-batches-openmp-1thread" 1 "$np" "$openmp" "$batch"
done
# one process on the 8 physical cores (the stlpar backend of g++ uses TBB,
# which takes all the hardware threads of the reserved cores)
run "$case_name-np1-batches-openmp-8threads" 8 1 "$openmp" "$batch"
run "$case_name-np1-batches-stlpar-8cores" 8 1 "$stlpar" "$batch"
# two processes of four cores (and the reference point by point on two
# processes of one core)
run "$case_name-np2-point-by-point" 1 2 "$generic"
run "$case_name-np2-batches-openmp-4threads" 4 2 "$openmp" "$batch"
echo "results in $out"
