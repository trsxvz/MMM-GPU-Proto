#!/bin/bash
# End to end measurements of the assembly of the tangent operator by
# kernels (positions of the entries in the sparse matrix, matrices without
# forming B): thermal case r5 and beam r4 on 8 processes, tangent operator
# assembled by kernels, with the build of the cases linked to the mfem-mgis
# installed in install/mfem-mgis (before, build/cases) and to the one
# installed in install/mfem-mgis-registers (after, build/cases-registers),
# plus the computation point by point as reference. To be run on a quiet
# machine, nothing else running on the cores (the processors N and N+8
# share a core).
#
# usage: cases/measure-positions.sh [results directory]
work_package=positions
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/positions-$(date +%Y-%m-%d-%H%M)}
batch=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
for case_name in thermal beam; do
  source "$P/cases/measure-common.sh" "$case_name" "$out"
  if [ "$case_name" = thermal ]; then
    linear_solver=gmres-amg
    linear_tolerance=1e-8
  fi
  run "$case_name-np8-point-by-point" 1 8 "$generic"
  run "$case_name-np8-kernels-before" 1 8 "$generic" "${batch[@]}"
  B=$P/build/cases-registers/$case_name
  executable=./beam
  run "$case_name-np8-kernels-after" 1 8 "$generic" "${batch[@]}"
done
echo "results in $out"
