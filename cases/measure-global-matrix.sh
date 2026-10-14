#!/bin/bash
# End to end measurements of the global matrix with the pattern and the
# exchanges kept (work package 7k): thermal case r5 and beam r4 on 8
# processes, tangent operator assembled by kernels, with the build of the
# cases linked to the mfem-mgis installed in install/mfem-mgis-jacobians
# (before, build/cases-jacobians) and in install/mfem-mgis-global (after,
# build/cases-global). To be run on a quiet machine, nothing else running on
# the cores (the processors N and N+8 share a core).
#
# usage: cases/measure-global-matrix.sh [results directory]
work_package=global
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/global-matrix-end-to-end-$(date +%Y-%m-%d-%H%M)}
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
for case_name in thermal beam; do
  source "$P/cases/measure-common.sh" "$case_name" "$out"
  if [ "$case_name" = thermal ]; then
    linear_solver=gmres-amg
    linear_tolerance=1e-8
  fi
  # measured versions, each one with its build directory
  # build/cases-<version>
  for version in jacobians global; do
    B=$P/build/cases-$version/$case_name
    run "$case_name-np8-kernels-$version" 1 8 "$generic" "${kernels[@]}"
  done
done
echo "results in $out"
