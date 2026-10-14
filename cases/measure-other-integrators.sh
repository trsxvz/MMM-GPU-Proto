#!/bin/bash
# End to end measurements of the other integrators of the form added to the
# computations by kernels (work package 7m): thermal case r5 on 8
# processes, convective exchange at the second end (option -hc, boundary
# face integrator like the Robin boundary condition of the example ex8),
# everything by kernels, assembled tangent operator, with the build of the
# cases linked to the mfem-mgis installed in install/mfem-mgis-global
# (before: loops of MFEM over all the elements, build/cases-global) and in
# install/mfem-mgis-others (after, build/cases-others), and, for comparison,
# the same computation with the imposed temperature (without other
# integrator). To be run on a quiet machine, nothing else running on the
# cores (the processors N and N+8 share a core).
#
# usage: cases/measure-other-integrators.sh [results directory]
work_package=others
P=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-$P/results/other-integrators-$(date +%Y-%m-%d-%H%M)}
source "$P/cases/measure-common.sh" thermal "$out"
# measured versions of mfem-mgis, each one with its build directories
# build/cases-<version> and build/mfem-mgis-<version>
versions="global others"
# the behaviour is the one of build/cases/thermal for both versions (it does
# not depend on mfem-mgis)
for version in $versions; do
  if [ ! -x "$P/build/cases-$version/thermal/thermal" ]; then
    echo "build/cases-$version/thermal/thermal missing" >&2
    exit 1
  fi
  source_dir=$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' \
    "$P/build/mfem-mgis-$version/CMakeCache.txt")
  commit=$(git -C "$source_dir" log -1 --format=%h 2> /dev/null ||
    echo "unknown, $source_dir missing")
  echo "mfem-mgis-$version: $commit" | tee -a "$out/summary.txt"
done
kernels=(MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
linear_solver=gmres-amg
linear_tolerance=1e-8
for version in $versions; do
  B=$P/build/cases-$version/thermal
  options="-hc 1"
  run "thermal-np8-robin-$version" 1 8 "$generic" "${kernels[@]}"
  options=""
  run "thermal-np8-dirichlet-$version" 1 8 "$generic" "${kernels[@]}"
done
echo "results in $out"
