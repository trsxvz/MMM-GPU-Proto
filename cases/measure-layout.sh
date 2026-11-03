#!/bin/bash
# Layout of the data at the integration points (parameter
# QuadratureDataLayout of mfem-mgis, given here by the environment variable
# MFEM_MGIS_QUADRATURE_DATA_LAYOUT): each computation is made with the
# arrays of the material stored point by point (PointByPoint) and component
# by component (ComponentByComponent):
#
# - beam: beam r4 (111,843 unknowns, small strain, tangent operator blocks
#   of 36 components), 10 steps, Newton, on the GPU with the tangent
#   operator assembled by kernels (PCG of hypre and BoomerAMG) or applied
#   without assembly (conjugate gradient of MFEM and Jacobi), and on the
#   CPU (PCG of hypre and BoomerAMG);
# - polycrystal: polycrystal of mm-opera-hpc (finite strain, tangent
#   operator blocks of 81 components), 10 steps, on the GPU and on the CPU,
#   the macroscopic stresses of the component by component layout being
#   compared to the ones of the point by point layout.
#
# On the GPU: `--execution gpu`, one process (behaviour integrated by
# batches on the GPU, generic-parallel interface, cuda backend, finite
# elements computed by kernels). On the CPU: `--execution cpu`,
# `CPU_PROCESSES` processes (8 by default, on a quiet machine). The
# stiffness matrices are assembled in a fixed order
# (MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered). The time of each phase is
# given by `cases/summary-profiling.py <directory>/*.txt`.
#
# The variable `executable` gives the executable of the case, built against
# a mfem-mgis which has the parameter QuadratureDataLayout: by default
# build/cases-layout/beam/beam and
# build/mm-opera-hpc-layout/polycrystal/uniaxial-polycrystal (see
# cases/README.md). The behaviours are the ones of the usual builds.
#
# usage: cases/measure-layout.sh [beam|polycrystal] [results directory]
work_package=layout
case_name=${1:-beam}
proto=$(cd "$(dirname "$0")/.." && pwd)
case "$case_name" in
  beam)
    executable=${executable:-$proto/build/cases-layout/beam/beam}
    ;;
  polycrystal)
    executable=${executable:-$proto/build/mm-opera-hpc-layout/polycrystal/uniaxial-polycrystal}
    ;;
  *)
    echo "usage: $0 [beam|polycrystal] [results directory]" >&2
    exit 1
    ;;
esac
source "$(dirname "$0")/measure-common.sh" "$case_name" "${2:-}"
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_tolerance=1e-4
np=${CPU_PROCESSES:-8}
for layout in PointByPoint ComponentByComponent; do
  settings=(MFEM_MGIS_QUADRATURE_DATA_LAYOUT=$layout
    MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered)
  # GPU, assembled tangent operator
  options="--execution gpu"
  linear_solver=pcg-amg
  reference="layout-gpu-assembled-PointByPoint"
  run "layout-gpu-assembled-$layout" 1 1 "$cuda" "${settings[@]}"
  # GPU, tangent operator applied without assembly (beam)
  if [ "$case_name" = beam ]; then
    linear_solver=cg-jacobi
    run "layout-gpu-matrix-free-$layout" 1 1 "$cuda" "${settings[@]}" \
      MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree
  fi
  # CPU
  options="--execution cpu"
  linear_solver=pcg-amg
  reference="layout-cpu-${np}p-PointByPoint"
  run "layout-cpu-${np}p-$layout" 1 "$np" "$generic" "${settings[@]}"
done
echo "results in $out"
