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
# With `counters=1`, the computations on the GPU only are made under Nsight
# Compute (`ncu`, whose hardware counters must be allowed to the user), which
# measures the first `counters_launches` (60 by default) kernels of mfem-mgis
# and of the integration of the behaviour, with the caches flushed before
# each kernel: sectors per request of the global loads and stores (8 for
# coalesced accesses to double values, a warp reading 256 bytes, that is 8
# sectors of 32 bytes; 32 if the values of consecutive threads are 32 bytes
# or more apart), bytes read and written in the DRAM, duration, throughput
# of the DRAM and of the SMs relative to their peak. Those counters do not
# depend on the size of the cache of the GPU, unlike the times. Reports in
# `<name>.ncu-rep` (not versioned), values in `<name>-counters.csv`, summary
# by kernel and by layout by `cases/summary-counters.py <directory>`.
#
# The variable `executable` gives the executable of the case, by default the
# one of the usual builds (build/cases/beam, build/mm-opera-hpc/polycrystal).
#
# usage: [counters=1] cases/measure-layout.sh [beam|polycrystal]
#        [results directory]
work_package=layout
case_name=${1:-beam}
case "$case_name" in
  beam | polycrystal) ;;
  *)
    echo "usage: $0 [beam|polycrystal] [results directory]" >&2
    exit 1
    ;;
esac
source "$(dirname "$0")/measure-common.sh" "$case_name" "${2:-}"
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_tolerance=1e-4
np=${CPU_PROCESSES:-8}
counters=${counters:-0}
# metrics of Nsight Compute (counters=1)
metrics=gpu__time_duration.sum
metrics+=,l1tex__t_requests_pipe_lsu_mem_global_op_ld.sum
metrics+=,l1tex__t_sectors_pipe_lsu_mem_global_op_ld.sum
metrics+=,l1tex__t_requests_pipe_lsu_mem_global_op_st.sum
metrics+=,l1tex__t_sectors_pipe_lsu_mem_global_op_st.sum
# bytes read and written in the DRAM: dram__bytes_read and dram__bytes_write
# on the A100 and H100 GPUs, dram__bytes_op_read and dram__bytes_op_write on
# the GPUs of the Blackwell generation; their sum exists on all of them
metrics+=,dram__bytes.sum
metrics+=,dram__throughput.avg.pct_of_peak_sustained_elapsed
metrics+=,sm__throughput.avg.pct_of_peak_sustained_elapsed
if [ "$counters" = 1 ]; then
  ncu=${NCU:-$(command -v ncu || true)}
  if [ -z "$ncu" ]; then
    echo "ncu not found (variable NCU or PATH)" | tee -a "$out/summary.txt"
    exit 1
  fi
fi
# usage: run_gpu name [VARIABLE=value...]
# computation on the GPU, under Nsight Compute if counters=1
run_gpu() {
  local name=$1
  shift
  if [ "$counters" = 1 ]; then
    launcher=("$ncu" --metrics "$metrics" --cache-control all
      --clock-control base --kernel-name-base demangled
      --kernel-name "regex:mfem_mgis::batch|generic_parallel::gpu"
      --launch-count "${counters_launches:-60}" -f -o "$out/$name")
  fi
  run "$name" 1 1 "$cuda" "$@"
  launcher=()
  if [ "$counters" = 1 ] && [ -f "$out/$name.ncu-rep" ]; then
    "$ncu" --import "$out/$name.ncu-rep" --csv --page raw \
      --metrics "$metrics" > "$out/$name-counters.csv"
  fi
}
for layout in PointByPoint ComponentByComponent; do
  settings=(MFEM_MGIS_QUADRATURE_DATA_LAYOUT=$layout
    MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered)
  # GPU, assembled tangent operator
  options="--execution gpu"
  linear_solver=pcg-amg
  reference="layout-gpu-assembled-PointByPoint"
  run_gpu "layout-gpu-assembled-$layout" "${settings[@]}"
  # GPU, tangent operator applied without assembly (beam)
  if [ "$case_name" = beam ]; then
    linear_solver=cg-jacobi
    run_gpu "layout-gpu-matrix-free-$layout" "${settings[@]}" \
      MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree
  fi
  # CPU
  if [ "$counters" != 1 ]; then
    options="--execution cpu"
    linear_solver=pcg-amg
    reference="layout-cpu-${np}p-PointByPoint"
    run "layout-cpu-${np}p-$layout" 1 "$np" "$generic" "${settings[@]}"
  fi
done
if [ "$counters" = 1 ]; then
  python3 "$P/cases/summary-counters.py" "$out" | tee -a "$out/summary.txt"
fi
echo "results in $out"
