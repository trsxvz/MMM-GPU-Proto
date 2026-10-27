#!/bin/bash
# Strong scaling of the beam on the GPUs and on the cores of one node (on
# Jean Zay: 4 H100 and 96 cores per node of gpu_p6, job
# deps/machines/jean-zay-h100/jobs/measure-multi-gpu.slurm), with the best
# method of the work package L5 (NewtonK0, conjugate gradient preconditioned
# by the aggressive multigrid of K0, prediction, tolerance of
# Eisenstat-Walker 1, as cases/measure-beam-scaling.sh), all the kernels:
# - on 1, 2 and 4 GPUs, one process per GPU (each process takes the GPU of
#   its rank on the node), the exchanges of MPI through the host, then, on
#   more than one GPU, with the arrays of the GPUs given to MPI
#   (--gpu-aware-mpi). With srun (scripts/mpirun-srun), each process runs on
#   the NUMA node of its GPU (`--cpu-bind=map_ldom`, the NUMA node of each GPU
#   being read in `nvidia-smi topo -m`), all the GPUs staying visible to all
#   the processes: with one core per process, the second process ran on the
#   other socket than its GPU, and the copies from the GPU to the host at
#   the end of each step were 9 times slower; with `--gpu-bind=closest`,
#   which hides the GPUs of the other socket, the GPU-aware MPI of Open MPI
#   with PSM2 fails (cuIpcOpenMemHandle) (Jean Zay, 26/10/2026);
# - on 24, 48 and 96 processes on the CPU.
# The beam is refined `refinement` times (5 by default: 840,000 unknowns; 6:
# about 6.5 million). The lists of processes are given by `gpu_processes`
# and `cpu_processes`. Each computation gives the time of the loop over the
# steps and the profiling tables (cases/summary-profiling.py).
#
# usage: [refinement=6] cases/measure-multi-gpu.sh [results directory]
work_package=multi-gpu
source "$(dirname "$0")/measure-common.sh" beam "$@"
refinement=${refinement:-5}
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
linear_tolerance=1e-4
base="-nls fixed-point -fpm newton-k0 -k0 -pk0 -alst EisenstatWalker1 -alsti 0.1"
amg="-els amg -amgo AggressiveCoarseningLevels=1"
echo "refinement: $refinement" | tee -a "$out/summary.txt"
device=cuda
# NUMA node of each GPU (column "NUMA Affinity" of `nvidia-smi topo -m`,
# after the CPU affinity), empty if unknown
read -r -a gpu_numa <<< "$(nvidia-smi topo -m 2> /dev/null |
  awk -F'\t' '/^GPU[0-9]/ {
    for (i = 2; i < NF; ++i) {
      if (($i ~ /^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$/) &&
          ($(i + 1) ~ /^ *[0-9]+ *$/)) {
        printf "%d ", $(i + 1); break
      }
    }
  }')"
echo "NUMA nodes of the GPUs: ${gpu_numa[*]:-unknown}" | tee -a "$out/summary.txt"
for np in ${gpu_processes:-1 2 4}; do
  # process r on the NUMA node of the GPU r (the GPU of its rank, see
  # `getDeviceIdentifier` in mfem-mgis)
  unset PROTO_SRUN_OPTIONS
  if [ "${#gpu_numa[@]}" -ge "$np" ]; then
    ldoms=$(IFS=,; echo "${gpu_numa[*]:0:$np}")
    export PROTO_SRUN_OPTIONS="--cpu-bind=map_ldom:$ldoms"
  fi
  options="$base $amg"
  run "beam-r$refinement-gpu-np$np" 1 "$np" "$cuda" "$batch" "${kernels[@]}"
  if [ "$np" != 1 ]; then
    options="$base $amg --gpu-aware-mpi"
    run "beam-r$refinement-gpu-np$np-gpu-aware-mpi" 1 "$np" "$cuda" \
      "$batch" "${kernels[@]}"
  fi
done
unset PROTO_SRUN_OPTIONS
device=""
for np in ${cpu_processes:-24 48 96}; do
  options="$base $amg"
  run "beam-r$refinement-cpu-np$np" 1 "$np" "$generic" "$batch" \
    "${kernels[@]}"
done
echo "results in $out"
