#!/bin/bash
# Strong scaling of the beam (or of the polycrystal of mm-opera-hpc, with
# multi_gpu_case=polycrystal) on the GPUs and on the cores of one node (on
# Jean Zay: 4 H100 and 96 cores per node of gpu_p6, job
# deps/machines/jean-zay-h100/jobs/measure-multi-gpu.slurm), with the best
# method of the work package L5 (NewtonK0, conjugate gradient preconditioned
# by the aggressive multigrid of K0, prediction, tolerance of
# Eisenstat-Walker 1, as cases/measure-beam-scaling.sh), all the kernels:
# - on 1, 2 and 4 GPUs, one process per GPU (each process takes the GPU of
#   its rank on the node), the exchanges of MPI through the host, then, on
#   more than one GPU, with the arrays of the GPUs given to MPI
#   (--gpu-aware-mpi). The processes are spread over the cores of the node
#   (cores / processes each, the cores being the largest number of
#   processes on the CPU), and the GPUs are ordered in CUDA_VISIBLE_DEVICES
#   (PROTO_CUDA_VISIBLE_DEVICES, set in each process by scripts/mpirun-srun,
#   srun replacing CUDA_VISIBLE_DEVICES) so that the GPU of each rank is on
#   the NUMA node, or at least on the socket, of its cores
#   (scripts/gpus-close-to-ranks.py), all the GPUs staying visible to all the
#   processes. With one core per process, the
#   second process ran on the other socket than its GPU, and the copies from
#   the GPU to the host at the end of each step were 9 times slower; with
#   `--gpu-bind=closest`, which hides the GPUs of the other socket, the
#   GPU-aware MPI of Open MPI with PSM2 fails (cuIpcOpenMemHandle); and srun
#   refuses `--cpu-bind=map_ldom` in a step which does not use the whole node
#   (Jean Zay, 26-27/10/2026);
# - on 24, 48 and 96 processes on the CPU.
# The beam is refined `refinement` times (5 by default: 840,000 unknowns; 6:
# about 6.5 million). The polycrystal (5 grains, single crystal behaviour of
# UO2, periodic, uniaxial tension, 10 steps; NewtonK0 and Eisenstat-Walker 1,
# the best method of results/L5-polycrystal-2026-10-10-0502) is refined
# `refinement` times (2 by default: about 1.3 million unknowns; 0: 20,553);
# its macroscopic stresses are compared to the ones of the computation on one
# GPU. The lists of processes are given by `gpu_processes` and
# `cpu_processes`. Each computation gives the profiling tables
# (cases/summary-profiling.py), and, for the beam, the time of the loop over
# the steps.
#
# usage: [multi_gpu_case=beam|polycrystal] [refinement=n] \
#   cases/measure-multi-gpu.sh [results directory]
work_package=multi-gpu
source "$(dirname "$0")/measure-common.sh" "${multi_gpu_case:-beam}" "$@"
if [ "$case_name" = polycrystal ]; then
  refinement=${refinement:-2}
  # number of unknowns of the polycrystal
  summary="$summary|number of finite element unknowns"
else
  refinement=${refinement:-5}
fi
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
linear_solver=pcg-amg
linear_tolerance=1e-4
base="-nls fixed-point -fpm newton-k0 -k0 -pk0 -alst EisenstatWalker1 -alsti 0.1"
amg="-els amg -amgo AggressiveCoarseningLevels=1"
if [ "$case_name" = polycrystal ]; then
  base="--non-linear-solver NewtonK0"
  base="$base --adaptive-linear-solver-tolerance EisenstatWalker1"
  base="$base --adaptive-linear-solver-initial-tolerance 0.1 -r $refinement"
  amg=""
  reference=polycrystal-r$refinement-gpu-np1
fi
echo "refinement: $refinement" | tee -a "$out/summary.txt"
device=cuda
cpu_list=(${cpu_processes:-24 48 96})
cores=${cpu_list[-1]}
for np in ${gpu_processes:-1 2 4}; do
  threads=$((cores / np))
  # GPU of the rank r close to its cores (cores r * threads and above)
  order=$(python3 -I "$P/scripts/gpus-close-to-ranks.py" "$np" "$threads")
  if [ -n "$order" ]; then
    export PROTO_CUDA_VISIBLE_DEVICES=$order
  fi
  echo "$np processes: $threads cores each, CUDA_VISIBLE_DEVICES=${order:-unchanged}" |
    tee -a "$out/summary.txt"
  options="$base $amg"
  run "$case_name-r$refinement-gpu-np$np" "$threads" "$np" "$cuda" "$batch" \
    "${kernels[@]}"
  if [ "$np" != 1 ]; then
    options="$base $amg --gpu-aware-mpi"
    run "$case_name-r$refinement-gpu-np$np-gpu-aware-mpi" "$threads" "$np" \
      "$cuda" "$batch" "${kernels[@]}"
  fi
done
unset PROTO_CUDA_VISIBLE_DEVICES
device=""
for np in ${cpu_processes:-24 48 96}; do
  options="$base $amg"
  run "$case_name-r$refinement-cpu-np$np" 1 "$np" "$generic" "$batch" \
    "${kernels[@]}"
done
echo "results in $out"
