#!/bin/bash
# polycrystal of mm-opera-hpc (5 grains, periodic, uniaxial tension,
# 20 553 unknowns, 10 steps), all by kernels, one process on the GPU (host
# cores given by CORES, 0-3 by default), cuda law: Newton with HyprePCG and
# BoomerAMG, NewtonK0 with the tolerance of Eisenstat and Walker (best
# variant measured on 10/10), Newton with cuDSS (analysis reused), hypre on
# the GPU or on the host; two passes; macroscopic stresses in *.res
P=$(cd "$(dirname "$0")/../.." && pwd)
source $P/scripts/env.sh > /dev/null
out=$P/results/cudss-polycrystal-2026-10-12
L=$P/build/cases/polycrystal/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
cd $P/build/mm-opera-hpc-cudss/polycrystal
export OMP_NUM_THREADS=1 MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch \
  MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels \
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels \
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
ew="--adaptive-linear-solver-tolerance EisenstatWalker1"
ew="$ew --adaptive-linear-solver-initial-tolerance 0.1"
for pass in 1 2; do
  for variant in newton newton-k0-ew1 newton-cudss newton-cudss-hypre-host; do
    case $variant in
      newton) extra="" ;;
      newton-k0-ew1) extra="--non-linear-solver NewtonK0 $ew" ;;
      newton-cudss) extra="--linear-solver CuDSSSolver --reuse-analysis" ;;
      newton-cudss-hypre-host)
        extra="--linear-solver CuDSSSolver --reuse-analysis --hypre-on-host" ;;
    esac
    echo "== pass $pass, $variant, load $(cut -d ' ' -f 1 /proc/loadavg)"
    taskset -c ${CORES:-0-3} timeout 900 mpirun -n 1 ./uniaxial-polycrystal \
      --disable-post-processing --nstep 10 --duration 10 -l $L \
      --macroscopic-stress-output-file $out/polycrystal-gpu-$variant-$pass.res \
      --device cuda $extra > $out/polycrystal-gpu-$variant-$pass.txt 2>&1
    echo "   exit $?"
  done
done
