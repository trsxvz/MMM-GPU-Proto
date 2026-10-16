#!/bin/bash
# ex8 of mfem-mgis-examples (thermomechanical plate, 5 steps): mechanics
# solved by MUMPS against the conjugate gradient of MFEM preconditioned by
# MueLu (rigid body modes), on 2 CPU processes and on the GPU (one
# process), everything by kernels, then the mesh refined once (5 steps on
# the GPU, first step on the GPU and on 2 CPU processes). The example must be built against an
# installation of mfem-mgis built with Trilinos (for instance
# install/mfem-mgis-trilinos, build build/examples-trilinos-ex8).
#
# usage: measure-ex8-muelu.sh [build of ex8] [results directory]
# Environment: RUNS (regular expression selecting the runs), RUN_TIMEOUT
# (seconds, 900 by default).
set -u
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
B=$(cd "${1:-$P/build/examples-trilinos-ex8}" && pwd)
out=${2:-$P/results/ex8-muelu-$(date +%Y-%m-%d-%H%M)}
if [ ! -x "$B/rjh_plate" ]; then
  echo "no example ex8 in $B" >&2
  exit 1
fi
mkdir -p "$out"
out=$(cd "$out" && pwd)
MPI=$(dirname "$(dirname "$(/usr/bin/grep -m1 '^MPI_CXX_COMPILER:' "$B/CMakeCache.txt" | sed 's/.*=//')")")
CUDA=${PROTO_CUDA:-$HOME/Bureau/spack-mmm-externals/cuda-13.0.2}
export PATH=$CUDA/bin:$MPI/bin:/usr/local/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=$P/install/h-interface/tfel/lib:$P/install/h-interface/mgis/lib:$CUDA/lib64
export OMPI_MCA_btl=self,vader
runs=${RUNS:-.}
env=(-x OMP_NUM_THREADS=1 -x MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
  -x MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  -x MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  -x MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels
  -x LD_LIBRARY_PATH -x OMPI_MCA_btl)
{
  echo "date: $(date '+%F %T')"
  echo "example: ${B#"$P"/}"
  echo "mfem-mgis: $(/usr/bin/grep -m1 '^MFEMMGIS_DIR:' "$B/CMakeCache.txt" | sed "s|.*=||; s|$P/||")"
  nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null
  LC_ALL=C lscpu | /usr/bin/grep "Model name"
} >> "$out/context.txt"
# name, number of processes, options of the example
run() {
  local name=$1 np=$2
  shift 2
  echo "$name" | /usr/bin/grep -Eq "$runs" || return 0
  (cd "$B" && timeout "${RUN_TIMEOUT:-900}" taskset -c 0,1 mpirun -np "$np" \
    --bind-to none "${env[@]}" ./rjh_plate "$@") < /dev/null \
    > "$out/$name.txt" 2>&1
  echo "$name: status $? $(/usr/bin/grep -a -E -- '--> (IterativeCouplingScheme::computeNextState|Thermal|Mechanics) ' "$out/$name.txt" | awk -F'|' '{gsub(/ +/," ",$3); gsub(/ +/,"",$5); print $3 $5}' | tr '\n' ' ')" | tee -a "$out/summary.txt"
}
common=(-m assemblage_hexa.msh -no-pp -ns 5)
muelu=(-lsMc CGSolver -pcMc MueLu)
# coarse mesh compared to the reference statistics of the example
coarse=(-m assemblage_hexa_coarse.msh -no-pp
  -rf assemblage_hexa_coarse-statistics.ref)
run check-coarse-cg-muelu 2 "${coarse[@]}" "${muelu[@]}"
run check-coarse-chebyshev 2 "${coarse[@]}" "${muelu[@]}" -mls chebyshev
run check-coarse-minres-chebyshev 2 "${coarse[@]}" -lsMc MINRESSolver \
  -pcMc MueLu -mls chebyshev
run check-coarse-jacobi-0.3 2 "${coarse[@]}" "${muelu[@]}" -mls jacobi -mld 0.3
run cpu-np2-mumps 2 "${common[@]}" -lsMc MUMPSSolver
run cpu-np2-cg-muelu 2 "${common[@]}" "${muelu[@]}"
run cpu-np2-cg-muelu-chebyshev 2 "${common[@]}" "${muelu[@]}" -mls chebyshev
run gpu-np1-mumps-spd 1 "${common[@]}" -lsMc MUMPSSolver \
  -msMc SymmetricPositiveDefinite --execution mixed --linear-solver-device cpu
run gpu-np1-cg-muelu 1 "${common[@]}" "${muelu[@]}" --execution gpu
run gpu-np1-cg-muelu-chebyshev 1 "${common[@]}" "${muelu[@]}" \
  -mls chebyshev --execution gpu
# the Jacobi smoother damped by 0.6 (default) and 0.3: not positive definite
run cpu-np2-cg-muelu-jacobi 2 "${common[@]}" "${muelu[@]}" -mls jacobi
run gpu-np1-cg-muelu-jacobi 1 "${common[@]}" "${muelu[@]}" -mls jacobi \
  --execution gpu
run gpu-np1-cg-muelu-jacobi-0.3 1 "${common[@]}" "${muelu[@]}" -mls jacobi \
  -mld 0.3 --execution gpu
r1=(-m assemblage_hexa.msh -r 1 -no-pp -ns 5 -pcMc MueLu -mls chebyshev
  --execution gpu)
run gpu-np1-r1-cg-muelu-chebyshev 1 "${r1[@]}" -lsMc CGSolver
run gpu-np1-r1-minres-muelu-chebyshev 1 "${r1[@]}" -lsMc MINRESSolver
run gpu-np1-r1-gmres-muelu-chebyshev 1 "${r1[@]}" -lsMc GMRESSolver
# first step of the refined mesh (same step as the first of the 5 steps)
step1=(-m assemblage_hexa.msh -r 1 -no-pp -et 2e4 -tr 1e5 -ns 1)
run gpu-np1-r1-step1-cg-muelu-chebyshev 1 "${step1[@]}" "${muelu[@]}" \
  -mls chebyshev --execution gpu
run cpu-np2-r1-step1-cg-muelu 2 "${step1[@]}" "${muelu[@]}"
# iterations of the conjugate gradient (verbosity 1)
run cpu-np2-cg-muelu-v1 2 "${common[@]}" "${muelu[@]}" -v 1
