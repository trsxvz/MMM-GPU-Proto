#!/bin/bash
# Variants of the options of MueLu on the GPU, beam of 421,443 unknowns
# (SIZE=32), one process: setup and resolution of a second operator.
# usage: measure-muelu-options.sh [build directory] [output file]
#   build directory: by default build/ginkgo/mfem-mgis-trilinos-host;
#   output file: one line per run (key=value, see LinearSolversBenchmark),
#   to be summarized by cases/summary-solvers.py.
set -u
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
B=${1:-$P/build/ginkgo/mfem-mgis-trilinos-host}
out=${2:-muelu-options.txt}
MPI=$(dirname "$(dirname "$(/usr/bin/grep -m1 '^MPI_CXX_COMPILER:' "$B/CMakeCache.txt" | sed 's/.*=//')")")
CUDA=${PROTO_CUDA:-$HOME/Bureau/spack-mmm-externals/cuda-13.0.2}
export PATH=$CUDA/bin:$MPI/bin:/usr/local/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=$P/install/h-interface/tfel/lib:$P/install/h-interface/mgis/lib:$CUDA/lib64
export OMPI_MCA_btl=self,vader OMP_NUM_THREADS=1
unset OMP_PROC_BIND OMP_PLACES
it='RelativeTolerance=1e-8,AbsoluteTolerance=0.,MaximumNumberOfIterations=2000'
jacobi='"smoother: type"=RELAXATION,"smoother: params"={"relaxation: type"=Jacobi,"relaxation: sweeps"=2,"relaxation: damping factor"=0.6}'
variants=(
  "jacobi|$jacobi"
  "jacobi-triple-product|$jacobi,\"rap: triple product\"=true"
  "jacobi-implicit-transpose|$jacobi,\"transpose: use implicit\"=true"
  "jacobi-unfiltered|$jacobi,\"sa: use filtered matrix\"=false"
  "jacobi-drop-0.02|$jacobi,\"aggregation: drop tol\"=0.02"
  "jacobi-coarse-10000|$jacobi,\"coarse: max size\"=10000"
  "jacobi-one-sweep|\"smoother: type\"=RELAXATION,\"smoother: params\"={\"relaxation: type\"=Jacobi,\"relaxation: sweeps\"=1,\"relaxation: damping factor\"=0.6}"
  "chebyshev-2|\"smoother: type\"=CHEBYSHEV,\"smoother: params\"={\"chebyshev: degree\"=2}"
  "chebyshev-2-eigen-30|\"smoother: type\"=CHEBYSHEV,\"smoother: params\"={\"chebyshev: degree\"=2,\"chebyshev: eigenvalue max iterations\"=30,\"chebyshev: boost factor\"=1.2}"
)
: > "$out"
for problem in elasticity heterogeneous; do
  case $problem in
    elasticity) args="-p elasticity -n 32" ;;
    heterogeneous) args="-p elasticity -n 32 -c 1000 -i 4" ;;
  esac
  for v in "${variants[@]}"; do
    IFS='|' read -r label options <<< "$v"
    # shellcheck disable=SC2086
    timeout 300 taskset -c 0,1 mpirun -np 1 --bind-to none -x LD_LIBRARY_PATH \
      -x OMPI_MCA_btl -x OMP_NUM_THREADS "$B/benchmarks/LinearSolversBenchmark" $args \
      -s CGSolver -sp "{$it,Preconditioner={Name=MueLu,Options={Parameters={$options}}}}" \
      -lb "$problem/muelu-$label/gpu" --device cuda 2>&1 | /usr/bin/grep -a "^label" >> "$out"
  done
done
