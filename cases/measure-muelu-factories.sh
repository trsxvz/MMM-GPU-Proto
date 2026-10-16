#!/bin/bash
# MueLu with the rigid body modes: which factory written with Kokkos (the
# default of MueLu on all the nodes of Tpetra but the serial one) degrades
# the coarse space? Linear elasticity of the benchmark LinearSolversBenchmark
# of mfem-mgis, conjugate gradient of MFEM preconditioned by MueLu with a
# Jacobi smoother (two sweeps damped by 0.6), one process, on the GPU
# (SIZE_GPU, 32 by default: 421,443 unknowns) then on the host with the
# serial node of Tpetra (SIZE_HOST, 16 by default: 56,355 unknowns). Each
# variant adds parameters of MueLu: factories written with Kokkos or not
# (`use kokkos refactor`), aggregation on the host or with Kokkos
# (`aggregation: backend`), local QR of the tentative prolongator
# (`tentative: calculate qr`), filtered matrix of the prolongator smoothing.
#
# usage: measure-muelu-factories.sh [build of mfem-mgis] [results directory]
#   build: by default build/ginkgo/mfem-mgis-trilinos-host (mfem-mgis built
#   with Trilinos); results directory: by default
#   results/muelu-factories-<date>.
# Environment: VERBOSITY (verbosity of MueLu, low by default; medium prints
# the aggregates and the levels), FILTER (regular expression selecting the
# variants), PREFIX (prefix of the names of the logs).
set -u
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build=$(cd "${1:-$P/build/ginkgo/mfem-mgis-trilinos-host}" && pwd)
out=${2:-$P/results/muelu-factories-$(date +%Y-%m-%d-%H%M)}
B=$build/benchmarks/LinearSolversBenchmark
if [ ! -x "$B" ]; then
  echo "no benchmark $B" >&2
  exit 1
fi
mkdir -p "$out"
CUDA=${PROTO_CUDA:-$HOME/Bureau/spack-mmm-externals/cuda-13.0.2}
export LD_LIBRARY_PATH=$P/install/h-interface/tfel/lib:$P/install/h-interface/mgis/lib:$CUDA/lib64
export OMPI_MCA_btl=self,vader
export OMP_NUM_THREADS=1
unset OMP_PROC_BIND OMP_PLACES
verbosity=${VERBOSITY:-low}
filter=${FILTER:-.}
prefix=${PREFIX:-}
it='RelativeTolerance=1e-8,AbsoluteTolerance=0.,MaximumNumberOfIterations=2000'
smoother='"smoother: type"=RELAXATION,"smoother: params"={"relaxation: type"=Jacobi,"relaxation: sweeps"=2,"relaxation: damping factor"=0.6}'
kokkos='"use kokkos refactor"=true'
classic='"use kokkos refactor"=false'
# name:parameters of MueLu added to the smoother
variants=(
  "classic:$classic"
  "kokkos:$kokkos"
  "kokkos-host-aggregation:$kokkos,\"aggregation: backend\"=host"
  "classic-kokkos-aggregation:$classic,\"aggregation: backend\"=kokkos"
  "kokkos-random:$kokkos,\"aggregation: deterministic\"=false"
  "kokkos-phase1-serial:$kokkos,\"aggregation: phase 1 algorithm\"=Serial"
  "classic-droptol:$classic,\"aggregation: drop tol\"=1e-14"
  "classic-noqr:$classic,\"tentative: calculate qr\"=false"
  "kokkos-host-aggregation-noqr:$kokkos,\"aggregation: backend\"=host,\"tentative: calculate qr\"=false"
  "kokkos-unfiltered:$kokkos,\"sa: use filtered matrix\"=false"
  "kokkos-no-lumping:$kokkos,\"filtered matrix: use lumping\"=false"
  "classic-unfiltered:$classic,\"sa: use filtered matrix\"=false"
  "kokkos-host-aggregation-unfiltered:$kokkos,\"aggregation: backend\"=host,\"sa: use filtered matrix\"=false"
  "classic-kokkos-aggregation-unfiltered:$classic,\"aggregation: backend\"=kokkos,\"sa: use filtered matrix\"=false")
{
  echo "date: $(date '+%F %T')"
  echo "build: ${build#"$P"/}"
  nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null
  LC_ALL=C lscpu | /usr/bin/grep "Model name"
} >> "$out/context.txt"
for target in gpu host; do
  if [ $target = gpu ]; then
    n=${SIZE_GPU:-32}
    device=(--device cuda)
  else
    n=${SIZE_HOST:-16}
    device=()
  fi
  for v in "${variants[@]}"; do
    name=${v%%:*}
    extra=${v#*:}
    echo "$name" | /usr/bin/grep -Eq "$filter" || continue
    log=$out/$prefix$target-$name.log
    timeout 300 taskset -c 0,1 "$B" -p elasticity -n "$n" -s CGSolver \
      -sp "{$it,Preconditioner={Name=MueLu,Options={Parameters={$smoother,$extra,verbosity=$verbosity}}}}" \
      -lb "elasticity/muelu-factories/$name/$target" "${device[@]}" \
      < /dev/null > "$log" 2>&1
    echo "$target $name: status $? $(/usr/bin/grep -a -o 'iterations=[0-9]*\|residual=[0-9.e+-]*\|converged=[01]' "$log" | tr '\n' ' ')" | tee -a "$out/summary.txt"
  done
done
