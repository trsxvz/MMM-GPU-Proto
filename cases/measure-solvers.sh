#!/bin/bash
# Campaign of measurements of the linear solvers of mfem-mgis (phase 5 of
# docs/plan-solver-libraries.md): each configuration of the list below is
# run by the benchmark LinearSolversBenchmark of mfem-mgis (benchmarks/,
# CMake option BENCHMARKS) on assembled systems: the linear elasticity of a
# cantilever beam, homogeneous or with stiff cubic inclusions (contrast
# 1000), and the Laplace equation on a cube. Each configuration is run on
# the CPU (8 processes, one per core) and, for the entries which compute on
# the GPU, on the GPU (one process, device cuda of MFEM).
#
# usage: measure-solvers.sh [build directory of mfem-mgis] [results directory]
#   build directory: by default build/ginkgo/mfem-mgis-trilinos (mfem-mgis
#   of the environment mmm-trilinos, with all the libraries);
#   results directory: by default results/solvers-<date>.
# Environment: SIZE (elements across the section, 16 by default), PROBLEMS
# (default "elasticity heterogeneous laplace"), TARGETS (default "cpu gpu"),
# FILTER (regular expression selecting the labels of the configurations),
# RUN_TIMEOUT (maximum duration of a run in seconds, 1800 by default).
#
# Outputs: one line per run in <results>/runs.txt (key=value, see the
# benchmark), the errors in <results>/errors.txt, a table in
# <results>/summary.md (cases/summary-solvers.py).
set -u
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
B=${1:-$P/build/ginkgo/mfem-mgis-trilinos}
out=${2:-$P/results/solvers-$(date +%Y-%m-%d-%H%M)}
exe=$B/benchmarks/LinearSolversBenchmark
size=${SIZE:-16}
problems=${PROBLEMS:-elasticity heterogeneous laplace}
targets=${TARGETS:-cpu gpu}
filter=${FILTER:-.}
run_timeout=${RUN_TIMEOUT:-1800}
if [ ! -x "$exe" ]; then
  echo "no benchmark $exe (configure mfem-mgis with -DBENCHMARKS=ON)" >&2
  exit 1
fi
mkdir -p "$out"
# environment of the build: MPI given to CMake, CUDA of scripts/env.sh
MPI=$(dirname "$(dirname "$(/usr/bin/grep -m1 '^MPI_CXX_COMPILER:' "$B/CMakeCache.txt" | sed 's/.*=//')")")
CUDA=${PROTO_CUDA:-$HOME/Bureau/spack-mmm-externals/cuda-13.0.2}
T=$P/install/h-interface/tfel
M=$P/install/h-interface/mgis
export PATH=$CUDA/bin:$MPI/bin:/usr/local/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=$T/lib:$M/lib:$CUDA/lib64
export OMPI_MCA_btl=self,vader
export OMP_NUM_THREADS=1
unset OMP_PROC_BIND OMP_PLACES

# configurations: label|placement|solver|parameters[|arguments]. The
# placement is `any` (CPU and GPU), `cpu` or `gpu`; the optional arguments
# are given to the benchmark (for instance `-bv`, unknowns ordered by vector
# dimension, which the elasticity options of BoomerAMG require). The
# tolerances are relative (1e-8), with at most 2000 iterations.
it='RelativeTolerance=1e-8,AbsoluteTolerance=0.,MaximumNumberOfIterations=2000'
configurations=(
  # Krylov solvers of MFEM and hypre, preconditioners of hypre
  "cg-jacobi|any|CGSolver|{$it,Preconditioner={Name=OperatorJacobiSmoother}}"
  "cg-boomeramg|any|CGSolver|{$it,Preconditioner={Name=HypreBoomerAMG}}"
  "cg-boomeramg-elasticity|any|CGSolver|{$it,Preconditioner={Name=HypreBoomerAMG,Options={Strategy=Elasticity}}}|-bv"
  "pcg-boomeramg|any|HyprePCG|{VerbosityLevel=0,Tolerance=1e-8,MaximumNumberOfIterations=2000,Preconditioner={Name=HypreBoomerAMG}}"
  "gmres-mgr|any|GMRESSolver|{$it,KDim=50,Preconditioner={Name=HypreMGR,Options={RelaxationType=2}}}"
  "cg-nodal-block-jacobi|any|CGSolver|{$it,Preconditioner={Name=NodalBlockJacobi}}"
  # AmgX (GPU)
  "cg-amgx|gpu|CGSolver|{$it,Preconditioner={Name=AmgX}}"
  "cg-amgx-aggregation|gpu|CGSolver|{$it,Preconditioner={Name=AmgX,Options={Multigrid=Aggregation}}}"
  "amgx-pcg|gpu|AmgXPCG|{$it}"
  # PETSc (host); ordered by vector dimension, GAMG receives the rigid body
  # modes
  "petsc-cg-gamg|cpu|PETScLinearSolver|{$it,Options={ksp_type=cg,pc_type=gamg}}|-bv"
  # the same, stopped on the true residual (the conjugate gradient of PETSc
  # stops on the preconditioned residual by default)
  "petsc-cg-gamg-true-residual|cpu|PETScLinearSolver|{$it,Options={ksp_type=cg,pc_type=gamg,ksp_norm_type=unpreconditioned}}|-bv"
  # Ginkgo
  "ginkgo-cg-jacobi|any|GinkgoCG|{$it,Preconditioner={Name=GinkgoJacobi}}"
  "ginkgo-cg-multigrid|any|GinkgoCG|{$it,Preconditioner={Name=GinkgoMultigrid}}"
  # Trilinos
  "belos-cg-relaxation|any|BelosSolver|{Solver=CG,$it,Preconditioner={Name=Ifpack2}}"
  # Chebyshev polynomial of degree 3 (degree 1, the default of Ifpack2, is
  # a damped Jacobi relaxation)
  "belos-cg-chebyshev|any|BelosSolver|{Solver=CG,$it,Preconditioner={Name=Ifpack2,Options={Type=CHEBYSHEV,Parameters={\"chebyshev: degree\"=3}}}}"
  "belos-cg-muelu|any|BelosSolver|{Solver=CG,$it,Preconditioner={Name=MueLu}}"
  "belos-cg-frosch|any|BelosSolver|{Solver=CG,$it,Preconditioner={Name=FROSch}}"
  "cg-muelu|any|CGSolver|{$it,Preconditioner={Name=MueLu}}"
  # direct solvers
  "mumps|cpu|MUMPSSolver|{Symmetric=true,PositiveDefinite=true}"
  "superlu|cpu|SuperLUSolver|{}"
  "strumpack|cpu|STRUMPACKSolver|{}"
  "mkl-cpardiso|cpu|MKLClusterSparseSolver|{}"
  "petsc-pastix|cpu|PETScLinearSolver|{Options={ksp_type=preonly,pc_type=lu,pc_factor_mat_solver_type=pastix}}"
  "amesos2-klu2|cpu|Amesos2Solver|{}"
  # Tacho (Cholesky factorization of Trilinos, on the GPU when Tpetra
  # computes there)
  "amesos2-tacho|any|Amesos2Solver|{Solver=Tacho}"
  "cudss|gpu|CuDSSSolver|{Symmetric=true,PositiveDefinite=true}"
)

# arguments of the benchmark for a problem
problem_arguments() {
  case "$1" in
    elasticity) echo "-p elasticity -n $size" ;;
    heterogeneous) echo "-p elasticity -n $size -c 1000 -i 4" ;;
    laplace) echo "-p laplace -n $((4 * size))" ;;
  esac
}

: > "$out/runs.txt"
: > "$out/errors.txt"
{
  echo "date: $(date)"
  echo "build: ${B#"$P"/}"
  echo "size: $size"
  echo "mfem-mgis: $(git -C "$P/scratch/ginkgo/mfem-mgis" describe --always --dirty 2>/dev/null)"
  nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null
  LC_ALL=C lscpu | /usr/bin/grep "Model name"
} > "$out/context.txt"

for problem in $problems; do
  for configuration in "${configurations[@]}"; do
    IFS='|' read -r label placement solver parameters arguments <<< "$configuration"
    echo "$label" | /usr/bin/grep -Eq "$filter" || continue
    for target in $targets; do
      if [ "$placement" != any ] && [ "$placement" != "$target" ]; then
        continue
      fi
      if [ "$target" = cpu ]; then
        launcher=(taskset -c 0-7 mpirun -np 8 --bind-to core
          -x LD_LIBRARY_PATH -x OMPI_MCA_btl -x OMP_NUM_THREADS)
        device=()
      else
        # one process, on the cores 0 and 1 (the other cores may be used by
        # timed computations of other work)
        launcher=(taskset -c 0,1 mpirun -np 1 --bind-to none
          -x LD_LIBRARY_PATH -x OMPI_MCA_btl -x OMP_NUM_THREADS)
        device=(--device cuda)
      fi
      log=$out/$problem-$label-$target.log
      # shellcheck disable=SC2046
      timeout "$run_timeout" "${launcher[@]}" "$exe" $(problem_arguments "$problem") \
        ${arguments:-} -s "$solver" -sp "$parameters" \
        -lb "$problem/$label/$target" "${device[@]}" > "$log" 2>&1
      status=$?
      if /usr/bin/grep -q "^label=" "$log"; then
        /usr/bin/grep "^label=" "$log" | sed "s/$/ target=$target status=$status/" >> "$out/runs.txt"
      else
        echo "$problem/$label/$target: exit $status: $(tail -c 300 "$log" | tr '\n' ' ')" >> "$out/errors.txt"
      fi
      tail -1 "$out/runs.txt"
    done
  done
done
python3 "$P/cases/summary-solvers.py" "$out/runs.txt" > "$out/summary.md"
