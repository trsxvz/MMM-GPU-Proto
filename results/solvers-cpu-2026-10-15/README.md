# Campaign of the linear solvers, without Trilinos (15/10/2026)

Phase 5 of [docs/plan-solver-libraries.md](../../docs/plan-solver-libraries.md),
first pass. Script: `cases/measure-solvers.sh` (benchmark
`LinearSolversBenchmark` of mfem-mgis, patch 0231), size 16. Build
`build/ginkgo/mfem-mgis-solvers` (environment `mmm-solvers`, mfem-mgis
8d829b4c26 with the module of Trilinos in progress, not enabled). The
solvers of Trilinos are not in this pass: the MFEM of `mmm-trilinos` is
not built yet. This directory: 8 processes on the cores 0 to 7, quiet
machine. GPU part: `../solvers-gpu-2026-10-15` (one process on the cores
0 and 1, RTX 5080 free). Tables per target: `summary.md` of each
directory.

Systems (assembled by MFEM, order 1, hexahedra):

- elasticity: cantilever beam of 64 x 16 x 16 elements, 56,355 unknowns,
  ordered by nodes;
- heterogeneous: the same beam, the coefficients of a checkerboard of cubic
  inclusions (4 across the section) multiplied by 1000;
- laplace: cube of 64^3 elements, 274,625 unknowns.

Relative tolerance 1e-8 (the residuals of the table are the true relative
residuals; GMRES of MFEM and the solvers preconditioned on the left stop on
another norm), at most 2000 iterations.

## CPU (8 processes) against GPU (1 process)

Setup plus first resolution, in seconds (the preconditioners of hypre are
built at their first application), then the following resolutions with the
same operator.

| System | Configuration | CPU, setup + first | GPU, setup + first | CPU, solve | GPU, solve |
|---|---|---|---|---|---|
| elasticity | CG + Jacobi (MFEM) | 0.094 | 0.037 | 0.091 | 0.035 |
| elasticity | CG + BoomerAMG (MFEM, hypre) | 0.203 | 0.111 | 0.119 | 0.042 |
| elasticity | PCG + BoomerAMG (hypre) | 0.205 | 0.076 | 0.136 | 0.039 |
| elasticity | direct: MUMPS (CPU), cuDSS (GPU) | 1.22 | 0.405 | 0.021 | 0.004 |
| heterogeneous | CG + Jacobi | 0.279 | 0.084 | 0.286 | 0.083 |
| heterogeneous | PCG + BoomerAMG | 0.318 | 0.107 | 0.207 | 0.071 |
| heterogeneous | direct: MUMPS (CPU), cuDSS (GPU) | 1.18 | 0.353 | 0.020 | 0.004 |
| laplace | CG + BoomerAMG | 0.120 | 0.051 | 0.046 | 0.015 |
| laplace | direct: MUMPS (CPU), cuDSS (GPU) | 5.66 | 3.24 | 0.098 | 0.016 |

- On these sizes, the GPU is 2 to 3 times faster than the 8 processes for
  the iterative solvers, and cuDSS factorises 2 to 3 times faster than
  MUMPS and solves 5 times faster.
- BoomerAMG divides the iterations by 6 on the elasticity (54 against 319
  on the GPU, 30 against 319 on the CPU: BoomerAMG is not the same on both)
  but, on a system this small, Jacobi stays competitive: one cycle of
  BoomerAMG costs about as much as 6 iterations of Jacobi.
- Fastest direct solver on the CPU: MKL Cluster Sparse Solver, then MUMPS,
  STRUMPACK, PaStiX (through PETSc) and SuperLU_DIST; MUMPS solves fastest
  (0.02 s).
- Do not converge or are much slower: the multigrid of Ginkgo on the
  elasticity (aggregation without the rigid body modes, 2000 iterations,
  known), MGR of hypre with an F-relaxation by AMG (GMRES, 362 iterations
  on the elasticity, 1147 with inclusions, against about 35 on the small
  beam of `HypreMGRTest`: its iterations grow with the size), GAMG of PETSc
  on the elasticity (177 to 284 iterations, probably because the rigid
  body modes are not given to it, not checked).
- `cg-boomeramg-elasticity` did not use the elasticity options: they
  require the unknowns ordered by vector dimension, and the benchmark
  ordered them by nodes (mfem-mgis then falls back to the systems options,
  with a warning; on the GPU, hypre does not support them anyway). The
  configuration now passes `-bv`; to be measured again on a quiet machine.
- Laplace: MGR is refused (one component), as expected.
