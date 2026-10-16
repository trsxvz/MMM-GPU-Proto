# Choosing a linear solver in mfem-mgis

Status: recommendations drawn on 16/10/2026 from the campaign of phase 5 of
[plan-solver-libraries.md](plan-solver-libraries.md), on one machine
(Ryzen 7 9800X3D, 8 cores, RTX 5080). They are a starting point for each
kind of problem, to be confirmed on the target machine.

## How the measurements were made

`cases/measure-solvers.sh` runs the benchmark `LinearSolversBenchmark` of
mfem-mgis on assembled systems: the linear elasticity of a cantilever beam,
homogeneous or with stiff cubic inclusions (contrast 1000), and the Laplace
equation on a cube. Each configuration runs on 8 processes on the CPU and,
when it computes on the GPU, on one process on the GPU. The relative
tolerance is 1e-8.

The time that matters for a nonlinear computation is the one of a second
operator of the same sparsity: setup plus resolution, paid at each Newton
iteration. The first setup also includes costs paid once.

Results: [first pass](../results/solvers-cpu-2026-10-15/README.md) (56,355
unknowns for the elasticity, 274,625 for the Laplace equation, all the
configurations), [second pass](../results/solvers-2026-10-15-byvdim/README.md)
(GAMG with the rigid body modes),
[large elasticity](../results/solvers-2026-10-15-large/elasticity/README.md)
(421,443 unknowns),
[large Laplace](../results/solvers-2026-10-15-large/laplace/README.md)
(2,146,689 unknowns), solvers of Trilinos at the
[size of the first passes](../results/solvers-2026-10-16-trilinos/README.md)
and at the [large sizes](../results/solvers-2026-10-15-large/muelu/README.md).
Trilinos built with CUDA computes on the GPU in all the runs (its runs
labelled `cpu` are 8 processes sharing the GPU).

## Recommendations by case

### Mechanics on the CPU

- Conjugate gradient of PETSc with GAMG, the unknowns ordered by vector
  dimension (`FiniteElementSpaceOrdering` set to `byVDIM`), so that GAMG
  receives the rigid body modes, and the option `ksp_norm_type
  unpreconditioned` so that the tolerance applies to the true residual:
  1.76 seconds per operator at 421,443 unknowns, twice as fast as
  BoomerAMG.
- Without PETSc, or with the unknowns ordered by nodes: conjugate gradient
  with BoomerAMG and its systems options (the default of mfem-mgis for
  vector unknowns).
- The elasticity options of BoomerAMG (`Strategy` set to `Elasticity`) do
  not reduce the iterations on the beam and cost three times more.

### Mechanics on the GPU

- Conjugate gradient with BoomerAMG of hypre on the GPU (systems options):
  0.33 second per operator at 421,443 unknowns, 5 times faster than the
  best configuration on 8 CPU cores, 6 times with the inclusions.
- The conjugate gradient with the Jacobi preconditioner is as fast on a
  homogeneous material, but its iterations grow quickly with the contrast
  (1797 against 109 with the inclusions).
- The classical multigrid of AmgX ignores the components of the unknowns:
  ten times slower than BoomerAMG on the elasticity (505 iterations at
  421,443 unknowns). Its aggregation multigrid converges on the elasticity
  only with a symmetric cycle and a Jacobi smoother damped by 0.6, the
  default of mfem-mgis since 16/10/2026 (141 iterations at 56,355
  unknowns, against 54 for BoomerAMG); the configuration
  `PCG_AGGREGATION_JACOBI` of AmgX diverges. Giving the matrix to AmgX by
  nodal blocks (`NodalBlocks`) does not reduce the iterations. AmgX suits
  scalar problems.
- MueLu of Trilinos with the rigid body modes needs fewer iterations than
  BoomerAMG (20 to 23 against 37 to 59 at 421,443 unknowns), and its cycles
  are the cheapest on the GPU with a Jacobi smoother (a resolution in
  0.07 second instead of 3.5 with its default symmetric Gauss-Seidel), but
  its setup takes 1.2 seconds: worth it when the hierarchy is reused for
  many resolutions (modified Newton methods), not when it is built again
  at each iteration. The Chebyshev polynomial of Ifpack2 (degree 3) is
  close to BoomerAMG on small systems and 2.5 times slower at 421,443
  unknowns.

### Scalar problems (heat transfer)

- Conjugate gradient with BoomerAMG, on the CPU and on the GPU (13 to 14
  iterations at 274,625 unknowns). On the GPU, AmgX converges in 11 to 13
  iterations but costs 2.5 to 3 times more per operator.
- At 2,146,689 unknowns on the GPU, the conjugate gradient with the Jacobi
  preconditioner is the fastest (0.18 second per operator, 183
  iterations), before BoomerAMG (0.42 second): the setup of a multigrid
  costs more than the iterations it saves on this well conditioned
  problem. On 8 cores, BoomerAMG (1.93 second) is 10 times slower than
  the best on the GPU.

### Direct solvers and repeated resolutions

- On the GPU, cuDSS factorizes 3.4 times faster than MUMPS on 8 cores
  (5.9 against 19.8 seconds at 421,443 unknowns) and solves 4 times
  faster. It is the choice for the small systems on the GPU (fewer than
  about 10,000 unknowns) and when a factorization is reused for many
  resolutions (modified Newton methods).
- On the CPU, MUMPS and MKL Cluster Sparse Solver factorize in similar
  times (20 to 22 seconds at 421,443 unknowns), and MUMPS solves fastest
  (0.2 second against 0.6).
- Tacho of Trilinos (on the GPU) factorizes the beam of 56,355 unknowns in
  0.72 second, 2.4 times slower than cuDSS; KLU2 (sequential) is not
  suited to 3D finite element matrices (38 seconds).

### To avoid

- The multigrid of Ginkgo on the elasticity (aggregation without the rigid
  body modes: no convergence).
- MGR of hypre with an AMG F-relaxation on the elasticity: its iterations
  grow with the size (362 at 56,355 unknowns, 1147 with the inclusions).
- The default stop of the conjugate gradient of PETSc when comparing with
  other solvers: it stops on the preconditioned residual (true residual
  40 to 4000 times larger than the tolerance in these measurements).
- FROSch with its default subdomain solver (KLU2) on 3D problems: minutes
  per setup; a solver of Amesos2 suited to finite element matrices must be
  given to it. FROSch does not support processes owning no unknown.
- The default smoother of MueLu on a GPU (symmetric Gauss-Seidel).

## To be completed

- The solvers of Trilinos on the CPU (a Trilinos built without CUDA).
- FROSch with a subdomain solver suited to 3D problems (Tacho).
- The setup of MueLu on the GPU, which dominates its cost.
