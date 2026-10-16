# Linear solvers, large elasticity systems (15/10/2026)

Pass of `cases/measure-solvers.sh` at a larger size (`SIZE=32`): the
cantilever beam of 128 x 32 x 32 elements (421,443 unknowns), homogeneous
(`elasticity`) and with stiff cubic inclusions of contrast 1000
(`heterogeneous`). Selection of configurations: most of the iterative ones which
converged in the first passes (not MGR, nor the elasticity options of
BoomerAMG, three times slower), and three of the fastest direct solvers
(MUMPS, MKL Cluster Sparse Solver, cuDSS). 8
processes on the CPU (Ryzen 7 9800X3D), one process on the GPU (RTX 5080).
Build `solvers` of mfem-mgis c11ccc241a (patches up to 0238; the commit
fbd1a6619d of `context.txt` only changes the documentation), quiet
machine, one repetition. Table: `summary.md`. The Laplace part of this pass
is in `../laplace` when measured.

## Main results

Time for a second operator of the same sparsity (setup plus resolution),
what each Newton iteration pays, in seconds:

| problem | best on the CPU (8 processes) | best on the GPU | GPU gain |
|---|---|---|---|
| elasticity | GAMG, true residual: 1.76 (26 it.) | BoomerAMG: 0.33 (59 it.) | 5.3 |
| heterogeneous | GAMG, true residual: 3.03 (59 it.) | BoomerAMG: 0.51 (109 it.) | 5.9 |
| elasticity, direct | MUMPS: 19.8 (solve 0.20) | cuDSS: 5.9 (solve 0.049) | 3.4 (4 for the solve) |

- On the CPU, GAMG with the rigid body modes (unknowns ordered by vector
  dimension) is twice as fast as BoomerAMG with its systems options (3.8
  seconds per operator, 37 iterations), stopped on the true residual. Its
  default stop, on the preconditioned residual, leaves a true residual of
  1e-5 (19 iterations, 1.57 seconds).
- On the GPU, BoomerAMG of hypre (systems options) and the conjugate
  gradient with the Jacobi preconditioner cost the same on the homogeneous
  beam (0.33 second), but the Jacobi preconditioner needs 1797 iterations
  with the inclusions (0.92 second against 0.51).
- The classical multigrid of AmgX (default configuration) needs 505 to
  956 iterations on the elasticity: it ignores the components of the
  unknowns. It is ten times slower than BoomerAMG on the GPU here.
- Ginkgo with the Jacobi preconditioner does not converge within 2000
  iterations with the inclusions (the conjugate gradient of MFEM with the
  same preconditioner needs 1797).
- The direct solvers factorize in about 20 seconds on the CPU (MUMPS,
  MKL Cluster Sparse Solver) and 5.9 seconds with cuDSS; they solve in
  0.05 to 0.6 seconds. They pay off when the factorization is reused for
  many resolutions (modified Newton methods, several right hand sides).
- From 56,355 to 421,443 unknowns (7.5 times more), the time per operator
  of BoomerAMG on the CPU grows 19 times and that of GAMG 15 times,
  probably because the matrix of the small system (50 to 70 MB) fits in
  the 96 MB of cache of the processor (not checked). On the GPU,
  the time of BoomerAMG grows 5 times.
