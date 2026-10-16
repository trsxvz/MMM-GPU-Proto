# Linear solvers, large Laplace systems (16/10/2026)

Pass of `cases/measure-solvers.sh` at a larger size (`SIZE=32`): the
Laplace equation on a cube of 128³ elements (2,146,689 unknowns), with the
iterative configurations only (the direct solvers would factorize a 3D
system of two million unknowns). 8 processes on the CPU, one process on
the GPU (RTX 5080). Build `trilinos` of mfem-mgis (environment
`mmm-trilinos`, all the libraries), commits fc9ba3b494 and eb3fecda10,
quiet machine, one repetition. Table: `summary.md`. The runs of MueLu on
one process ran out of memory on the GPU here, before the algorithm of the
products of sparse matrices was changed: they are in `../muelu`. Those
labelled `cpu` (4.19 and 5.02 seconds in `summary.md`) are 8 processes
sharing the GPU, Trilinos built with CUDA computing on the GPU.

## Main results

| configuration | CPU, 8 processes (s) | GPU (s) | iterations CPU / GPU |
|---|---|---|---|
| CG + BoomerAMG | 1.93 | 0.42 | 14 / 16 |
| CG + GAMG (PETSc), true residual | 3.03 | (host only) | 11 |
| CG + Jacobi | 3.49 | 0.18 | 183 / 183 |
| CG + aggregation multigrid of AmgX | | 0.55 | 47 |
| CG + classical multigrid of AmgX | | 0.99 | 12 |
| CG + multigrid of Ginkgo | 6.37 | 1.96 | 52 / 46 |
| CG + MueLu (`../muelu` on the GPU) | 4.19 (8 processes on the GPU) | 2.34 | 9 / 8 |

Times for a second operator of the same sparsity (setup plus resolution).

- On the GPU, the conjugate gradient with the Jacobi preconditioner is the
  fastest at this size (0.18 seconds, 183 iterations): the setup of a
  multigrid costs more than the iterations it saves on this problem.
  BoomerAMG on the GPU follows (0.42 seconds), 4.6 times faster than on 8
  cores.
- On the CPU, BoomerAMG is the fastest (1.93 seconds), ahead of GAMG of
  PETSc (3.03 seconds stopped on the true residual).
- From the best on 8 cores (1.93 seconds) to the best on the GPU (0.18
  seconds): a factor of 10.
- The classical multigrid of AmgX needs few iterations (12) but its setup
  costs 0.9 seconds; its aggregation multigrid (symmetric cycle of
  mfem-mgis) sets up faster and needs 47 iterations.
