# Solvers of Trilinos, large systems (16/10/2026)

Configurations of Trilinos at the larger size (`SIZE=32`): the beam of
421,443 unknowns, homogeneous and with stiff inclusions, and the Laplace
equation of 2,146,689 unknowns. Build `trilinos` of mfem-mgis. Table:
`summary.md`, which merges:

- `runs.txt` and `repeat-elasticity/runs.txt`: MueLu, with the products of
  sparse matrices of Kokkos Kernels (`SPGEMM_KK`, default of mfem-mgis
  since this pass);
- `../elasticity-trilinos/runs.txt`: CG preconditioned by the Chebyshev
  polynomial of Ifpack2 (degree 3), built from the commits fc9ba3b494 and
  eb3fecda10. In this first run, MueLu ran out of memory on the GPU (16
  GB) in the products of sparse matrices of cuSPARSE, the default of
  Tpetra on CUDA (`../elasticity-trilinos/errors.txt`).

Trilinos built with CUDA computes on the GPU in all the runs: the runs
labelled `cpu` are 8 processes sharing the GPU.

## Main results (GPU, one process)

| problem | configuration | iterations | second operator (s) |
|---|---|---|---|
| elasticity | CG + Chebyshev (Ifpack2, degree 3) | 227 | 0.84 |
| elasticity | CG of MFEM + MueLu | 20 | 4.87 |
| elasticity | CG of Belos + MueLu | 23 | 6.46 |
| heterogeneous | CG + Chebyshev | 661 | 1.42 |
| heterogeneous | CG of MFEM + MueLu | 38 | 8.2 |
| laplace | CG of MFEM + MueLu | 8 | 2.34 |

- MueLu needs few iterations (20 on the elasticity, 8 on the Laplace
  equation) but its cycles are slow on the GPU with its default options:
  0.18 seconds per iteration on the elasticity, against 0.003 for
  BoomerAMG of hypre (0.33 seconds per operator, 59 iterations, in
  `../elasticity`). Its setup takes 1.2 seconds.
- The Chebyshev polynomial of degree 3 needs 227 iterations on the
  elasticity, 0.84 seconds per operator: 2.5 times slower than BoomerAMG on
  the GPU, and its iterations grow with the contrast of the inclusions.
- The products of sparse matrices of cuSPARSE (default of Tpetra on CUDA)
  need too much memory for the multigrid of MueLu at these sizes: those of
  Kokkos Kernels are used instead.
