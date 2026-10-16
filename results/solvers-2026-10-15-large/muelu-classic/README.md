# MueLu with its factories not written with Kokkos, GPU (16/10/2026)

Configurations of MueLu of `cases/measure-solvers.sh` with mfem-mgis
d08224ba8f, which sets `use kokkos refactor` to false by default: MueLu
uses its factories not written with Kokkos (`TentativePFactory`,
`CoalesceDropFactory`) instead of those written with Kokkos, its default
on a GPU. One process on the GPU (RTX 5080), `SIZE=32`, whole machine
reserved, one repetition. Table: `summary.md`. Compare with
`../gpu-conversion` (factories written with Kokkos).

| problem | configuration | iterations (Kokkos factories) | setup (s) | resolution (s) | second operator (s) |
|---|---|---|---|---|---|
| elasticity | Jacobi smoother | 12 (22) | 0.51 (0.81) | 0.040 (0.070) | 0.58 (0.90) |
| elasticity | Jacobi, hierarchy reused | 12 (23) | 0.14 (0.13) | 0.040 (0.074) | 0.18 (0.20) |
| elasticity | Chebyshev of degree 2 | stopped after 7 (stopped after 6) | 0.51 | 0.025 | 0.56 |
| elasticity | symmetric Gauss-Seidel | 10 | 0.54 | 1.86 | 2.42 |
| heterogeneous | Jacobi, hierarchy reused | 25 (47) | 0.14 | 0.080 (0.148) | 0.22 (0.32) |
| heterogeneous | Chebyshev of degree 2 | 20 (40) | 0.54 | 0.067 | 0.60 |
| laplace | Jacobi, hierarchy reused | 10 (10) | 0.070 | 0.054 | 0.12 (0.13) |

- With the rigid body modes, the factories not written with Kokkos halve
  the iterations and shorten the setup on the GPU. On the Laplace
  equation (scalar unknowns), nothing changes.
- MueLu with a Jacobi smoother and its hierarchy reused: 0.18 seconds per
  operator on the elasticity, against 0.33 to 0.36 for BoomerAMG of hypre
  on the GPU, and 0.22 with the inclusions (0.51 to 0.57).
- The Chebyshev smoother of degree 2 still stops the conjugate gradient on
  the homogeneous beam (residual 2.5e-2 after 7 iterations: the
  preconditioner is not positive definite, see `../muelu-tuning`), and
  converges in 20 iterations with the inclusions.
