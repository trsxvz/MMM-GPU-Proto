# Trilinos: conversion of the matrices on the GPU (16/10/2026)

Since the commit 3a3909cec8 of mfem-mgis, the matrices of hypre are
converted to Tpetra by kernels of Kokkos, on the GPU when Tpetra and hypre
compute there. Before, they were read on the host, converted there, then
copied on the GPU. One process on the GPU (RTX 5080), `SIZE=32`, GPU
reserved for these runs, two repetitions (the logs of the second one end
with `-rep2`); `summary.md` gives the minimum times. The directory
`before` holds the runs of `belos-cg-relaxation` with the previous
conversion (build `trilinos`, f904a70d36).

The relaxation of Ifpack2 (one Jacobi sweep) has an almost free setup: its
setup is the conversion.

| problem | setup of `belos-cg-relaxation`, before (s) | after (s) |
|---|---|---|
| elasticity, 421,443 unknowns | 0.486 | 0.064 |
| heterogeneous, 421,443 unknowns | 0.484 | 0.064 |
| laplace, 2,146,689 unknowns | 0.450 | 0.029 |

MueLu with the Jacobi smoother of `../muelu-jacobi`, with and without the
reuse of the hierarchy (`reuse: type` set to `RP`, `../muelu-reuse`):

| problem | configuration | setup of the second operator (s) | second operator (s) | iterations |
|---|---|---|---|---|
| elasticity | hierarchy built again | 0.82 (1.28 before) | 0.90 | 22 |
| elasticity | hierarchy reused | 0.17 (0.60 before) | 0.24 | 23 |
| heterogeneous | hierarchy reused | 0.17 | 0.32 | 47 |
| laplace | hierarchy built again | 0.24 | 0.29 | 10 |
| laplace | hierarchy reused | 0.074 | 0.13 | 10 |

- The conversion is 7.6 times faster on the elasticity and 15 times on the
  Laplace equation.
- With the reuse of its hierarchy, MueLu is now the fastest preconditioner
  measured on the GPU per operator: 0.24 seconds on the elasticity
  against 0.33 to 0.36 for BoomerAMG of hypre (`../elasticity`), 0.32
  against 0.51 to 0.57 with the stiff inclusions, and 0.13 seconds on the
  Laplace equation against 0.18 for the conjugate gradient with the
  Jacobi preconditioner (`../laplace`).
- Variants of MueLu (`cg-muelu-chebyshev`: Chebyshev smoother of degree 2;
  `cg-muelu-jacobi-unsmoothed`: tentative prolongators): the unsmoothed
  aggregation saves little setup (0.74 seconds against 0.82) and triples
  the iterations. With the Chebyshev smoother, the conjugate gradient of
  MFEM stops after 6 iterations on the homogeneous elasticity (not
  converged: the preconditioner is not positive definite, probably
  because of an underestimated largest eigenvalue, not checked), and
  converges in 40 iterations with the inclusions.
