# MueLu with a Jacobi smoother, large systems (16/10/2026)

Configuration `cg-muelu-jacobi` of `cases/measure-solvers.sh`: the
conjugate gradient of MFEM preconditioned by MueLu, smoothed by two sweeps
of a Jacobi relaxation damped by 0.6 (`smoother: type` set to
`RELAXATION`), against `cg-muelu`, MueLu with its default smoother
(symmetric Gauss-Seidel). One process on the GPU (RTX 5080), `SIZE=32`:
the beam of 421,443 unknowns, homogeneous and with stiff inclusions, and
the Laplace equation of 2,146,689 unknowns. Build `trilinos` of mfem-mgis
f904a70d36. Two repetitions (the logs of the second one end with
`-rep2`), with the GPU reserved for these runs (compilations ran on the
other CPU cores); `summary.md` gives the minimum times. A first run, during
which another computation used the GPU, was discarded.

| problem | smoother | iterations | setup (s) | resolution (s) | second operator (s) |
|---|---|---|---|---|---|
| elasticity | symmetric Gauss-Seidel | 20 | 1.25 | 3.84 | 5.25 |
| elasticity | Jacobi, 2 sweeps | 22 to 23 | 1.25 | 0.071 | 1.32 |
| heterogeneous | symmetric Gauss-Seidel | 39 | 1.25 | 7.56 | 8.85 |
| heterogeneous | Jacobi, 2 sweeps | 47 | 1.25 | 0.148 | 1.40 |
| laplace | symmetric Gauss-Seidel | 8 | 0.53 | 1.67 | 2.40 |
| laplace | Jacobi, 2 sweeps | 10 | 0.53 | 0.054 | 0.66 |

- With the Jacobi smoother, a resolution is 54 times faster on the
  elasticity (0.071 seconds against 3.84) for two or three more
  iterations, and 31 times faster on the Laplace equation. The setup of
  MueLu (1.25 seconds on the elasticity) then makes most of the time per
  operator.
- The resolution alone is faster than BoomerAMG of hypre on the GPU (0.16
  to 0.19 seconds on the elasticity, 59 iterations, in `../elasticity`),
  but the time per operator is four times larger (1.32 seconds against
  0.33 to 0.36): MueLu pays off when its hierarchy is reused for several
  resolutions.
- On the Laplace equation, the resolution takes 0.054 seconds, against
  0.18 seconds per operator for the conjugate gradient with the Jacobi
  preconditioner, the fastest configuration at this size (`../laplace`).
- The number of iterations with the Jacobi smoother on the elasticity
  differs between the repetitions (22 and 23), although the aggregation is
  deterministic, probably because the computations on the GPU are not
  reproducible to the last bit (not checked).
