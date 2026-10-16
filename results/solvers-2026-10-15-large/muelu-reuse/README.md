# MueLu: reuse of the hierarchy, large systems (16/10/2026)

Configuration `cg-muelu-jacobi-reuse` of `cases/measure-solvers.sh`: the
conjugate gradient of MFEM preconditioned by MueLu with the Jacobi smoother
of `../muelu-jacobi`, the prolongators and the restrictions of the first
operator being kept for the second one (`reuse: type` set to `RP`, the
coarse operators and the smoothers being computed again), against
`cg-muelu-jacobi`, whose hierarchy is built again. The second operator of
the benchmark is the first one times two. One process on the GPU (RTX
5080), `SIZE=32`. Build `trilinos-host` of mfem-mgis (`context.txt`), GPU
reserved for these runs, one repetition. Table: `summary.md`.

| problem | hierarchy | setup of the second operator (s) | second operator (s) | iterations |
|---|---|---|---|---|
| elasticity | built again | 1.28 | 1.36 | 22 |
| elasticity | reused (`RP`) | 0.60 | 0.67 | 23 |
| heterogeneous | built again | 1.26 | 1.41 | 47 |
| heterogeneous | reused (`RP`) | 0.59 | 0.74 | 47 |
| laplace | built again | 0.62 | 0.68 | 10 |
| laplace | reused (`RP`) | 0.46 | 0.51 | 10 |

- Keeping the prolongators and the restrictions halves the setup on the
  elasticity (0.60 seconds instead of 1.28), for the same iterations.
- A run by hand keeping the whole hierarchy (`reuse: type` set to `full`)
  set up the second operator of the elasticity in 0.49 seconds: this is an
  upper bound of the conversion of the matrix from hypre to Tpetra, made on
  the host, and of the update of the finest level. The resolution then
  failed (residual 0.015 after 6.2 seconds), as the coarse operators of
  the first operator were kept: `full` only suits an unchanged operator.
- With the reuse, a Newton iteration on the elasticity costs 0.67 seconds
  with MueLu, against 0.33 to 0.36 seconds with BoomerAMG of hypre on the
  GPU (`../elasticity`). The conversion of the matrix on the host is the
  next cost to remove: done since, see `../gpu-conversion` (0.24 seconds
  per operator with the reuse).
