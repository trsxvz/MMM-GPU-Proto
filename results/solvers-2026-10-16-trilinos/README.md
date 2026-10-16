# Solvers of Trilinos, size of the first passes (16/10/2026)

Pass of `cases/measure-solvers.sh` restricted to the configurations of
Trilinos (Belos, Ifpack2, MueLu, FROSch, Amesos2 with KLU2 and Tacho), on
the problems of the first passes: the cantilever beam of 56,355 unknowns,
homogeneous and with stiff inclusions, and the Laplace equation of 274,625
unknowns. Build `trilinos` of mfem-mgis (environment `mmm-trilinos`,
Trilinos 17.2.1 with Kokkos 5.2.1 and CUDA), commits fc9ba3b494 and
eb3fecda10. One repetition, quiet machine. Table: `summary.md`
(`runs.txt` and `laplace-rest/runs.txt`).

**Trilinos built with CUDA computes on the GPU in all the runs**: Tpetra
places its matrices and vectors on the GPU whatever the device of MFEM.
The runs labelled `cpu` are 8 processes sharing the GPU, not computations
on the CPU. A comparison with the CPU needs a Trilinos built without CUDA.

Not measured: FROSch on the Laplace equation (its default subdomain
solver, KLU2, took more than 4 minutes on 8 subdomains of the cube, and
FROSch with one subdomain on the GPU is a factorization by KLU2 on the
host: 49 seconds on the beam), KLU2 on the Laplace equation (sequential
factorization of a 3D problem).

## Main results (GPU, one process)

| problem | configuration | iterations | second operator (s) |
|---|---|---|---|
| elasticity | CG + Chebyshev (Ifpack2, degree 3) | 114 | 0.081 |
| elasticity | CG + Jacobi relaxation (Ifpack2) | 319 | 0.084 |
| elasticity | CG + MueLu | 22 to 25 | 0.63 to 0.79 |
| elasticity | Tacho (Amesos2) | direct | 0.72 (solve 0.003) |
| heterogeneous | CG + Chebyshev (degree 3) | 299 | 0.126 |
| laplace | CG + Chebyshev (degree 3) | 32 | 0.071 |
| laplace | CG + MueLu | 7 to 8 | 0.30 to 0.39 |
| laplace | Tacho | direct | 5.42 |

- The Chebyshev polynomial of Ifpack2 (degree 3) is close to BoomerAMG of
  hypre on the GPU on these sizes (0.070 second per operator on the
  elasticity, 54 iterations), and needs 2.8 times fewer iterations than
  the Jacobi relaxation with the inclusions.
- MueLu with its default options (smoothed aggregation with the rigid
  body modes, symmetric Gauss-Seidel smoother, coarse problem solved on
  the host) reduces the iterations as much as BoomerAMG, but each cycle
  is about 15 times more expensive on the GPU: its options for the GPU are
  still to be studied.
- Tacho factorizes the elasticity in 0.72 second on the GPU (second
  operator), close to MUMPS on 8 cores (0.59) and 2.4 times slower than
  cuDSS (0.30); on the Laplace equation, 5.4 seconds against 3.2 with
  cuDSS.
