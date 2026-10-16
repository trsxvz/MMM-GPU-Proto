# Solvers of Trilinos, large elasticity systems, first run (16/10/2026)

First run of the configurations of Trilinos on the beam of 421,443
unknowns (build `trilinos`, commits fc9ba3b494 and eb3fecda10): the
Chebyshev polynomial of Ifpack2 ran (`runs.txt`), MueLu ran out of memory
on the GPU in the products of sparse matrices of cuSPARSE
(`errors.txt`). The runs of MueLu, with the products of Kokkos Kernels,
and the tables are in `../muelu`.
