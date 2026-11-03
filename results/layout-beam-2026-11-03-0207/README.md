# Layout of the data at the integration points on the beam (03/11/2026)

Beam r4 of the prototype (111,843 unknowns, J2 plasticity in small strain,
tangent operator blocks of 36 components), 10 steps, Newton. mfem-mgis
b283483f4b, script `cases/measure-layout.sh beam`, executable built
against `install/mfem-mgis-layout`. One RTX 5080, 16 cores, quiet machine
(load 0.7). Stiffness matrices assembled in a fixed order.

The arrays of the material are stored point by point (AoS) or component by
component (SoA). The two layouts give the same norm of the solution
(48.83704586) and the same numbers of iterations.

| Computation | AoS | SoA |
|---|---|---|
| GPU, tangent operator assembled (PCG of hypre and BoomerAMG) | 5.51 s | 5.51 s |
| GPU, tangent operator applied without assembly (CG and Jacobi) | 17.06 s | 16.17 s |
| CPU, 8 processes (PCG of hypre and BoomerAMG) | 11.72 s | 14.04 s |

Total times (`cases/summary-profiling.py`). Detail of the phases which
depend on the layout:

| Computation | Phase | AoS | SoA |
|---|---|---|---|
| GPU, without assembly | resolution (action of the tangent operator at each iteration) | 15.45 s | 14.56 s |
| GPU, both | integration of the behaviour | 0.05 s | 0.04 s |
| CPU | integration of the behaviour | 0.18 s | 0.56 s |
| CPU | stiffness matrices | 0.42 s | 1.66 s |

On the GPU, with the current kernels of the finite elements (one thread per
element), the component by component layout brings 6 % on the action of
the tangent operator, nothing elsewhere. On the CPU, it makes the
integration of the behaviour 3 times slower and the stiffness matrices 4
times slower. Both agree with the micro benchmarks of
[results/layout-microbenchmarks-2026-11-02](../layout-microbenchmarks-2026-11-02/README.md).
