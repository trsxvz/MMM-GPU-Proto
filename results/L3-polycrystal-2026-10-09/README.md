# Measurements of work package L3 on the polycrystal (09/10/2026)

Polycrystal of mm-opera-hpc (5 grains, single crystal behaviour of UO2,
finite strain, orthotropic, periodic problem), 10 steps, on CPU, `generic`
library. Gradients, internal forces and element matrices of the tangent
computed element by element or by kernels (mfem-mgis ed91564643). Script:
`cases/measure-L3.sh polycrystal`.

| Case | Total | Integration | Residual | Tangent | Solve |
|---|---|---|---|---|---|
| 1 process, point by point | 149.73 | 104.52 | 4.16 | 8.53 | 32.20 |
| 1 process, batches | 148.03 | 102.74 | 4.16 | 8.34 | 32.47 |
| 1 process, gradients by kernels | 144.97 | 99.52 | 4.19 | 8.43 | 32.51 |
| 1 process, residual by kernels | 145.46 | 102.55 | 1.49 | 8.61 | 32.49 |
| 1 process, tangent by kernels | 151.04 | 102.70 | 4.14 | 11.50 | 32.38 |
| 1 process, everything by kernels | 145.37 | 99.55 | 1.49 | 11.53 | 32.47 |
| 8 processes, point by point | 30.35 | 14.45 | 0.59 | 1.69 | 13.38 |
| 8 processes, batches | 30.21 | 14.31 | 0.58 | 1.68 | 13.40 |
| 8 processes, everything by kernels | 29.89 | 13.85 | 0.21 | 2.16 | 13.43 |

Macroscopic stresses: identical to the point by point computation on the
same number of processes up to 5e-9 (relative to the maximum of each
column), at the level of the tolerance of the iterative solver.

Findings: the residual by kernels is 2.8 times faster. The gradients by
kernels save 3 s out of 103 (the behaviour dominates: see L2 for its
threads of execution). The element matrices by kernels stay slower on CPU
(+38 %), as on the beam.
