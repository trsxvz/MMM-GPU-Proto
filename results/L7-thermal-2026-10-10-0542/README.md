# Steady-state thermal by kernels, first measurement (10/10/2026)

Nonlinear steady-state thermal problem in the beam (`cases/thermal`, r5,
279,873 unknowns, 5 steps), conductivity depending on the temperature,
temperature imposed at both ends. mfem-mgis 6cae30af99. Script:
`cases/measure-L7-thermal.sh`. One RTX 5080, 16 cores, load 1.9.

| Case | Time | Newton iterations |
|---|---|---|
| CPU 8 processes, point by point | 9.17 s | 18 |
| CPU 8 processes, kernels | 4.89 s | 18 |
| GPU, kernels | 3.01 s | 18 |
| GPU, kernels, Eisenstat-Walker 1 | 4.03 s | 25 |
| GPU, kernels, matrix free (GMRES + Jacobi) | 183.40 s | 18 |

On the GPU, the first integration cost 1.76 s: one-off computation of the
geometry of the elements (`computeElementsGeometry`), done on the host by
a finite element space of order one and a projection of the coordinates
(1.0 s and 0.83 s, measured separately). This cost is removed by the
computation of the jacobians by a kernel: see
`results/L7-thermal-2026-10-10-0606/`.
