# Inexact Newton on the beam, jacobians by a kernel (10/10/2026)

Same cases as `results/L5-beam-2026-10-10-0429/` (beam r4, 111,843
unknowns, 10 steps, all the kernels, script `cases/measure-L5-ew.sh`,
initial tolerance of Eisenstat-Walker 0.5), with the computation of the
jacobians by a kernel (mfem-mgis 8b0b47acfb; `summary.txt` gives the base
commit). Load 1.1 at the start.

| Case | GPU, before | GPU | CPU 8 processes, before | CPU 8 processes |
|---|---|---|---|---|
| Newton | 3.12 s | 2.98 s | 13.59 s | 13.57 s |
| Newton, Eisenstat-Walker 1 | 3.01 s | 2.74 s | 15.73 s | 15.76 s |
| Newton, Eisenstat-Walker 2 | 3.32 s | 3.10 s | 17.70 s | 17.82 s |
| NewtonK0 | 2.57 s | 2.56 s | 11.56 s | 11.59 s |
| NewtonK0, Eisenstat-Walker 1 | 2.14 s | **1.90 s** | 10.93 s | **11.00 s** |
| NewtonK0, Eisenstat-Walker 2 | 2.19 s | 1.99 s | 12.30 s | 12.26 s |

On the GPU, the gain goes from 0 to 0.27 s depending on the case, one run
per case: the geometry of 197,000 tetrahedra probably cost about 0.2 s on
the host (one eighth of the 1.8 s measured in r5). With 8 processes, each
process only computes one eighth of it: no measurable difference.
