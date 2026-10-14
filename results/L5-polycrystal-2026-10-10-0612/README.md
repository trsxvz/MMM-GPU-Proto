# Polycrystal, check after the jacobians by a kernel (10/10/2026)

Same cases as `results/L5-polycrystal-2026-10-10-0502/`, with the
computation of the jacobians by a kernel (mfem-mgis 8b0b47acfb;
`summary.txt` gives the base commit). Load 4.8 at the start (end of the
measurements of the beam).

| Case | CPU, 8 processes | GPU |
|---|---|---|
| Newton | 29.65 s | 11.86 s |
| Newton, EW1 | 34.08 s | 11.60 s |
| NewtonK0 | 25.12 s | 11.25 s |
| NewtonK0, EW1 | 23.50 s | 9.58 s |

Total times (`cases/summary-profiling.py`). Everything is equal to the
previous measurement within 0.5 %: the mesh is small, the geometry did not
weigh there. Macroscopic stresses equal to Newton up to 2.4e-13 without
Eisenstat-Walker, from 9.6e-7 to 3.5e-6 with it (within the tolerance of
the microscopic Newton, see the previous measurement).
