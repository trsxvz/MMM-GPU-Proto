# Polycrystal, check after the rework of the geometry (10/10/2026)

Same cases as `results/L5-polycrystal-2026-10-10-0502/`, with mfem-mgis
3e3cc5b0f2: jacobians computed by a kernel from the nodes or the vertices
of the elements, stiffness matrices computed per element on the host
(62908a890d). Load 2.6 at the start.

| Case | CPU, 8 processes | GPU |
|---|---|---|
| Newton | 29.35 s | 11.85 s |
| Newton, EW1 | 33.46 s | 11.57 s |
| NewtonK0 | 24.69 s | 11.23 s |
| NewtonK0, EW1 | **22.83 s** | **9.56 s** |

Total times (`cases/summary-profiling.py`).

- Macroscopic stresses equal to Newton up to 1.1e-14 without
  Eisenstat-Walker, from 9.6e-7 to 3.5e-6 with it: as before.
- On CPU, the assembly of the tangent goes from 2.01 s to 1.60 s (Newton,
  98 tangents) thanks to the computation of the matrices per element on
  the host; NewtonK0 with Eisenstat-Walker 1 from 23.50 s to 22.83 s.
- On the GPU, nothing changes (the kernel of the matrices keeps one thread
  per pair of nodes there).
