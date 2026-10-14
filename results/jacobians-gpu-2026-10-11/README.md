# Inverses of the Jacobians: end to end on the GPU (11/10/2026)

Script: `cases/measure-jacobians-gpu.sh`. Beam r4 (10 steps) and thermal
r5 (5 steps), one process and the RTX 5080, everything by kernels,
`generic-parallel` behaviour compiled for CUDA. Before: mfem-mgis
4a007e85f2 (`install/mfem-mgis-registers`); after: bf7673f4f8
(`install/mfem-mgis-jacobians`, inverses of the Jacobians computed once;
the merged version, 7af829b3a5, only changes the documentation, the names
of two functions and a test). Quiet machine.

Times in seconds (`cases/summary-profiling.py`, total with the setup;
"time" of `summary.txt`: loop over the steps):

| Case | Loop over the steps | Integration | Tangent | Solution | Prediction |
|---|---|---|---|---|---|
| beam, before | 3.18 | 0.05 | 0.18 | 2.33 | 0.56 |
| beam, after | 3.04 | 0.05 | 0.14 | 2.23 | 0.55 |
| thermal, before | 1.32 | 0.28 | 0.39 | 0.56 | |
| thermal, after | 1.08 | 0.26 | 0.18 | 0.56 | |

- Thermal: tangent ×2.2, loop over the steps ×1.22.
- Beam: tangent ×1.3; the solution (conjugate gradient and BoomerAMG of
  hypre) dominates, total gain 4%.
