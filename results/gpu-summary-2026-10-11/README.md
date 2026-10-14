# Summary on the GPU of the work packages of 11/10/2026 (beam r4)

Script: `cases/measure-gpu-summary.sh`. Beam r4 (111,843 unknowns, 10
steps), one process and the RTX 5080, everything by kernels,
`generic-parallel` behaviour compiled for CUDA. Before: mfem-mgis
8a6325335e (`install/mfem-mgis`); after: e111df8454
(`install/mfem-mgis-global`: positions of the entries, matrices and action
without forming B, inverses of the Jacobians, global matrix kept; the
latter is not used with hypre on the GPU). Quiet machine.

Loop over the steps ("time" of `summary.txt`) and phases
(`cases/summary-profiling.py`), in seconds:

| Variant | Before | After | Tangent before / after | Solution before / after |
|---|---|---|---|---|
| Newton, assembled tangent | 3.17 | 2.92 | 0.22 / 0.14 | 2.25 / 2.12 |
| NewtonK0 + Eisenstat-Walker, assembled tangent | 1.90 | 1.71 | 0.22 / 0.15 | 1.11 / 1.04 |
| NewtonK0 + Eisenstat-Walker, without assembly | 2.12 | 1.87 | 0 / 0 | 1.48 / 1.26 |

- Gains of 9 to 14%: the solution (conjugate gradient and BoomerAMG of
  hypre on the GPU) dominates. Without assembly, the solution gains 15%
  through the kernel of the action (B not formed, inverses of the
  Jacobians kept).
- Same solutions within the Newton tolerance (norms identical within 1e-9,
  Newton iterations within one for the variants with adapted tolerance).
