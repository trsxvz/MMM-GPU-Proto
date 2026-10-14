# Measurements of work package L3 on the beam (09/10/2026)

Beam r4 (111,843 unknowns), PCG + BoomerAMG (tolerance 1e-4), 10 steps, on
CPU. Batch integration with the gradients, the internal forces (residual)
and the element matrices of the tangent computed element by element or by
kernels (mfem-mgis ed91564643). Script: `cases/measure-L3.sh beam`.

Caveat: an independent review was running tests at the same time. The
point by point computation with 1 process, started first, is slower
everywhere (solve included) and is not comparable. The comparisons are made
with "batches".

| Case | Total | Integration | Residual | Tangent | Solve |
|---|---|---|---|---|---|
| 8 processes, batches | 15.69 | 0.51 (41) | 0.45 (41) | 2.24 (31) | 9.56 |
| 8 processes, gradients by kernels | 15.17 | 0.20 (41) | 0.45 (41) | 2.22 (31) | 9.37 |
| 8 processes, residual by kernels | 15.34 | 0.50 (41) | 0.21 (41) | 2.23 (31) | 9.48 |
| 8 processes, tangent by kernels | 14.93 | 0.47 (38) | 0.40 (38) | 2.51 (28) | 8.61 |
| 8 processes, everything by kernels | 14.32 | 0.19 (38) | 0.20 (38) | 2.51 (28) | 8.47 |
| 1 process, batches | 65.71 | 3.18 (38) | 2.96 (38) | 11.79 (28) | 33.85 |
| 1 process, gradients by kernels | 62.48 | 0.94 (38) | 2.99 (38) | 11.64 (28) | 33.25 |
| 1 process, residual by kernels | 64.04 | 3.21 (38) | 1.40 (38) | 11.71 (28) | 33.78 |
| 1 process, tangent by kernels | 68.04 | 3.22 (38) | 2.98 (38) | 14.62 (28) | 33.38 |
| 1 process, everything by kernels | 64.12 | 1.01 (38) | 1.35 (38) | 14.55 (28) | 33.39 |

In parentheses: number of calls. Same solution everywhere (norm
48.83704586). With 8 processes, the tangent by kernels changes the rounding
and the number of Newton iterations (28 instead of 31), hence a lower total
time that does not come from the kernels.

## Findings

- **Gradients by kernels**: integration phase divided by 2.5 (8 processes)
  to 3.4 (1 process). The call of the behaviour remains.
- **Residual by kernels**: divided by 2.1, while the loop of MFEM over the
  elements still runs (it only returns zeros). It is removed when the form
  has no other integrator (next commit).
- **Tangent by kernels**: slower on CPU (+24 % per call). Computing all the
  element matrices, storing them and then copying them into the assembly of
  MFEM costs more than computing them in its loop. It will only pay off on
  GPU with an assembly on the device, still to be written
  (`ElementRestriction::FillJAndData` of MFEM only handles scalar spaces),
  or matrix free (option a of the plan).
- **First trial, before the fix** (directory `scratch/`): the kernels used
  `Mesh::GetGeometricFactors`, which gives nodes to the mesh
  (`EnsureNodes`). All the element transformations then went through these
  nodes, and the phases left element by element were 50 % slower. The
  geometric factors are now computed on a copy of the vertices.
