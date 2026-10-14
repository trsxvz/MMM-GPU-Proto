# cuDSS on the beam r4 (12/10/2026)

Beam r4 of the prototype (111,843 unknowns, von Mises with linear
hardening, 10 steps), everything by kernels, cuda behaviour, one process
on the GPU (cores 0 to 3), hypre on the GPU, elastic prediction. mfem-mgis
45afaa0653 (`install/mfem-mgis-mumps`), build `build/cases-cudss` (only the
executable `beam`, linked to this installation). Options added to the case
`beam`: `-ls cudss` (solver `CuDSSSolver`, symmetric positive definite),
`-mra` (reuse of the analysis of the direct solvers) and `-els cudss`
(elastic operator of the fixed point factorised by cuDSS). Two passes per
variant (`summary.txt`, commands in `measure.sh` and `measure-k0.sh`).
Time of the loop over the steps, in seconds:

| Method | Solver | Time |
|---|---|---|
| Newton | PCG of hypre and BoomerAMG (tolerance 1e-4) | 3.12 and 3.04 |
| Newton | cuDSS, analysis reused | 9.60 and 9.56 |
| NewtonK0 (K0 assembled once, prediction by K0) | one cycle of BoomerAMG on K0, Eisenstat and Walker tolerance | 1.87 and 1.85 |
| NewtonK0 | K0 factorised once by cuDSS, tolerance 1e-4 | 3.01 and 3.03 |
| NewtonK0 | K0 factorised once by cuDSS, Eisenstat and Walker tolerance | 2.73 and 2.97 |

- At this size, the direct solver does not pay off. With Newton, one
  factorisation of cuDSS takes 0.235 s even with the analysis reused
  (6.1 s for 26 factorisations, and 3.0 s for the 10 predictions, which
  factorise the elastic operator), against 68 ms for one solve by PCG and
  BoomerAMG. One solve by cuDSS then takes only 8 ms.
- With NewtonK0, K0 is factorised only once (0.6 s), but each tangent
  solve by the conjugate gradient preconditioned by K0 costs 67 ms
  (triangular forward and backward substitutions at each iteration of the
  conjugate gradient), against 22 ms with one cycle of BoomerAMG.
- The best variant stays NewtonK0 with one cycle of BoomerAMG and the
  Eisenstat and Walker tolerance (1.85 s; 2.14 s on 10/10 with mfem-mgis
  517c0c0b08, `results/L5-beam-2026-10-10-0429`).
- cuDSS is the right choice for small systems (ex7, 2,634 unknowns; ex8,
  9,516 and 28,548 unknowns: `results/cudss-2026-10-12`), where the
  iterative solves on the GPU are dominated by the kernel launches, not
  for large three-dimensional systems, where the fill-in of the
  factorisation costs more than the AMG.
