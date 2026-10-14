# Other integrators of the form in the computations by kernels (11/10/2026)

Script: `cases/measure-other-integrators.sh`. Thermal r5 (279,873
unknowns, 5 steps), 8 single-threaded processes, everything by kernels,
assembled tangent, GMRES and BoomerAMG. Second end: convective exchange
(option `-hc 1`, boundary face integrator of the nonlinear form, like the
Robin condition of the example ex8) or imposed temperature (without any
other integrator, for comparison). Before: mfem-mgis e111df8454
(`install/mfem-mgis-global`, loops of MFEM over all the elements as soon
as another integrator is present); after: 517e33cde4
(`install/mfem-mgis-others`). Quiet machine, three runs (`summary.txt`,
`repetition-2`, `repetition-3`).

Loop over the steps (line `time`), in seconds, median and extremes of the
three runs:

| Case | Before | After |
|---|---|---|
| convective exchange | 4.17 (4.10 to 4.17) | 2.70 (2.66 to 2.96) |
| imposed temperature | 2.98 (2.86 to 3.04) | 2.86 (2.85 to 2.91) |

Phases of the first run (`cases/summary-profiling.py`, number of calls in
parentheses):

| Case | Residual | Tangent | Solution |
|---|---|---|---|
| convective exchange, before | 0.49 (20) | 1.62 (15) | 1.53 |
| convective exchange, after | 0.25 (20) | 0.49 (15) | 1.68 |
| imposed temperature, after | 0.23 (23) | 0.47 (18) | 1.57 |

- With the convective exchange, the tangent goes from 108 to 33 ms per
  call (×3.3) and the residual from 24.5 to 12.5 ms (×2); loop over the
  steps ×1.54 (medians). Per call, about 25% more than without any other
  integrator remains (tangent 33 against 26 ms, residual 12.5 against
  10 ms): evaluation of the boundary faces on the host and prolongation of
  the unknowns.
- Same norms of the solution and same Newton iterations before and after
  (15 for the convective exchange, 18 for the imposed temperature).
- Without any other integrator, the path does not change: the gaps of the
  imposed temperature (less than 4%) are measurement noise.
- The linear solution now dominates the computation (1.5 to 1.7 s out of
  2.7 to 3 s).
