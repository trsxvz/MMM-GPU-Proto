# Newton solver at the floor of the residual, beam r4 in float (14/10/2026)

Beam of the prototype (`cases/beam`), 4 refinements (111,843 unknowns), 10
time steps, J2 plasticity with linear hardening, Newton, PCG of hypre
preconditioned by BoomerAMG (linear tolerance 1e-4), two processes on the
cores 0 and 1, quiet machine. Float build: spack environment
`mmm-float-cpu` (MFEM, hypre and MUMPS in single precision, hypre with its
internal LAPACK), mfem-mgis 9c2dc8c4e5 with the option
`StagnationRelativeTolerance` of the Newton solver (default value, 10⁴ ε =
1.2e-3): a residual which has not decreased by more than 10 % during the
last two iterations, below this tolerance relative to the reference
residual, is accepted as converged, and the best iterate of the
resolution is kept. Files: output of each computation.

| Computation | Newton criterion | Time (s) | Newton iterations | Norm of the displacements | Gap to double at 1e-8 |
|---|---|---|---|---|---|
| double (reference, `results/precision-beam-2026-10-14-0514`) | 1e-8 | 34.27 | 28 | 48.83704586 | |
| float, before the stagnation criterion (`results/precision-beam-2026-10-14-0506`) | 1e-8 | | 50 at the first stalled step, no convergence | | |
| float | 1e-8 | 30.32 | 27 | 48.83361435 | -7.0e-5 |
| float | 1e-4 | 20.92 | 16 | 48.83312988 | -8.0e-5 |

- With the criterion of the case (1e-8, bounded to 100 ε = 1.2e-5 in
  float), the residual reaches its floor in two or three iterations: 1.36e-5
  to 2.02e-5 of the reference residual, growing with the load. The Newton
  solver accepts this stagnation in 5 steps out of 10, without restoring a
  previous iterate (the last one is the best one), and the computation
  converges. The first version of the criterion (mfem-mgis 11fda6426e,
  residual not halved) gave the same iterations and the same solution.
- The solution is slightly closer to the solution in double than with a
  criterion of 1e-4 (7.0e-5 against 8.0e-5 on the norm of the
  displacements).
- In double, the criterion is not triggered: no stagnation message in the
  logs of the test suites of mfem-mgis, of the examples and of mm-opera-hpc
  (CI simulation), and the results of mm-opera-hpc are the same bit for
  bit.
