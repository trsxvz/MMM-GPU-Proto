# Newton solver at the floor of the residual, beam r4 in float (14/10/2026)

Beam of the prototype (`cases/beam`), 4 refinements (111,843 unknowns), 10
time steps, J2 plasticity with linear hardening, Newton, PCG of hypre
preconditioned by BoomerAMG (linear tolerance 1e-4), two processes on the
cores 0 and 1. Float build: spack environment `mmm-float-cpu` (MFEM, hypre
and MUMPS in single precision, hypre with its internal LAPACK), mfem-mgis
with the option `StagnationRelativeTolerance` of the Newton solver (default
value, 10⁴ ε = 1.2e-3). Files: output of each computation.

| Computation | Newton criterion | Time (s) | Newton iterations | Norm of the displacements | Gap to double at 1e-8 |
|---|---|---|---|---|---|
| double (reference, `results/precision-beam-2026-10-14-0514`) | 1e-8 | 34.27 | 28 | 48.83704586 | |
| float, before the stagnation criterion (`results/precision-beam-2026-10-14-0506`) | 1e-8 | | 50 at the first stalled step, no convergence | | |
| float | 1e-8 | 30.26 | 27 | 48.83361435 | -7.0e-5 |
| float | 1e-4 | 20.79 | 16 | 48.83312988 | -8.0e-5 |

- With the criterion of the case (1e-8, bounded to 100 ε = 1.2e-5 in
  float), the residual reaches its floor in two or three iterations: 1.6e-5
  to 2.0e-5 of the reference residual, growing with the load. The Newton
  solver accepts this stagnation after two iterations without halving the
  smallest residual (5 steps out of 10), and the computation converges.
- The solution is slightly closer to the solution in double than with a
  criterion of 1e-4 (7.0e-5 against 8.0e-5 on the norm of the
  displacements).
- In double, the criterion is not triggered: no stagnation message in the
  logs of the test suites of mfem-mgis, of the examples and of mm-opera-hpc
  (CI simulation), and the results of mm-opera-hpc are the same bit for
  bit.
