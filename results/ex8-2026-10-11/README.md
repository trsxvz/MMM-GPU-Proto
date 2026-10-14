# ex8 point by point and everything by kernels, CPU and GPU (11/10/2026)

Script: `cases/measure-ex8.sh` (case `ex8` of `cases/measure-common.sh`).
Example ex8 of mfem-mgis-examples (thermomechanical plate: thermal with a
Robin condition on boundary faces, mechanics with swelling, iterative
coupling), original mesh (8424 hexahedra), 5 steps, MUMPS for the
mechanics, GMRES and BoomerAMG for the thermal. Separate build
`build/examples-others` linked to mfem-mgis 517e33cde4
(`install/mfem-mgis-others`), mfem-mgis-examples 3a907f1. CPU: 2
single-threaded processes, as in the test of the example; GPU: one process,
`--device cuda`, behaviours of the example compiled with the `generic`
interface (integrated on the host). Two runs (`summary.txt`,
`repetition-2`); during the first one, the load average was still 6 (end
of a previous attempt), during the second one 1; both agree within 2%.

Coupling over the 5 steps, in seconds:

| Case | Point by point | Everything by kernels | Gain |
|---|---|---|---|
| CPU, 2 processes | 30.97 and 30.62 | 22.35 and 22.25 | ×1.38 to ×1.39 |
| GPU, 1 process | 167.66 and 166.04 | 38.32 and 37.59 | ×4.4 |

Phases of the first run (`cases/summary-profiling.py`, two problems
summed):

| Case | Integration | Residual | Tangent | Factorisation (MUMPS) and preconditioner | Solution |
|---|---|---|---|---|---|
| CPU, point by point | 13.01 | 2.66 | 5.77 | 7.74 | 1.76 |
| CPU, kernels | 10.19 | 0.33 | 2.26 | 7.75 | 1.80 |
| GPU, point by point | 66.87 | 48.47 | 36.33 | 11.94 | 4.00 |
| GPU, kernels | 19.53 | 1.61 | 1.44 | 11.79 | 3.92 |

- Same statistics of the fields (temperature, norm of the displacement,
  swelling) in the four cases, and same numbers of Newton iterations (64:
  35 in thermal, 29 in mechanics, for 114 evaluations of the residual).
- On CPU, the integration of the mechanical behaviour (U3Si2) and the
  factorisation of MUMPS dominate; the residual (×8) and the tangent
  (×2.6) by kernels no longer weigh.
- On the GPU, point by point is much slower than on CPU (residual 48.5 s
  against 2.7 s, cause not analysed); with everything by kernels, what
  remains is the integration of the behaviours on the host (19.5 s) and
  MUMPS (11.8 s). The next lever is to compile the behaviours of ex8 for
  the GPU (`generic-parallel` interface).
- The mesh refined once is out of reach for this measurement: more than 20
  min per step with MUMPS on 8 processes, and 1,100 s of PCG with
  BoomerAMG (systems options) for 5 steps on 2 processes.
