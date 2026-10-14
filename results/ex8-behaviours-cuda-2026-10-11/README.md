# ex8 entirely on the GPU, behaviours compiled for the GPU (11/10/2026)

Script: `cases/measure-ex8.sh`. Example ex8 of mfem-mgis-examples
(thermomechanical plate: thermal with a Robin condition, mechanics with
swelling, iterative coupling), original mesh (8424 hexahedra), 5 steps,
MUMPS for the mechanics, GMRES and BoomerAMG for the thermal. Separate
build `build/examples-pointwise-ex8` (ex8 alone, mfem-mgis-examples
6505b78) linked to mfem-mgis 8011fe9eaa (`install/mfem-mgis-pointwise`:
point-wise model integrated on the device). Behaviours of the example
compiled with the `generic` interface (integrated on the host), or with the
`generic-parallel` interface for the GPU (cuda backend, cases/ex8,
integrated by batches on the GPU, swelling model included). CPU: 2
single-threaded processes; GPU: one process. Quiet machine (load 0.95),
GPU free.

Coupling over the 5 steps, in seconds:

| Case | Coupling |
|---|---|
| CPU, 2 processes, point by point | 30.51 |
| CPU, 2 processes, everything by kernels | 22.13 |
| GPU, point by point | 165.99 |
| GPU, everything by kernels, `generic` behaviours | 37.78 |
| GPU, everything by kernels, `generic-parallel` behaviours (cuda) | 18.45 (repeated: 18.53) |

Phases (`cases/summary-profiling.py`, two problems summed):

| Case | Integration | Residual | Tangent | Factorisation (MUMPS) and preconditioner | Solution |
|---|---|---|---|---|---|
| GPU, `generic` behaviours | 19.57 | 1.59 | 1.40 | 11.63 | 3.56 |
| GPU, cuda behaviours | 0.66 | 1.51 | 1.11 | 11.61 | 3.52 |

- Same statistics of the fields as point by point on CPU (temperature,
  norm of the displacement, swelling), same numbers of Newton iterations
  (64 iterations, 114 evaluations of the residual).
- With the behaviours compiled for the GPU, the integration goes from 19.6
  to 0.66 s (×30); ex8 runs in 18.5 s on the GPU, against 166 s point by
  point on the GPU and 22.1 s everything by kernels on 2 CPU processes.
- What remains is MUMPS (factorisation and solutions on the host): 15.1 s
  out of 18.6 s.
- Without the synchronisation of the mechanical gradients added to
  FieldUpdaterModel (mfem-mgis-examples 6505b78), the thermal reads stale
  gradients on the host and the statistics differ from the reference by
  about 3e-5 in relative terms (attempt on the coarse mesh).

Note added on the evening of 11/10: these computations on the GPU ran with
the host memory of Umpire, the default choice of MFEM compiled with Umpire,
which slows down the loops done on the host. With the standard host memory
(`results/host-memory-2026-10-11`): point by point 164.6 -> 57.3 s;
everything by kernels with the cuda behaviours 18.6 -> 16.9 s.
