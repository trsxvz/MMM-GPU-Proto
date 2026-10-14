# Float against double precision on the GPU, beam (14/10/2026)

Beam r4 of the prototype (111,843 unknowns), 10 steps, Newton (criterion
1e-4, and the one of the case, 1e-8, in double), PCG of hypre
preconditioned by BoomerAMG (tolerance 1e-4), everything on the GPU
(`--execution gpu`, behaviour of the `cuda` backend), one process, NVIDIA
RTX 5080 GPU free (the neighbouring session paused). Script:
`DEVICE=gpu cases/measure-precision-beam.sh`. Double: build/cases-9b
(mfem-mgis gpu-proto installed separately, environment mmm-gpu). Float:
build/cases-float-gpu (environment mmm-float-gpu, MFEM and hypre in single
precision with CUDA). Same version of mfem-mgis.

| Computation | Criterion | Min. time (s) | Median (s) | Iterations | Norm of the displacements |
|---|---|---|---|---|---|
| double | 1e-8 | 2.91 | 3.02 | 32 | 48.83704586 |
| double | 1e-4 | 1.86 | 1.86 | 16 | 48.83715008 |
| float | 1e-4 | 1.42 | 1.42 | 16 | 48.85159302 and 48.84963608 (2 computations out of 3) |

Time per phase (first run at 1e-4): linear solution 1.13 s in double,
0.81 s in float; prediction 0.55 and 0.49 s; the rest is small.

- Float is 1.31 times faster when it converges, the gain coming mostly
  from the linear solution.
- But the float path on the GPU is not reliable: one computation out of
  three diverged (Newton at the tenth step), the converged computations
  differ from each other (4e-5) and from double (3e-4, against 8e-5 for
  float on the CPU). Analysis (`scratch`, 14/10): 4 failures out of 8 with
  the atomic assembly, 8 out of 8 by kernels with the ordered assembly and
  an iterative solver (a computation with MUMPS converges); the
  computation converges if the finite elements are on the CPU, while the
  behaviours and the solver stay on the GPU, and diverges if the finite
  elements are on the GPU, whatever the placements of the behaviours and
  of the solver (PCG and BoomerAMG). The same kernels in float on the CPU
  converge. The convergence of Newton degrades from step to step on the
  GPU.
  In double, with the criterion of the case (1e-8), the three computations
  are bitwise identical; with the criterion 1e-4, they differ by 4e-7
  (atomic assembly).
- Diagnosis (results/diagnostic-float-gpu-2026-10-14): the contractions
  into fused multiply-adds by nvcc in the stiffness kernels made Newton
  diverge in float. Fixed in mfem-mgis (in float, stiffness kernels
  compiled for the GPU without contraction). Measurement redone:
  results/precision-beam-gpu-2026-10-14-0843 (float 1.37 times faster than
  double, 3 computations out of 3 converged, within 3.4e-5 of double).
