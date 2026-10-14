# Float against double precision on the GPU, beam, after the fix (14/10/2026)

Same measurement as results/precision-beam-gpu-2026-10-14-0604, after the
fix of mfem-mgis: in float, the kernels of the stiffness matrices are
compiled for the GPU without contraction into fused multiply-adds; in
double, the code is unchanged (see results/diagnostic-float-gpu-2026-10-14).

Beam r4 of the prototype (111,843 unknowns), 10 steps, Newton (criterion
1e-4, and the one of the case, 1e-8, in double), PCG of hypre
preconditioned by BoomerAMG (tolerance 1e-4), everything on the GPU
(`--execution gpu`, atomic assembly, behaviour of the `cuda` backend), one
process on core 0, NVIDIA RTX 5080 GPU free. No other computation during
the measurement (the load of 3.5 over the previous minute comes from the
test suites finished just before). Script:
`DEVICE=gpu cases/measure-precision-beam.sh`. Double: build/cases-9b
(environment mmm-gpu). Float: build/cases-float-gpu (environment
mmm-float-gpu, MFEM and hypre in single precision with CUDA).

| Computation | Criterion | Min. time (s) | Median (s) | Iterations | Norm of the displacements |
|---|---|---|---|---|---|
| double | 1e-8 | 2.89 | 2.94 | 31 to 32 | 48.83704586 |
| double | 1e-4 | 1.85 | 1.85 | 16 | 48.83715811 to 48.83716422 |
| float | 1e-4 | 1.35 | 1.35 | 16 | 48.83538818 to 48.83585739 |

Time per phase (first run at 1e-4): linear solution 1.12 s in double,
0.80 s in float; prediction 0.55 and 0.42 s; assembly of the tangent 97
and 83 ms for 16 assemblies (100 and 84 to 85 ms before the fix).

- The three float computations converge, with the same iterations as
  double. Before the fix, one out of three diverged.
- Float is 1.37 times faster than double (1.31 before the fix, on the
  computations that converged).
- The norm of the displacements in float differs by 2.7e-5 to 3.6e-5 from
  double with the same criterion, and by 3.4e-5 at most from converged
  double (3e-4 before the fix; double with the criterion 1e-4 is at 2.4e-6
  from it).
- The computations in double differ from each other by 1.3e-7, the
  computations in float by 9.6e-6: the atomic assembly does not guarantee
  the order of the additions (the ordered assembly,
  `MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered`, gives bitwise identical
  computations).
