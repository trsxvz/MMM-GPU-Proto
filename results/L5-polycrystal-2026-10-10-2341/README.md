# AmgX against BoomerAMG on the polycrystal (10/10/2026)

Polycrystal of mm-opera-hpc: 5 grains, single crystal behaviour of UO2,
periodic problem, uniaxial tension, 20,553 unknowns, on the GPU.
mfem-mgis de2e76a913, all the kernels. Script:
`cases/measure-L5-polycrystal-amgx.sh` (option
`--elastic-operator-linear-solver AmgX` of the polycrystal, parameters of
AmgX by `--amgx-parameters`, single precision by
`--amgx-single-precision`). One run per case, load 3.

NewtonK0 with Eisenstat-Walker 1 (initial tolerance 0.1), the multigrid of
the elastic operator being BoomerAMG or AmgX. "Tuned":
`strength_threshold=0.8, interpolator=MULTIPASS`, the best setting on the
beam.

| Multigrid of K0 | Solve | Integrations | Linear solves |
|---|---|---|---|
| Newton (reference, PCG of hypre and BoomerAMG) | 11.66 s | | 98, 5.71 s |
| **BoomerAMG** | **9.04 s** | 225 | 164, 0.78 s |
| default AmgX | 13.60 s | 274 | 213, 3.34 s |
| tuned AmgX | 10.77 s | 245 | 184, 1.41 s |
| tuned AmgX, single | 10.70 s | 245 | 184, 1.19 s |

Macroscopic stresses identical to those of Newton (relative gap at most
3.6e-6).

The integration of the behaviour dominates (75 % of the time): a less
efficient multigrid adds Newton iterations, hence integrations, and AmgX
loses 1.7 s, while the linear solves themselves only weigh 0.8 to 1.4 s.
On this case, BoomerAMG stays the best choice. On the beam, where the
linear solve dominates, tuned AmgX in single precision gains 12 to 15 %
(results/L5-beam-2026-10-10-2259).
