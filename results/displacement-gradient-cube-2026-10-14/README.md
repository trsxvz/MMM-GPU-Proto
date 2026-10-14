# Displacement gradient on the cube of the uniaxial test (14/10/2026)

Test `UniaxialTensileTest` of mfem-mgis (cube of one element, cyclic
tension and compression up to 2.1 %, 100 steps, one process), in double
and in single precision, with the deformation gradient F (default) and
with the displacement gradient H = F - I (option
`--displacement-gradient`, parameter `DisplacementGradient` of the
behaviour integrators). Two finite strain behaviours of the tests of
mfem-mgis:

- `LogarithmicStrainElasticity`: isotropic elasticity in the Hencky strain
  framework, whose strain is computed from H without cancellation;
- `SaintVenantKirchhoffElasticity`: written with the deformation gradient,
  which receives F = I + H (no gain expected).

Files: first two components of the deformation gradient, first component
of the first Piola-Kirchhoff stress and equivalent strain at the first
integration point, at each step (format of the reference files of
mfem-mgis).

## Double

`double-LogarithmicStrainElasticity-{F,H}.txt`: same results with F and H,
within 1e-13 on the gradients and 7e-14 relative on the stress (the
reference `LogarithmicStrainElasticity.ref` of mfem-mgis is the run with
F).

## Float, default tolerance of the test (1e-3 relative)

`float-*-{F,H}.txt`, compared with the double reference:

| Behaviour | F | H |
|---|---|---|
| LogarithmicStrainElasticity | 3.9e4 Pa (4.4e-4) | 4.2e4 Pa (3.2e-4) |
| SaintVenantKirchhoffElasticity | 4.9e4 Pa (5.5e-4) | 2.8e4 Pa (3.6e-4) |

Largest gap on the stress (largest relative gap where the stress exceeds
1e7 Pa). At this tolerance, the gap is set by the convergence criterion of
the Newton solver, not by the strain.

## Float, reachable tolerance

`tolerances.txt`: Newton solver with a relative tolerance of 1e-4, 1e-5
and 1e-6, the acceptance of a stagnation of the residual being disabled
(`--stagnation-relative-tolerance 0`).

- `LogarithmicStrainElasticity` with H converges at every step down to
  1e-8 (`tolH-1e-{6,7,8}.txt`): gap to double 4.7e3 Pa (4.8e-5 relative),
  ten times smaller than with the default tolerance.
- With F, it fails from 1e-4 on (first steps): the residual stagnates at
  the floor set by F, known within 6e-8.
- `SaintVenantKirchhoffElasticity` fails in both cases: it is written with
  F, which the interface computes as I + H.
