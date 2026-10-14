# NewtonK0 on the polycrystal (10/10/2026)

Polycrystal of mm-opera-hpc: 5 grains, single crystal behaviour of UO2
(orthotropic, finite strain, implicit), periodic problem, uniaxial
tension, 10 steps. mfem-mgis 9249013150. Script:
`cases/measure-L5-polycrystal.sh`. One RTX 5080, 16 cores. Quiet machine
at the start (load 1.09).

All the kernels (batch integration, gradients, internal forces and tangent
assembled by kernels). On CPU: 8 processes, `generic` behaviour. On GPU:
one process, `generic-parallel` behaviour compiled for `cuda`.

- **Newton**: PCG of hypre preconditioned by BoomerAMG, tolerance 1e-12,
  new setup of BoomerAMG at each Newton iteration (on the consistent
  tangent operator).
- **NewtonK0**: option `--non-linear-solver NewtonK0` of the polycrystal,
  method `NewtonK0` of the `FixedPoint` solver of mfem-mgis. Conjugate
  gradient of MFEM (tolerance 1e-12) preconditioned by one cycle of
  BoomerAMG built only once, on the elastic operator K0 obtained by one
  integration of the behaviour (`ElasticOperatorComputation` set to
  `Integration`, `ElasticOperatorUpdate` set to `Once`).

| Case | Total | Integration | Tangent | Solve |
|---|---|---|---|---|
| CPU, 8 processes, Newton | 29.60 | 13.84 (158) | 2.02 (98) | 13.42 (98) |
| CPU, 8 processes, NewtonK0 | **24.92** | 14.08 (159) | 2.04 (99) | 8.53 (98) |
| GPU, 1 process, Newton | 11.88 | 5.09 (158) | 0.20 (98) | 5.11 (98) |
| GPU, 1 process, NewtonK0 | **11.26** | 5.16 (159) | 0.20 (99) | 4.45 (98) |

Times in seconds (maximum over the processes), number of calls in
parentheses. The residual (0.07 s on CPU, 0.02 s on GPU) and the end of
the steps are negligible. Macroscopic stresses of NewtonK0 equal to those
of Newton in the same configuration, up to 2.4e-13 (CPU) and 1.8e-14 (GPU)
of the stress scale.

## Findings

- **Same Newton iterations**: 98 solves in both cases. NewtonK0 does one
  more integration and one more assembly (K0, at the first step).
- **CPU**: the solve goes from 13.42 s to 8.53 s (-36 %), the computation
  from 29.60 s to 24.92 s (×1.19). The setup of BoomerAMG at each
  iteration thus weighs about one third of the solve of Newton.
- **GPU**: the solve goes from 5.11 s to 4.45 s (-13 %), the computation
  from 11.88 s to 11.26 s (×1.06). The setup of BoomerAMG is cheaper on
  the GPU. The integration of the behaviour (5.1 s, 32 ms per call) now
  dominates.
- First trial on GPU (same day, less quiet machine): 11.90 s against
  11.24 s, gap 1.3e-14: reproducible.
- As on the beam (`results/L5-beam-2026-10-10-0250/README.md`), NewtonK0
  is the only stable gain on GPU among the variants of the fixed point. It
  does not touch the cost of the behaviour.
