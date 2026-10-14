# Polycrystal in float: gap to double and limits (14/10/2026)

Polycrystal of mm-opera-hpc (`uniaxial-polycrystal`, mesh `5crystals.msh`,
behaviour `Mono_UO2_Cosh_Jaco3` in finite strain), a single step of
duration d (0.1, 0.4 or 1), point by point integration on the CPU, two
processes. Float: build/mm-opera-hpc-float-gpu (environment
mmm-float-gpu). Double: build/mm-opera-hpc-9b. The behaviour is generated
by the mfront of the prototype with the `generic` interface, in float
(`--dsl-option=real_type:float`) or in double, and passed by `--library`.
Variants of the behaviour: as written (`ref`), with `cosh(x) - 1` and
`sinh(x)` computed from `expm1(x)` without cancellation (`expm1`), and
with `@Epsilon 1e-6` instead of 1e-11 (`eps6`, float only). `.res` files:
time, F11, F22, F33, then σ11, σ22, σ33 (Pa) at the end of the
computation.

| d | Double σ33 (Pa) | Float σ33 (Pa), behaviour as written | Gap | Remark |
|---|---|---|---|---|
| 0.1 | 8,997,894 | 9,053,098 | +6.1e-3 | F33 - 1 is 4.9949e-5 in float instead of 5e-5 (rounding of F) |
| 0.4 | 35,338,761 | 35,373,936 | +1.0e-3 | |
| 1 | 60,410,661 | stop at 0.86 s | | repeated failures of the integration of the behaviour (51 messages on the two processes) and two failures of the fixed point at 1e4 Pa (see the end), step cut down to 0.002 s, maximum number of substeps reached; same behaviour in a first attempt (stop at 0.85 s) |

- The gap in float decreases when the strain increases (6.1e-3 at 5e-5 of
  axial strain, 1.0e-3 at 2e-4).
- The `expm1` form gives the same result at d = 0.1 and almost the same at
  d = 0.4 (gap 1.0e-3): the cancellation of `cosh(x) - 1` is not the
  cause. In double, it changes nothing. At d = 1, it fails too (stop at
  0.53 s, with one failure of the fixed point at 1e4 Pa, see the end).
- With `@Epsilon 1e-6`, the behaviour no longer stops, but the gap at
  d = 0.1 becomes -1.5e-2, and at d = 1 the fixed point on the macroscopic
  stress stagnates at 78,000 Pa, above its threshold of 1e4 Pa (stop at
  0.5 s).

Interpretation: F is stored in float near the identity. It only moves by
steps of one ulp (6e-8 below 1, 1.2e-7 above), that is 1e-3 of a strain
of 5e-5, and a stress quantum of about 1.5e4 Pa for a modulus of
2.5e11 Pa, of the same order as the tolerance of the fixed point of these
measurements (1e4 Pa).
The imposed strain itself is rounded. These measurements argue for an
interface in H = F - I (displacement gradient, small, keeping all the
relative precision of float) for the finite strain behaviours in float;
decision to be taken (see docs/plan-precision.md).

These measurements were done with a tolerance of 1e4 Pa for the fixed
point on the lateral stresses. Since mm-opera-hpc 05358bd, it takes at
least twice its resolution, that is 30,804 Pa in float here
(results/cermet-float-2026-10-14). Redone with this rule (files
`*-resolution.res`):

- d = 0.1 and 0.4: bitwise identical files (the fixed point stops at the
  same iterations);
- d = 1, behaviour as written: stop at 0.82 s, by the failures of the
  integration of the behaviour alone (50 messages), without failure of the
  fixed point (before, two failures of the fixed point were added to
  them); `expm1` form: stop at 0.57 s, same cause;
- `@Epsilon 1e-6`: at d = 0.1, the step of 0.1 s still fails and it is
  cut into two steps of 0.05 s; the second one is no longer cut in turn
  (the fixed point stops at 12,316 Pa, below the tolerance), and the gap
  becomes -1.3e-2. At d = 1, the fixed point still stagnates around
  69,660 Pa, above the tolerance (stop at 0.66 s): a true
  non-convergence, due to the loose criterion of the behaviour (the stress
  only moves by about 20 Pa when F moves by at least two ulps, and the
  microscopic equilibrium does no Newton iteration, its residual of 55,165
  being below the absolute tolerance of 6e4).

The conclusion does not change: with a step of 1 s, the behaviour fails in
float.
