# Cermet in float: tolerance of the fixed point on the lateral stresses (14/10/2026)

Cermet of mm-opera-hpc (test `CermetNullStrainXXYY`: mesh `5grains.msh`,
two steps of 0.2 s, no post-processing). The behaviour is passed by
`--library`: float library compiled by gcc 13 (float build of the
prototype) or by gcc 14 (simulated CI), in both cases with the options of
`tfel-config --oflags`. Executables:

- "CI": cermet of the CI simulated locally (`build/ci-sim-float`,
  environment `ci-cpu-float`; `build/ci-sim-clean` in double), in an empty
  environment;
- "prototype": cermet of the float build of the prototype
  (`build/mm-opera-hpc-float`, environment `mmm-float-cpu`).

Files `<executable>-<n>p-<behaviour>-<tolerance>.res` (n processes) and
`double-*.res`: time, F11, F22, F33, then σ11, σ22, σ33 (Pa) at the end of
each step. `fixed-point.txt`: residual of the fixed point,
|res| = (σ11² + σ22²)^½, at each iteration of each computation.

## Rule adopted

The tolerance on the lateral stresses is 1e4 Pa, and the fixed point
takes at least twice its resolution (mm-opera-hpc 05358bd and ae5768b).
When the two lateral stresses are close to the same value S, the
correction of the lateral components of F is about
(1 + ν) (1 - 2 ν) S / E. It rounds to zero when it is smaller than half
the gap between two neighbouring floats of these components. With the
effective Kröner moduli of the UO2 grains (E = 173 GPa, ν = 0.321), used
by the cermet and the polycrystal, and lateral components just below 1
(gap of 6e-8), the resolution is 1.54e4 Pa, and the tolerance 30,804 Pa.
In compression, the lateral components exceed 1 and the resolution
doubles. In double, it is about 3e-5 Pa: nothing changes (bitwise
identical files).

The files `*-resolution.res` (rule adopted) are bitwise identical to the
files `*-3e4.res` of the same configuration.

## Sweep of the constant tolerances

Double, executable of the CI, one process (the test): σ33 = 17,219,167 Pa
at 0.2 s and 29,377,289 Pa at 0.4 s with steps of 0.2 s; 17,642,604 Pa
(+2.46 %) and 29,886,489 Pa (+1.73 %) with steps of 0.1 s.

Float: σ33 at 0.2 s then at 0.4 s, and gap to the double with steps of
0.2 s.

| Executable | Processes | Tolerance (Pa) | Behaviour compiled by gcc 13 | Behaviour compiled by gcc 14 |
|---|---|---|---|---|
| CI | 1 | 1e4 | 17,684,464 (step cut), +2.70 %; 29,501,400, +0.42 % | 17,684,506 (step cut), +2.70 %; 29,501,338, +0.42 % |
| CI | 1 | 2e4 | 17,247,006, +0.16 %; 29,387,932, +0.04 % | 17,247,020, +0.16 %; 29,387,976, +0.04 % |
| CI | 1 | 3e4 | 17,247,006, +0.16 %; 29,396,894, +0.07 % | 17,247,020, +0.16 %; 29,396,888, +0.07 % |
| CI | 1 | 5e4 | 17,262,850, +0.25 %; 29,397,188, +0.07 % | 17,262,862, +0.25 %; 29,397,298, +0.07 % |
| CI | 2 | 1e4 | 17,684,538 (step cut), +2.70 %; 29,501,340, +0.42 % | 17,684,532 (step cut), +2.70 %; 29,501,422, +0.42 % |
| CI | 2 | 2e4 | 17,246,848, +0.16 %; 29,387,992, +0.04 % | 17,246,870, +0.16 %; 29,387,964, +0.04 % |
| CI | 2 | 3e4 | 17,246,848, +0.16 %; 29,396,982, +0.07 % | 17,246,870, +0.16 %; 29,396,966, +0.07 % |
| CI | 2 | 5e4 | 17,262,702, +0.25 %; 29,397,332, +0.07 % | 17,262,718, +0.25 %; 29,397,140, +0.07 % |
| prototype | 2 | 1e4 | 17,245,144, +0.15 %; 29,393,436, +0.05 % | 17,690,582 (step cut), +2.74 %; 29,498,396, +0.41 % |
| prototype | 2 | 3e4 | 17,237,086, +0.10 %; 29,393,648, +0.06 % | 17,237,132, +0.10 %; 29,393,532, +0.06 % |

Reference of Cast3M for the test: 17,246,000 Pa at 0.2 s and
29,462,000 Pa at 0.4 s, within 1 %.

- With 1e4 Pa, the fixed point of the first step does not converge in five
  computations out of six: it stagnates around 10,684 to 10,830 Pa with
  the executable of the CI, and oscillates between 17,234.6 and
  15,454.4 Pa with the executable of the prototype and the behaviour of
  gcc 14; with this same executable, the behaviour of gcc 13 converges.
  The rounding of the whole computation, behaviour included, decides the
  point where the fixed point stops. The step is rejected and then cut
  into two steps of 0.1 s: the gap of 2.7 % to the double comes from this
  time step (with steps of 0.1 s, the float departs from the double by
  0.24 to 0.27 %).
- 2e4 Pa works here, but only 1.16 times above the highest oscillation
  measured (17,235 Pa). 5e4 Pa leaves larger lateral stresses (+0.25 % at
  0.2 s).

## Gap to the double with the rule adopted

+0.10 to +0.16 % at 0.2 s and +0.06 to +0.07 % at 0.4 s. At 0.2 s, on the
computations of the CI (+0.16 %), roughly:

- +0.07 %: absolute tolerance of the microscopic equilibrium in float
  (6e4), measured in double with this tolerance
  (`double-microscopic-tolerance-6e4.res`: +0.072 % at 0.2 s, +0.023 % at
  0.4 s);
- +0.03 %: lateral stresses left by the fixed point, which shift σ33 by
  ν (σ11 + σ22) (-0.045 % on the computations of the prototype);
- +0.017 %: rounding of the imposed strain (1 + 1e-4 is 1.000100017 in
  float);
- the rest, about +0.04 %, is the rounding of the computation in float.

## Polycrystal

It uses the same class and the same moduli. Its measurements in float
(results/precision-polycrystal-2026-10-14), redone with the rule adopted:
bitwise identical files for steps of 0.1 and 0.4 s; for the others, see
the README of that directory.
