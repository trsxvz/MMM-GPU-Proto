# Displacement gradient in the cermet and the polycrystal (14/10/2026)

Cermet (`CermetNullStrainXXYY` setup: mesh `5grains.msh`, duration 0.4 s,
two steps of 0.2 s or four steps of 0.1 s, no post-processing) and
polycrystal (`uniaxial-polycrystal`, mesh `5crystals.msh`, one step of
duration d = 0.1, 0.4 or 1 s) of mm-opera-hpc, two processes, in double
and in single precision, with the deformation gradient (F, default) and
with the displacement gradient H = F - I (option `--displacement-gradient`).
Script: `measure-opera.sh` in this folder.

Versions: TFEL 688c9c2cb, MGIS 9b1a03a, mfem-mgis be675e2c54,
mm-opera-hpc 457549d (the final commits, TFEL 5522bc432, MGIS 5a28283,
mfem-mgis a42e0038da and mm-opera-hpc 1322426, only add checks, tests
and documentation). With the option:

- the behaviours are given H. The metal of the cermet (`NortonCr`, Hencky
  strain framework) computes its strain from H without cancellation. The
  grains of both cases (`Mono_UO2_Cosh_Jaco3`, written with the
  deformation gradient) receive F = I + H, computed in float by TFEL;
- the macroscopic unknowns of mm-opera-hpc are H in both modes; with the
  option, mfem-mgis is given the macroscopic H, otherwise F = I + H
  rounded in float. The resolution of the fixed point on the lateral
  stresses is computed on the gradients given to mfem-mgis.

Files `<precision>-<case>-<steps>-<F|H>.res`: time, F11, F22, F33 (computed
in double from H), then σ11, σ22, σ33 (Pa) at the end of each step;
`.log`: output of the run.

## Double

H and F give the same results within 1e-12 relative, and the same values
as before the change (cermet 29,377,289 Pa at 0.4 s with two steps,
29,886,489 Pa with four steps; polycrystal 8,997,894, 35,338,761 and
60,410,661 Pa).

## Float

Gap of σ33 to the double at the end of the computation (Pa, relative):

| Case | F | H |
|---|---|---|
| cermet, 2 x 0.2 s | +24,941 (8.5e-4) | +8,243 (2.8e-4) |
| cermet, 4 x 0.1 s | +69,559 (2.3e-3) | +70,429 (2.4e-3) |
| polycrystal, d = 0.1 | +61,832 (6.9e-3) | +53,197 (5.9e-3) |
| polycrystal, d = 0.4 | +35,179 (1.0e-3) | +14,459 (4.1e-4) |
| polycrystal, d = 1 | stop at 0.877 s | stop at 0.848 s |

- Fixed point on the lateral stresses: with F, its tolerance is raised to
  twice its resolution (30,804 Pa) in 25 time steps of these runs; with H,
  never. The lateral stresses at the end are 2 to 30 times smaller with H
  (cermet with two steps: 1.7e3 and -1.1e3 Pa instead of 2.8e4 and
  -1.2e4 Pa; polycrystal at 0.4 s: 107 and -141 Pa instead of 2.3e3 and
  4.4e3 Pa).
- Cermet with four steps: both modes depart from the double by about 4e4
  Pa at the second step already (17,682,784 Pa with F, 17,678,366 Pa with
  H, 17,642,604 Pa in double at 0.2 s). The kinematics is not the cause:
  the grains, written with F and viscoplastic, are integrated in float
  from F = I + H, and the microscopic equilibrium stops at its absolute
  tolerance of 6e4 in float (stagnation of the residual set by the
  rounding of F in the grains).
- Polycrystal at d = 0.1: the imposed axial strain is now exact (H33 =
  5e-5), but the behaviour of the grains receives F33 = fl(1 + 5e-5), that
  is a strain of 4.9949e-5: the gap stays close to the one with F. The
  gain requires a behaviour written with the displacement gradient or in
  a strain framework.
- Polycrystal at d = 1: in both modes, repeated failures of the
  integration of the behaviour in float (80 messages with F, 59 with H),
  nine cuts of the time step, stop at the maximum number of substeps
  (known limit of this behaviour in float, see
  `results/precision-polycrystal-2026-10-14`).

## Conclusion

The displacement gradient removes the quantization of the macroscopic
loading and of the fixed point, and the cancellation in the behaviours
written in the Green-Lagrange or Hencky strain frameworks. In these two
cases, the float error is then set by the grains, written with the
deformation gradient: rewriting this behaviour with H (or in a strain
framework) is the next step to benefit from it.
