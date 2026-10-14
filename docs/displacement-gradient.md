# Displacement gradient for the finite strain behaviours

## Problem

In single precision, the deformation gradient F = I + H, close to the
identity, is only known within one ulp of 1 (6e-8): at a strain of 5e-5,
the displacement gradient H is only known within 1e-3 in relative value.
The strain measures computed from F lose as much (C - I = Fᵀ F - I by
cancellation). Measured effects: polycrystal 0.6 % away from the double at
5e-5 of strain, floor of the residual of 3e-4 to 8e-4 on a finite strain
cube, quantized fixed point on the lateral stresses of the cermet.

## Principle

The solver may give the displacement gradient H = F - I instead of F. The
behaviour computes its strain measure from H without cancellation:

- Green-Lagrange strain framework: E = (H + Hᵀ + Hᵀ H) / 2;
- Hencky strain framework: eigenvalues of C - I = 2 E, then
  log(1 + λ) / 2 by `log1p`;
- behaviours written with F (DSLs `DefaultFiniteStrain`,
  `ImplicitFiniteStrain`...): they receive F = I + H computed by the
  interface, as before (no gain, no loss).

F = I + H is still used by the conversions of the stresses and of the
tangent operators, which do not suffer from the cancellation. The tangent
operators are unchanged: the derivatives with respect to F are the
derivatives with respect to H.

## Encoding

The option is chosen at runtime by the caller, as the stress measure: the
second value of the `K` array (stress measure: 0 Cauchy, 1 PK2, 2 PK1) is
offset by 10 when the gradients are the displacement gradient, as the
first value is offset by 100 when the speed of sound is requested.

- New library, old caller: the offset is never given, nothing changes.
- Old library, new caller asking for H: the value is an invalid stress
  measure for the old library, which reports an error (no silent misuse).
  MGIS refuses the option earlier, with a clear message, when the library
  does not export the symbol `<function>_AcceptsDisplacementGradient`.
- The API version of the generic interface is unchanged (MGIS checks it
  strictly: a change would refuse all the libraries generated before).

## Plane stress

TFEL expects F with a zero axial component (computed by the behaviour),
whereas MGIS and mfem-mgis give 1 with the deformation gradient (mismatch
older than this work, to be reported upstream). With H, the caller gives a
zero axial component and the identity is not added to it: F = I + H has a
zero axial component, as TFEL expects, and the standard and Green-Lagrange
frameworks add the axial component computed by the behaviour. The axial
Green-Lagrange strain given to the behaviour is zero, as with MGIS and the
deformation gradient. The Hencky framework sets the axial component after
building the logarithmic strain handler: the handler is built with
F_zz = 1 and E_zz = 0.

## Implementation

Patches `patches/tfel/0121` to `0124`, `patches/mgis/0071` to `0073`,
`patches/mfem-mgis/0182` and `patches/mm-opera-hpc/0034` (TFEL 5522bc432,
MGIS 5a28283, mfem-mgis a42e0038da, mm-opera-hpc 1322426).

TFEL:

- `computeGreenLagrangeTensorFromDisplacementGradient` (TFEL/Math);
- `LogarithmicStrainHandler` built from F and from the Green-Lagrange
  strain; its divided differences of the logarithm are computed with
  `log1p` from the eigenvalues of C - I, also when it is built from F
  (the former formulas lost all precision in float when two eigenvalues
  are close: 1e-1 of error on the second Piola-Kirchhoff stress and 4.6e2
  on its derivative at 10 % of strain, with a stress not coaxial with the
  strain; now below 4e-6);
- `MFront/GenericBehaviour`: decoding of the option
  (`ComputationOptions::displacement_gradient`), F = I + H and the
  Green-Lagrange strain from H (`DisplacementGradient.hxx`), for the
  integration, the initialize functions and the post-processings of the
  three kinds of finite strain behaviours, `generic` and
  `generic-parallel` interfaces; the MTest file generated on failure
  contains F = I + H;
- symbol `<function>_AcceptsDisplacementGradient`.

MGIS:

- `FiniteStrainBehaviourOptions::gradient` (`DEFORMATION_GRADIENT` by
  default, `DISPLACEMENT_GRADIENT`), `usesDisplacementGradient`;
- `load` checks the symbol, offsets the stress measure and names the
  gradient `DisplacementGradient`; the initial states are zero;
- `FiniteStrainSupport` uses F = I + H;
- the options of the behaviour are now given to the initialize functions
  and post-processings (they were not: the material data managers even
  gave them a null `K` array);
- C, Fortran, Python and Julia bindings.

mfem-mgis:

- parameter `DisplacementGradient` of the behaviour integrators;
- the macroscopic gradients are the macroscopic values of the gradients
  of the behaviour: H, which the caller may keep without rounding;
- the mechanical post-processings use F = I + H;
- refused for the small strain behaviours and the regularisations.

mm-opera-hpc:

- option `--displacement-gradient` of the cermet and of the polycrystal;
- the macroscopic unknowns are the displacement gradient in both modes.

## Results

- Cube of the uniaxial test, Hencky behaviour in float: Newton converges
  down to 1e-8 with H, against a floor of 2e-5 with F
  (`results/displacement-gradient-cube-2026-10-14`); a test of mfem-mgis
  checks the convergence with a relative tolerance of 1e-5.
- Cermet and polycrystal (`results/displacement-gradient-opera-2026-10-14`):
  the quantization of the macroscopic loading and of the fixed point
  disappears; the gap to the double is divided by 3 for the cermet with
  two steps and by 2.4 for the polycrystal at 0.4 s. The remaining error
  comes from the grains, written with the deformation gradient.

## Next steps

- Behaviours of the grains written with H or in a strain framework.
- Plane stress with the deformation gradient: align MGIS and mfem-mgis
  (axial component 1) with TFEL (0).
