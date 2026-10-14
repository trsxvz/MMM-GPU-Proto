---
title: The generic behaviour interface
author: Thomas Helfer
date: 12/09/2018
lang: en-EN
link-citations: true
colorlinks: true
geometry:
  - margin=2cm
papersize: a4
figPrefixTemplate: "$$i$$"
tblPrefixTemplate: "$$i$$"
secPrefixTemplate: "$$i$$"
eqnPrefixTemplate: "($$i$$)"
---

The generic behaviour interface has been introduced to provide an
interface suitable for most needs. Whereas other interfaces target a
specific solver and thus are restricted by choices made by this specific
solver, this interface has been created for developers of homebrew
solvers who are able to modify their code to take full advantage of
`MFront` behaviours.

This interface is tightly linked with the
`MFrontGenericInterfaceSupport` project which is available on `github`:
<https://github.com/thelfer/MFrontGenericInterfaceSupport>. This project
has a more liberal licence which allows it to be included in both
commercial and open-source solvers/library. This licensing choice
explains why this project is not part of the `TFEL`.

The [`generic-parallel` interface](generic-parallel-behaviour-interface.html)
integrates a behaviour over a batch of integration points, on CPUs or
GPUs, with the same conventions.

# Features

## Supported behaviours

The generic behaviour interface supports:

- strain based behaviours. The definition of the strain depends on the
  strain measure used (see the `MFront`'`@StrainMeasure` keyword).
- finite strain behaviours.
- cohesive zone models
- general behaviours which relate a set of gradients to a set of
  thermodynamic forces.

## Supported modelling hypotheses

- Axisymmetrical generalised plane strain.
- AXISYMMETRICAL,
- PLANESTRESS,
- PLANESTRAIN,
- GENERALISEDPLANESTRAIN,
- TRIDIMENSIONAL,

## Finite strain behaviours {#sec:generic_behaviour_interface:finite_strain_behaviours}

The finite strain behaviours are a special case. For greater
flexibility, the generic behaviour interface allows choosing at runtime
the stress measure and the tangent operator to be returned.

The following stress measures are available:

- the Cauchy stress
- the first Piola-Kirchhoff stress
- the second Piola-Kirchhoff stress

The following tangent operators are available:

- the derivative of the Cauchy stress with respect to the deformation
  gradient.
- the derivative of the first Piola-Kirchhoff stress with respect to the
  deformation gradient.
- the derivative of the second Piola-Kirchhoff stress with respect to
  the Green-Lagrange strain.
- the derivative of the Kirchhoff stress with respect to the spatial
  increment of the deformation gradient.

### Displacement gradient {#sec:generic_behaviour_interface:displacement_gradient}

By default, the gradient of a finite strain behaviour is the deformation
gradient $\tenseur{F}$. The calling solver may give the displacement
gradient $\tenseur{H}=\tenseur{F}-\tenseur{I}$ instead. The
deformation gradient being close to the identity, its components are
only known within the precision of the numerical type near $1$: about
$6\,10^{-8}$ in single precision, i.e. $10^{-3}$ in relative value
for a strain of $5\,10^{-5}$. The displacement gradient keeps the
relative precision of small strains:

- the behaviours written in the Green-Lagrange strain framework receive
  $\tenseur{\varepsilon}^{\mathrm{GL}}=\left(\tenseur{H}+\tenseur{H}^{T}+\tenseur{H}^{T}\,\cdot\,\tenseur{H}\right)/2$,
  computed without the cancellation of
  $\left(\tenseur{F}^{T}\,\cdot\,\tenseur{F}-\tenseur{I}\right)/2$
  (see the `computeGreenLagrangeTensorFromDisplacementGradient`
  function);
- the behaviours written in the Hencky strain framework receive a
  logarithmic strain computed from the eigenvalues of
  $\tenseur{C}-\tenseur{I}=2\,\tenseur{\varepsilon}^{\mathrm{GL}}$ (the
  `LogarithmicStrainHandler` class can be built from the deformation
  gradient and the Green-Lagrange strain). The divided differences of
  the logarithm used by the conversions of the stresses and of the
  tangent operators are also computed from these eigenvalues, without
  cancellation, even when two of them are equal;
- the behaviours written with the deformation gradient receive
  $\tenseur{F}=\tenseur{I}+\tenseur{H}$, computed by the interface.

The tangent operators are unchanged: the derivatives with respect to the
deformation gradient are the derivatives with respect to the
displacement gradient.

In plane stress, the axial component of the displacement gradient given
by the caller is zero: the axial component of the deformation gradient
is computed by the behaviour. The interface adds it to the deformation
gradient for the behaviours written with the deformation gradient and
in the Green-Lagrange strain framework, and sets it in the
`LogarithmicStrainHandler` class in the Hencky strain framework.

This option is given at runtime: the stress measure, second value of
the `K` array, is offset by `10` (constant
`mfront::gb::displacement_gradient_stress_measure_offset`), as the type
of computation, first value of the `K` array, is offset by `100` to
request the speed of sound. A library generated before this option
reports an invalid stress measure. The behaviours which accept the
displacement gradient export the symbol
`<behaviour_function_name>_AcceptsDisplacementGradient`, an
`unsigned short` equal to `1`. This holds for the `generic` and
`generic-parallel` interfaces.

## Support of orthotropic behaviours

Orthotropic behaviours require to:

- rotate the gradients in the material frame to the global frame before
  the behaviour integration.
- rotate the thermodynamic forces and the tangent operator blocks from
  the material frame to the global frame.

By design, the `generic` behaviour interface does not automatically
perform those rotations as part of the behaviour integration but
generates additional functions to do it. This choice allows the calling
solver to use their own internal routines to handle the rotations
between the global and material frames.

However, the `generic` interface also generates helper functions which
can perform those rotations. Those functions are named as follows:

- `<behaviour_function_name>_<hypothesis>_rotateGradients`
- `<behaviour_function_name>_<hypothesis>_rotateThermodynamicForces`
- `<behaviour_function_name>_<hypothesis>_rotateTangentOperatorBlocks`

They all take three arguments: 

- a pointer to the location where the rotated variables will be stored.
- a pointer to the location where the original variables are stored.
- a pointer to the rotation matrix from the global frame to the material
  frame. The rotation matrix has 9 components stored in column-major
  format. For the function handling the thermodynamic forces and the
  tangent operators blocks, this rotation matrix is transposed
  internally to get the rotation matrix from the material frame to the
  global frame.

In place rotations are explicitly allowed, i.e. the first and second
arguments can be a pointer to the same location.

The rotated variables can be scalars (which are copied), vectors,
symmetric tensors and tensors, as well as arrays of such variables. All
the tangent operator blocks can be rotated. A block associated with
arrays of variables is rotated entry by entry. A block for which `TFEL`
has no type which can be rotated (derivative of a vector with respect to
a symmetric tensor or a tensor, and conversely, derivative of a tensor
with respect to a symmetric tensor) is rotated by
`mfront::gb::rotateTangentOperatorBlock`, which computes
\(Q_{1}\,K\,Q_{2}^{T}\), where \(Q_{i}\) is the matrix of the rotation
of the objects of each variable. The generated code can be used in device
code.

The three previous functions work for an integration point. Three other
functions are also generated:

- `<behaviour_function_name>_<hypothesis>_rotateArrayOfGradients`
- `<behaviour_function_name>_<hypothesis>_rotateArrayOfThermodynamicForces`
- `<behaviour_function_name>_<hypothesis>_rotateArrayOfTangentOperatorBlocks`

Those functions take an additional argument which is the number of
integration points to be treated.

### Finite strain behaviours

As devised in Section
@sec:generic_behaviour_interface:finite_strain_behaviours, finite strain
behaviours are a special case, because the returned stress measure and
the returned tangent operator can be chosen at run time. A specific
rotation function is generated for each supported stress measure and
each supported tangent operator.

For those behaviours, the `rotateThermodynamicForces` and
`rotateTangentOperatorBlocks` functions are not generated. Here is the
list of the generated functions:

-
  `<behaviour_function_name>_<hypothesis>_rotateThermodynamicForces_CauchyStress`.
  This function assumes that its first argument is the Cauchy stress in
  the material frame.
-
  `<behaviour_function_name>_<hypothesis>_rotateThermodynamicForces_PK1Stress`.
  This function assumes that its first argument is the first
  Piola-Kirchhoff stress in the material frame.
-
  `<behaviour_function_name>_<hypothesis>_rotateThermodynamicForces_PK2Stress`.
  This function assumes that its first argument is the second
  Piola-Kirchhoff stress in the material frame.
-
  `<behaviour_function_name>_<hypothesis>_rotateTangentOperatorBlocks_dsig_dF`.
  This function assumes that its first argument is the derivative of the
  Cauchy stress with respect to the deformation gradient in the material
  frame.
-
  `<behaviour_function_name>_<hypothesis>_rotateTangentOperatorBlocks_dPK1_dF`.
  This function assumes that its first argument is the derivative of the
  first Piola-Kirchhoff stress with respect to the deformation gradient
  in the material frame.
-
  `<behaviour_function_name>_<hypothesis>_rotateTangentOperatorBlocks_dPK2_dEGL`.
  This function assumes that its first argument is the derivative of the
  second Piola-Kirchhoff stress with respect to the Green-Lagrange
  strain in the material frame.
-
  `<behaviour_function_name>_<hypothesis>_rotateTangentOperatorBlocks_dtau_ddF`.
  This function assumes that its first argument is the derivative of the
  Kirchhoff stress with respect to the spatial increment of the
  deformation gradient in the material frame.

Each of those functions has a variant treating an array of integration
points, whose name is obtained by replacing `rotate` by `rotateArrayOf`
(for instance
`<behaviour_function_name>_<hypothesis>_rotateArrayOfThermodynamicForces_CauchyStress`).

## Numerical type

By default, the data exchanged with the generated functions are of type
`double`. The `real_type` DSL option selects `float` instead:

~~~~{.bash}
$ mfront --obuild --interface=generic --dsl-option=real_type:float \
    Plasticity.mfront
~~~~

The type of the data (`mfront_gb_real`, declared in
`MFront/GenericBehaviour/Types.h`) is then given by the
`MFRONT_GB_REAL_TYPE` macro, defined by the generated sources. Each
behaviour exports the size of its numerical type
(`<behaviour_function_name>_real_size`), and the version of its
application programming interface is shifted by 1000 for another type
than `double`: a caller expecting `double`, such as `MTest`, refuses the
library. With `float`, the convergence criteria (`@Epsilon`) are raised to
a tenth of the machine epsilon of `float` (about `1.2e-8`) and the
perturbation of the numerical jacobian to the machine epsilon (about
`1.2e-7`), with a warning. The same option is supported by the
[`generic-parallel` interface](generic-parallel-behaviour-interface.html)
and by the `generic` and `generic-parallel` interfaces for material
properties (`MFRONT_GMP_REAL_TYPE` macro). See the release notes of
version 5.2 for the writing of behaviours usable with `float`.
