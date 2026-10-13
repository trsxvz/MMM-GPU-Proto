.. _mfem_mgis_boundary_conditions:

===================
Boundary conditions
===================

.. contents::
    :depth: 3
    :local:

General purpose boundary conditions
===================================

The :code:`UniformDirichletBoundaryCondition` boundary condition
----------------------------------------------------------------

The :code:`UniformDirichletBoundaryCondition` allows imposing
a prescribed value on a boundary.

The following snippet sets to zero the first component of the unknowns
on the boundary labeled :code:`5`:

.. code:: c++

   problem.addBoundaryCondition(
       ctx, std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
                problem.getFiniteElementDiscretizationPointer(), 5, 0)) |
       or_die;

The boundary may also be given by a name (a regular expression selecting
the boundaries by name), and a last argument may give the imposed value
as a function of time.

The :code:`addUniformDirichletBoundaryCondition` method of the
:code:`NonLinearEvolutionProblem` class builds this boundary condition
from parameters:

- :code:`Boundary` or :code:`Boundaries`: the boundaries, given by their
  identifiers or by their names;
- :code:`Component`: the component of the unknowns;
- :code:`LoadingEvolution` (optional): the imposed value as a function of
  time (zero by default).

.. code:: c++

   problem.addUniformDirichletBoundaryCondition(
       ctx, {{"Boundary", "UpperBoundary"},
             {"Component", 1},
             {"LoadingEvolution", [](const auto t) { return 6e-3 * t; }}}) |
       or_die;

The same parameters can be given to the constructor
:code:`UniformDirichletBoundaryCondition(problem, parameters)`.

The :code:`ImposedDirichletBoundaryConditionAtClosestNode` boundary condition
-----------------------------------------------------------------------------

The :code:`ImposedDirichletBoundaryConditionAtClosestNode` boundary
condition imposes the value of a component of the unknowns at the node
closest to a given point (given by its 2 or 3 coordinates), for instance
to block the rigid body motions:

.. code:: c++

   problem.addBoundaryCondition(
       ctx,
       std::make_unique<mfem_mgis::ImposedDirichletBoundaryConditionAtClosestNode>(
           problem.getFiniteElementDiscretizationPointer(),
           std::array<mfem_mgis::real, 3u>{0, 0, 0}, 0)) |
       or_die;

A last argument may give the imposed value as a function of time.

Mechanical boundary conditions
==============================

The :code:`UniformImposedPressure` boundary condition
-----------------------------------------------------

The :code:`UniformImposedPressure` boundary condition computes the
residual associated with the following variational operator:

.. math::

   -\int_{\partial\Omega_{r}} p\,\vec{N}\,\cdot\,\,\vec{u}^{\star}\,\mathrm{d}S

where:

- :math:`\partial\Omega_{r}` denotes the boundary in the reference
  configuration on which the pressure is applied,
- :math:`p` is the imposed pressure,
- :math:`\vec{N}` is normal to the boundary in the reference configuration,
- :math:`\vec{u}^{\star}` is a virtual displacement.

The following snippet imposes a uniform pressure, which evolves linearly
with time, on the boundary labeled :code:`3`:

.. code:: c++
   
   problem.addBoundaryCondition(
       ctx,
       std::make_unique<mfem_mgis::UniformImposedPressureBoundaryCondition>(
           problem.getFiniteElementDiscretizationPointer(), 3,
           [](const mfem_mgis::real t) noexcept { return 150e6 * t; })) |
       or_die;

The boundary may also be given by its name. The constructor
:code:`UniformImposedPressureBoundaryCondition(problem, parameters)` takes
the parameters :code:`Boundary` or :code:`Boundaries`, and
:code:`LoadingEvolution` (required), the pressure as a function of time.

Heat transfer boundary conditions
=================================

The :code:`UniformHeatSource` boundary condition
------------------------------------------------

The :code:`UniformHeatSource` boundary condition is applied on materials,
not on boundaries. It computes the
residual associated with the following variational operator:

.. math::

   \int_{\Omega_{r}} q\,T^{\star}\,\mathrm{d}V

where:

- :math:`\Omega_{r}` denotes the domain of interest in the reference
  configuration,
- :math:`q` is the heat source,
- :math:`T^{\star}` is a virtual temperature.

The following snippet imposes a uniform heat source, which evolves
linearly with time, on the material named :code:`fuel`:

.. code:: c++
   
   problem.addBoundaryCondition(
       ctx, std::make_unique<mfem_mgis::UniformHeatSourceBoundaryCondition>(
                problem.getFiniteElementDiscretizationPointer(), "fuel",
                [&q](const auto t) { return q * t; })) |
       or_die;

The material may also be given by its identifier. The constructor
:code:`UniformHeatSourceBoundaryCondition(problem, parameters)` takes the
parameters :code:`Material` or :code:`Materials`, and
:code:`LoadingEvolution` (required), the heat source as a function of
time.

Boundary conditions independent of the unknowns
===============================================

The :code:`UniformImposedPressure` and :code:`UniformHeatSource` boundary
conditions do not depend on the unknowns. With the inner
forces computed by kernels (:code:`InnerForcesComputation` set to
:code:`Kernels`), they are not integrators of the non linear form: their
contributions to the residual are assembled once, for a unit value, and
added to the inner forces, on the device if a device is enabled. If the
non linear form contains no other integrator than the one of
:code:`mfem-mgis` (no boundary condition depending on the unknowns), the
tangent operator can then be applied without assembly
(:code:`TangentOperatorComputation` set to :code:`MatrixFree`). With the
inner forces computed element by element, they are integrators of the non
linear form, and the tangent operator can't be applied without assembly.
In both cases, the tangent operator can be assembled by kernels
(:code:`TangentOperatorComputation` set to :code:`Kernels`), the
contributions of the other integrators of the non linear form being added
to the matrix assembled by kernels (see :ref:`mfem_mgis_batch_computations`).
