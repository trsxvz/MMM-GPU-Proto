.. _mfem_mgis_release_notes_1_1:

============================
Release notes of Version 1.1
============================

This version inherits from all the features introduced in:

- Version 1.0.1
- Version 1.0.2
- Version 1.0.3
- :ref:`mfem_mgis_release_notes_1_0_4`

.. contents:: Table of Contents
   :depth: 3
   :local:
   :backlinks: none

Highlights
==========

- `NonLinearEvolutionProblem` are now able to make a prediction of the
  solution, removing major convergence issues when imposed displacements
  are imposed.
- The :math:`\bar{F}` method, or FBar formulation, has been implemented
  for finite strain behaviours to handle nearly incompressible materials.
- The regularization proposed by Faltus et al. in the context of the
  third medium contact has been implemented for plane strain, plane
  stress and tridimensional hypotheses.
- `MGIS`'s contexts now handle gathering computation time information to
  create the performance table instead of the previously used
  `CatchTimeSection`.
- Many methods have been deprecated to have a consistent error handling
  scheme based on `MGIS`'s one. As such, many methods and functions now
  take a `MGIS`'s :cxx:`Context` as their first argument. The deprecated
  methods of :cxx:`AbstractNonLinearEvolutionProblem` have been removed,
  see the known incompatibilities below.
- The name of the materials and boundaries are automatically retrieved
  from |MFEM|'s mesh.

Known incompatibilities
=======================

- In previous versions, the failure of Hypre's linear iterative solvers
  were discarded by `MFEM/MGIS`'s Newton solver due to the lack of
  methods to test their convergence in |MFEM|'s version prior to 4.10.
  The parameter `DiscardLinearSolverFailure` can be passed to
  `MFEM/MGIS`'s Newton solver to recover the behavior of previous
  versions.
- The deprecated methods of :cxx:`AbstractNonLinearEvolutionProblem`
  and of its derived classes which did not take a :cxx:`Context` have
  been removed, as well as the deprecated functions :cxx:`getGradient`,
  :cxx:`getThermodynamicForce` and :cxx:`getInternalStateVariable`
  which did not take a :cxx:`Context`. The overloads taking a
  :cxx:`Context` as first argument shall be used instead.
- Partial quadrature spaces are now associated with a location, i.e. a
  material or a boundary (see Issue 466):

  - the :cxx:`getId` method of :cxx:`PartialQuadratureSpace` has been
    removed, the :cxx:`getLocation` method shall be used instead.
  - the constructors of :cxx:`PartialQuadratureSpace` taking a material
    identifier have been removed. A :cxx:`LocationIdentifier`, built
    from a :cxx:`MaterialIdentifier` or a :cxx:`BoundaryIdentifier`,
    shall be given instead, e.g. :cxx:`MaterialIdentifier{.id = m}`.
  - a partial quadrature space defined on a boundary using a finite
    element discretization is now built on the boundary elements of the
    mesh rather than on the submesh of the boundary.
  - :cxx:`LocationIdentifier` is no longer an aggregate and can not be
    assigned anymore.
  - :cxx:`QPDependency::getMaterialIdentifier`,
    :cxx:`PartialQuadratureFunctionsSet::getMaterialIdentifiers` and
    the overload of :cxx:`PartialQuadratureFunctionsSet::get` taking a
    material identifier have been replaced by :cxx:`getLocation`,
    :cxx:`getLocations` and an overload taking a
    :cxx:`LocationIdentifier`.
  - the constructor of :cxx:`QPDependency` taking a material identifier
    now takes a :cxx:`LocationIdentifier`.
  - the static :cxx:`buildQuadratureSpace` methods of the behaviour
    integrators now take :cxx:`attributes::Throwing` as first argument.
  - :cxx:`PartialQuadratureSpace::treatInvalidElementIndex` now takes a
    :cxx:`LocationIdentifier` instead of a material identifier.
  - :cxx:`PartialQuadratureSpace::getLocationName` now returns a
    description of the location, e.g. ``material 'BEAM' (12)``. It
    previously returned the name alone for named materials and
    boundaries. The output of :cxx:`getInformation` has changed
    accordingly.
  - the overloads of :cxx:`updateGridFunction` taking a mesh or a
    submesh have been removed. :cxx:`updateGridFunction` requires a grid
    function whose finite element space is managed by the finite element
    spaces manager of the partial quadrature functions, as the ones
    created by :cxx:`makeGridFunction`.
  - :cxx:`updateL2Projection` now reports an error if a function is
    defined on a material or a boundary which is not part of the mesh of
    the result. Such functions were previously ignored.
  - the :cxx:`setMaterialProperty` and :cxx:`setExternalStateVariable`
    methods of behaviour integrators, called with an evaluator, require
    an evaluator defined on an equivalent partial quadrature space, i.e.
    built on a finite element discretization sharing the same finite
    element spaces manager. Those calls are now collective.
- :cxx:`updateGridFunction` now takes a :cxx:`Context` as first
  argument and returns a boolean (see Issue 475).
- :cxx:`assign_values` now reports an error if both views are not
  defined on the same mesh discretization, even if both discretizations
  are built from the same mesh (see Issue 476). Two problems
  exchanging values must be built on the same mesh discretization, for
  instance by sharing a finite element spaces manager.
- the deprecated functions :cxx:`getMPIrank` and :cxx:`getMPIsize` have
  been removed. :cxx:`getMPIRankForDefaultCommunicator` and
  :cxx:`getMPISizeForDefaultCommunicator` shall be used instead (see
  Issue 474).
- The location where a linear solver computes is given by its
  description (see :ref:`linear_solvers`), and no longer deduced from its
  name: a linear solver registered by the
  :cxx:`LinearSolverFactory::add` method without description is assumed
  to compute where the device of |MFEM| computes, whatever its name.
- The preconditioners are registered in the factory of preconditioners:
  the error message of an unknown preconditioner lists the available
  ones.

New features
============

Descriptions of the linear solvers and factory of preconditioners
------------------------------------------------------------------

Each linear solver and each preconditioner is registered with a
description (:cxx:`LinearSolverDescription`): library, kind (iterative,
direct, approximate solver or preconditioner), location of the
computations, accepted preconditioners and options. The preconditioners
are registered in the new factory of preconditioners
(:cxx:`PreconditionerFactory`). The location of each linear solver and
preconditioner is checked against the execution settings, and the
preconditioners which a linear solver does not accept are refused with a
message. The available linear solvers and preconditioners are listed in
:ref:`linear_solvers`; the solvers of each library are registered by a
module (see :ref:`linear_solver_modules`).

Linear solver and preconditioner overrides
------------------------------------------

The options :code:`--linear-solver-override` and
:code:`--preconditioner-override` (declared by :cxx:`declareDefaultOptions`)
or the environment variables :code:`MFEM_MGIS_LINEAR_SOLVER_OVERRIDE` and
:code:`MFEM_MGIS_PRECONDITIONER_OVERRIDE` replace the linear solvers and
the preconditioners set by the program by the ones described by their names
and their options, in the syntax of the data of TFEL, for the comparisons
of linear solvers and preconditioners on existing programs. See
:ref:`linear_solver_override`.

One option per optional library of MFEM
---------------------------------------

Each library of |MFEM| which mfem-mgis uses for its solvers and
preconditioners (PETSc, Ginkgo, MUMPS, SuperLU_DIST, STRUMPACK,
SuiteSparse, SUNDIALS, the sparse direct solvers of oneMKL, cuDSS and
MAGMA) has its own CMake option, :code:`enable-<library>`: :code:`AUTO` by
default (the library is used if |MFEM| is built with it), :code:`ON` (the
library is required, the configuration fails if |MFEM| is not built with
it) or :code:`OFF` (the library is not used, even if |MFEM| is built with
it). The result is given by the macros :code:`MFEM_MGIS_HAVE_<LIBRARY>`
and by the variables :code:`MFEMMGIS_HAVE_<LIBRARY>` of the installed
CMake configuration. A package manager can map each library to one switch,
for instance a variant of a Spack package. See
:ref:`mfem_mgis_optional_libraries`.

Linear solvers and preconditioners of PETSc
-------------------------------------------

The :cxx:`PETScLinearSolver` linear solver (a Krylov solver of PETSc)
and the :cxx:`PETScPreconditioner` preconditioner (a preconditioner of
PETSc) are available when |MFEM| is built with PETSc. They are
configured by the options of PETSc given by their :cxx:`Options`
parameter, and give access to the Krylov solvers and to most
preconditioners of PETSc, among which its algebraic multigrid (GAMG). A preconditioner of
mfem-mgis can be given to :cxx:`PETScLinearSolver`, and
:cxx:`PETScPreconditioner` can be used by the Krylov solvers of |MFEM|.
For the linear elasticity, with the unknowns ordered by vector dimension,
the rigid body modes are given to the preconditioners of PETSc (GAMG).
The solvers of PETSc compute on the host, the matrices of hypre being
copied without initializing the GPU. The :cxx:`Fields` option set to
:cxx:`Components` gives the field split preconditioner of PETSc one field
per component of the unknowns, whatever their ordering. See
:ref:`petsc_linear_solvers`.

Linear solvers and preconditioners of Trilinos
----------------------------------------------

The Krylov solvers of Belos (:cxx:`BelosSolver`), the direct solvers of
Amesos2 (:cxx:`Amesos2Solver`, :cxx:`Amesos2Preconditioner`), the
preconditioners of Ifpack2 (:cxx:`Ifpack2`), the smoothed aggregation
algebraic multigrid of MueLu (:cxx:`MueLu`, with the rigid body modes as
near null space in mechanics) and the two level Schwarz preconditioner of
FROSch (:cxx:`FROSch`) are available in parallel when mfem-mgis is
built with Trilinos (CMake option :code:`enable-trilinos`), in the
numerical type of mfem-mgis (:code:`float` requiring Tpetra instantiated
for it, a combination which is not tested). They compute on the GPU when
Trilinos is built with CUDA (tested) or HIP (not tested). When Trilinos
is built with them, the fast incomplete factorizations of ShyLU
(:cxx:`FAST_ILU`, :cxx:`FAST_ILU_B`, :cxx:`FAST_IC`, :cxx:`FAST_ILDL`
types of :cxx:`Ifpack2`) and the block preconditioners of Teko
(:cxx:`Teko`, one block per component of the unknowns by default, CMake
option :code:`enable-teko`) are also available. See
:ref:`trilinos_linear_solvers`.

Linear solvers and preconditioners of Ginkgo
--------------------------------------------

The Krylov solvers of Ginkgo (CG, FCG, pipelined CG, BiCG, BiCGSTAB, CGS,
GMRES, CB-GMRES, GCR, IDR, MINRES), its iterative refinement, its sparse
direct solver and its preconditioners (Jacobi, ILU, IC and their variants
with incomplete sparse approximate inverses, SOR, Gauss-Seidel, algebraic
multigrid) are available when |MFEM| is built with Ginkgo, on the CPU and
on the GPU (they compute where the device of |MFEM| computes), in
sequential and, for the distributed ones, in parallel, the local
preconditioners being wrapped in the Schwarz preconditioner of Ginkgo.
Each entry is described by a configuration of Ginkgo, completed by its
:cxx:`Configuration` option, and the generic :cxx:`GinkgoSolver` and
:cxx:`GinkgoPreconditioner` entries accept any configuration of Ginkgo.
See :ref:`ginkgo_linear_solvers`.

Solvers and preconditioners of hypre
------------------------------------

The BiCGSTAB, COGMRES and LGMRES Krylov solvers of hypre, its hybrid
solver (diagonal scaling, then BoomerAMG when the convergence is slow),
its factorized sparse approximate inverse (:cxx:`HypreFSAI`), its
multigrid reduction (:cxx:`HypreMGR`, the unknowns being split by their
components, see :ref:`hypre_mgr`) and its smoothers
(:cxx:`HypreSmoother`: Jacobi, l1-Jacobi, Gauss-Seidel, polynomial
smoothers) are available. The BoomerAMG preconditioner has new
options: truncation of the interpolation, aggressive coarsening, size of
the coarsest level, relaxation weight, order of the Chebyshev smoothers
and transposed interpolations kept. The values of the
:cxx:`StrengthThreshold` option of BoomerAMG out of :math:`[0, 1]` and the
Krylov dimensions of :cxx:`HypreFGMRES` lower than one are now refused, as
well as the :cxx:`Elasticity` strategy of BoomerAMG for unknowns whose
number of components is not the dimension of the mesh, or a mesh whose
dimension is not the space dimension (:code:`MFEM` aborted or computed
wrong rigid body modes, or the strategy was replaced by :cxx:`System`
with unknowns ordered by nodes or hypre on a GPU).
See :ref:`linear_solvers`.

Sparse direct solvers of oneMKL
-------------------------------

The PARDISO solver of oneMKL (:cxx:`MKLPardisoSolver`, sequential) and
its Cluster Sparse Solver (:cxx:`MKLClusterSparseSolver`, parallel) are
available when :code:`MFEM` is built with them. See
:ref:`mkl_linear_solvers`.

Sparse solvers of MAGMA
-----------------------

The Krylov solvers of the sparse part of MAGMA (sequential) and its
Jacobi, incomplete LU and incomplete Cholesky preconditioners (applied to
the diagonal block of each process in parallel), which compute on one
GPU, are available when :code:`MFEM` is built with MAGMA in double
precision and the library :code:`magma_sparse` of MAGMA is found. See
:ref:`magma_linear_solvers`.

Jacobi preconditioner by nodal blocks
-------------------------------------

The :cxx:`NodalBlockJacobi` preconditioner inverts the blocks of the
matrix which couple the components of the unknown at each node, with the
batched dense linear algebra of :code:`MFEM` (native kernels, cuBLAS or
hipBLAS, MAGMA). With the tangent operator applied without assembly, the
blocks are computed by kernels. See :ref:`nodal_block_jacobi`.

Nonlinear solver of KINSOL
--------------------------

The :cxx:`KINSOL` nonlinear solver uses KINSOL (:code:`SUNDIALS`) when
:code:`MFEM` is built with it: Newton's method with or without line
search, modified Newton method, Picard's method with Anderson's
acceleration. See :doc:`../nonlinear_algorithms`.

Nonlinear solvers of PETSc
--------------------------

The :cxx:`SNES` nonlinear solver uses the nonlinear solvers of PETSc
(SNES) when :code:`MFEM` is built with it, in parallel: Newton's method
with or without line search, trust region, quasi-Newton and the other
methods of SNES, selected by the options of PETSc, with the linear solver
of the problem or the Krylov solvers and preconditioners of PETSc. Unlike
the global mode of :cxx:`setPETSc`, it is selected per problem. See
:doc:`../nonlinear_algorithms`.

Preconditioners of the sparse matrices of MFEM, options of MUMPS
----------------------------------------------------------------

The Gauss-Seidel and Jacobi smoothers of the sparse matrices of MFEM
(:cxx:`GSSmoother`, :cxx:`DSmoother`), its block incomplete LU
factorization (:cxx:`BlockILU`) and the direct solvers of SuiteSparse
(:cxx:`UMFPack`, :cxx:`KLU`) are available as preconditioners, applied to
the diagonal block of each process in parallel, as well as the
:cxx:`KLUSolver` direct solver in sequential and the orthogonalizing
wrapper :cxx:`OrthoSolver` for the operators whose kernel is the constant
vector. The :cxx:`MUMPSSolver` direct solver accepts a reordering
strategy, a tolerance of its block low-rank factorization and a verbosity
level; its factorization is also available as a preconditioner
(:cxx:`MUMPSPreconditioner`). See :ref:`linear_solvers`.

Direct solvers of SuperLU_DIST and STRUMPACK
--------------------------------------------

The direct solver of SuperLU_DIST (:cxx:`SuperLUSolver`) and the direct
solvers of STRUMPACK, in double and in mixed precision, with optional
compression of the frontal matrices (:cxx:`STRUMPACKSolver`,
:cxx:`STRUMPACKMixedPrecisionSolver`, :cxx:`STRUMPACKPreconditioner`),
are available when |MFEM| is built with these libraries. See
:ref:`strumpack_linear_solvers`.

Prediction of the solution
--------------------------

By default, a nonlinear evolution problem uses the solution at the
beginning of the time step, modified by applying Dirichlet boundary
conditions, as the initial guess of the solution at the end of the time
step, see below for details.

This can be changed by using the :cxx:`setPredictionPolicy`
method, as follows:

.. code-block:: c++

   // use the elastic operator by default
   mechanics.setPredictionPolicy(
      {.strategy = PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});

Available strategies
^^^^^^^^^^^^^^^^^^^^

Default prediction (:cxx:`PredictionStrategy::DEFAULT_PREDICTION`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

By default, a nonlinear evolution problem uses the solution at the
beginning of the time step, modified by applying Dirichlet boundary
conditions, as the initial guess of the solution at the end of the time
step.

.. warning::

   In mechanics, this may lead to very high increments of the deformation
   gradients or the strain in the neighboring elements of boundaries where
   evolving displacements are imposed.
   
Prediction for the state at the beginning of the time step (:cxx:`PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The prediction for the state at the beginning of the time step strategy
determines the increment of the displacement :math:`\Delta\,\mathbb{u}`
by solving the following linear system:

.. math::

   \mathbb{K}\,\cdot\,\Delta\,\mathbb{u} = \ets{\mathbb{F}_{e}}-\bts{\mathbb{F}_{i}}

where:

- :math:`\mathbb{K}` denotes the prediction operator (see below).
- :math:`\ets{\mathbb{F}_{e}}` denotes the external forces at the end of
  the time step.
- :math:`\bts{\mathbb{F}_{i}}` denotes the inner forces at the beginning
  of the time step.
- :math:`\Delta\,\mathbb{u}` is submitted to the increment of the
  imposed Dirichlet boundary conditions.

.. note::

   Although the wording explicitly refers to mechanics, this equation
   applies to all physics.

The following prediction operators can be chosen:

- :cxx:`PredictionOperator::ELASTIC`: the elastic operator
- :cxx:`PredictionOperator::SECANT`: the secant operator is typically
  defined by the elastic operator
- :cxx:`PredictionOperator::TANGENT`: the tangent operator,
  defined by the time-continuous derivative of the thermodynamic force
  with respect to the gradients.
- :cxx:`PredictionOperator::LAST_ITERATE_OPERATOR`: this operator shall
  reuse the one computed at the last iteration of the previous time step.
  It is not implemented yet: the elastic operator is always used.

Prediction based on a behaviour integration with constant gradients (:cxx:`PredictionStrategy::CONSTANT_GRADIENTS_INTEGRATION_PREDICTION`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :cxx:`CONSTANT_GRADIENTS_INTEGRATION` strategy determines the
increment of the unknown :math:`\Delta\,{u}` by solving the following
linear system:
     
.. math::

   \tilde{\mathbb{K}}\,\cdot\,\Delta\,\mathbb{u} =
   \ets{\mathbb{F}_{e}}-\ets{\tilde{\mathbb{F}}_{i}}
     
where:
     
- :math:`\tilde{\mathbb{K}}` denotes the operator computed at the end of
  the behaviour integration.
- :math:`\ets{\mathbb{F}_{e}}` denotes the external forces at the end of
  the time step.
- :math:`\ets{\tilde{\mathbb{F}}_{i}}` denotes an approximation inner
  forces at the end of the time step computed by assuming that the
  gradients are constant over the time (and thus equal to their values at
  the beginning of the time step).
- :math:`\Delta\,\mathbb{u}` is submitted to the increment of the
  imposed Dirichlet boundary conditions.

The behaviour integration allows taking into account:

- the evolution of stress-free strain over the time step (thermal
  expansion, swelling, etc..),
- the viscoplastic relaxation of the stress. This relaxation can be
  discarded by integrating the behaviour with a null time step.

The following operators are available:

- :cxx:`IntegrationOperator::ELASTIC`: the elastic operator,
- :cxx:`IntegrationOperator::SECANT`: the secant operator is typically
  defined by the elastic operator affected by damage,
- :cxx:`IntegrationOperator::TANGENT`: the tangent operator, defined by
  the time-continuous derivative of the thermodynamic force with respect
  to the gradients,
- :cxx:`IntegrationOperator::CONSISTENT_TANGENT`: the consistent
  tangent operator, defined by the derivative of the thermodynamic force
  with respect to the gradients at the end of the time step. See
  :cite:`simo_consistent_1985` for details.

FBar formulation
----------------

The :math:`\bar{F}` method, or FBar formulation, is implemented following
:cite:`de_souza_neto_design_1996`. This formulation is designed to handle
nearly incompressible materials in large strain analysis by using a modified
deformation gradient :math:`\bar{\underline{F}}` that separates volumetric and
deviatoric responses.

The method replaces the standard deformation gradient :math:`\underline{F}`
with an assumed modified counterpart :math:`\bar{\underline{F}}` in the
computation of stresses. This modification is based on a multiplicative split
into volumetric and deviatoric parts:

.. math::

   \underline{F} = J^{1/3}\, \bar{\underline{F}}

where :math:`J = \det{\underline{F}}` is the Jacobian of the deformation
gradient. This formulation effectively avoids locking issues in nearly
incompressible materials while maintaining accuracy.

The :math:`\bar{F}` formulation is particularly suited for low-order finite
elements and is applicable to arbitrary material models. It ensures quadratic
rates of convergence in Newton-Raphson schemes.

This formulation is enabled by passing an additional parameter to the
:cxx:`Mechanics` behaviour integrator, as follows:

.. code:: c++

  const auto fbar_parameters = mfem_mgis::dict{
      {"Regularization", mfem_mgis::dict{{"FBar", mfem_mgis::dict{}}}}};
  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "beam", library,
                                   behaviour, fbar_parameters) | or_die;

Faltus 2026 regularization
---------------------------

The regularization proposed by Faltus et al. in the context of contact
mechanics using a third medium is implemented here :cite:`faltus_deformation_2026`. This
regularization only applies to finite strain behaviours. Currently, this
regularization is only available for isotropic behaviours.

This regularization adds a contribution to the standard variational
operator in finite strain and can be derived from an energy :math:`W`
which penalizes the difference between the deformation gradient
:math:`\underline{F}` at a given quadrature point and its value
:math:`\bar{\underline{F}}` at the centroid of the element:

.. math::

    W\left(\underline{F}, \bar{\underline{F}}\right) =
    \alpha\,\left(\underline{F}-\bar{\underline{F}}\right)\,\colon\,
    \left(\underline{F}-\bar{\underline{F}}\right)

where :math:`\alpha` is a penalization coefficient.

This regularization is enabled by passing an additional parameter to the
:cxx:`Mechanics` behaviour integrator, as follows:

.. code:: c++

  const auto faltus_parameters = mfem_mgis::Parameters{
      {"Regularization",
       mfem_mgis::Parameters{
           {"Faltus2026",
            mfem_mgis::Parameters{{"PenalizationCoefficient", 1e11}}}}}};
  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "ThirdMedium", library,
                                   behaviour2, faltus_parameters) | or_die;

The :cxx:`info` function
------------------------

The :cxx:`info` function allows displaying information about an object
in an output stream.


Retrieving information on a finite element discretization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example of usage
""""""""""""""""

.. code:: c++

   const auto& fed = problem.getFiniteElementDiscretization();
   const auto success = mfem_mgis::info(ctx, fed);

Retrieving information on a partial quadrature space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :cxx:`PartialQuadratureSpaceInformation` structure contains some
relevant information about a partial quadrature space:

- the identifier and a description of the underlying material or
  boundary, including its name, if defined,
- the total number of elements,
- the total number of integration points,
- the number of quadrature points per geometric type.

This structure is created by:

- :cxx:`getLocalInformation`, which returns the information relative to
  the current process.
- :cxx:`getInformation`, which returns the information gathered from all
  processes.

The :cxx:`PartialQuadratureSpaceInformation` structure can be printed to
an output stream using the :cxx:`info` function.

Example of usage
""""""""""""""""

.. code:: c++

   auto& bi = problem.getBehaviourIntegrator(ctx, 1, 0) | or_die;
   const auto& qspace = bi.getPartialQuadratureSpace();
   const auto success = mfem_mgis::info(ctx, std::cout, qspace);

Resolution of dependencies of a nonlinear evolution problem using another
-------------------------------------------------------------------------

The :cxx:`resolveBehaviourIntegratorsDependencies` function resolves the
material properties and external state variables of the behaviour
integrators of a nonlinear evolution problem using the gradients,
thermodynamic forces and internal state variables of another one.

Both problems must exchange data defined on equivalent partial
quadrature spaces. Hence, the finite element discretizations of both
problems must share the same finite element spaces manager, i.e. the
finite element discretization of one problem must be built from the
finite element spaces manager of the other.

The :cxx:`Simulation` class
---------------------------

The :cxx:`Simulation` class manages the time steps of a simulation over
several temporal sequences (see :ref:`simulation_class_overview`).

Line-search-like handling of behaviour integration failures
-----------------------------------------------------------

When a new estimate of the Newton solver leads to an integration failure,
the amplitude of the correction is divided by two (see
:ref:`mfem_mgis_nonlinear_algorithms`).

Profiling Toolkit
-----------------

Timers are now accessible by MGIS contexts, offering greater flexibility:

- Eliminates the need for a single global object
- Allows multiple contexts and the ability to distinguish different parts of the code managed through separate objects
- Prevents potential negative interactions with other code using the same timer type.

Example of usage:

.. code:: c++

    CatchTimeSection(ctx, "Class::FunctionName");

Refactoring of submesh and finite element space creation
--------------------------------------------------------

The creation of submeshes and associated finite element spaces has been refactored to be handled through the :cxx:`MeshDiscretization` and :cxx:`FiniteElementSpacesManager` classes. This refactoring introduces the following new classes:

- :cxx:`MeshDiscretization`: A fundamental class that handles the lifetime of meshes and provides utilities for managing mesh-related operations, including submesh creation based on material or boundary identifiers.
- :cxx:`FiniteElementDiscretization`: Extends :cxx:`MeshDiscretization` to handle the lifetime of finite element collections and spaces, providing a high-level interface for creating and managing finite element spaces.
- :cxx:`FiniteElementSpacesManager`: Manages similar finite element spaces (siblings) that share the same mesh and finite element collection but may have different vectorial dimensions, with automatic reuse of existing spaces.

This refactoring improves code organization, reduces duplication, and provides a more consistent and flexible API for working with submeshes and finite element spaces.

Partial quadrature functions defined on another mesh
----------------------------------------------------

A partial quadrature function can be updated from a grid function
defined on another mesh of the same discretization, i.e. the whole
mesh, a submesh of materials or a submesh of boundaries (see Issue
466). The correspondence between the elements of both meshes is given
by :cxx:`MeshDiscretization::getElementsMapping`.

A partial quadrature space defined on a boundary can be built on the
boundary elements of the whole mesh or of a submesh of materials.
Boundary elements oriented differently in a submesh of materials and
in the whole mesh are refused.

The values of a partial quadrature function can be copied by
:cxx:`assign_values` into a function defined on the same location of
another mesh of the same discretization. The values are then copied
element by element, and the call is collective in parallel (see Issue
476).

Search of the behaviour integrator providing a quantity
-------------------------------------------------------

The functions :cxx:`hasGradientProvider`,
:cxx:`hasThermodynamicForceProvider` and
:cxx:`hasInternalStateVariableProvider` search, among the behaviour
integrators defined on a given location, the one providing a given
quantity. The ``ParaviewExportIntegrationPointResultsAtNodes``
post-processing relies on them, so that a result can be exported on a
material carrying several behaviour integrators.

MPI utilities
-------------

The :cxx:`isMPIMainProcess`, :cxx:`getMPIRank` and :cxx:`getMPISize`
functions return the information for a given communicator or for the
communicator of a mesh discretization. The
:cxx:`isMPIMainProcessForDefaultCommunicator`,
:cxx:`getMPIRankForDefaultCommunicator` and
:cxx:`getMPISizeForDefaultCommunicator` functions use the default
communicator.

The collective :cxx:`synchronize` function gathers the error messages of
an execution context on all the processes of a communicator or of a mesh
discretization. It does nothing for a sequential mesh discretization.
A message shared by all the processes having one is kept unchanged.
Otherwise, each message is preceded by the rank of its process:

.. code-block:: text

   - Error on rank 1: <error message>
   - Error on rank 3: <error message>

Messages with several lines are preceded by a ``# Error on rank r`` title.

Issues fixed
============

- Issue 482: MPI functions on a MeshDiscretization and
  updateGridFunction with no function
- Issue 477: synchronize `Context` to merge error messages on MPI
  processes
- Issue 476: Fix :cxx:`assign_values`
- Issue 475: Error handling in :cxx:`updateGridFunction`
- Issue 474: Add function :cxx:`isMPIRootProcess` and remove tests
  :cxx:`getMPIrank() == 0`
- Issue 466: Refactor partial quadrature space to better support spaces
  on boundaries
- Issue 429: ex1 and ex3 sources shall only live in mm
- Issue 427: satoh, rve and mox examples are duplicates from mm-examples
- Issue 404: Remove deprecated methods of `AbstractNonLinearEvolutionProblem`
- Issue 292: Remove deprecated usage of `getMaterial` in
  `ParaviewExportIntegrationPointResultsAtNodesBase::getPartialQuadratureFunctionViews`
  and `ParaviewExportIntegrationPointResultsAtNodesBase::getResultDescription`
- Issue 290: Add support for partial quadratures spaces on boundaries
- Issue 288: master doesnt build without MPI
- Issue 286: updateGridFunction dilutes interface values when the
  functions do not cover all the materials of the mesh
- Issue 284: NonLinearModel does not export the initial state
- Issue 281: Nodal export of integration point results uses the node
  index as an integration point index
- Issue 279: ParaviewExportResults silently ignores unknown parameters
- Issue 277: NonLinearModel calls the deprecated solve and loses the
  caller's context
- Issue 275: UniformHeatSourceBoundaryCondition: only the last
  integration point contributes to the residual
- Issue 273: Fix forgotten change while parallel execution bug in
  behaviour integration
- Issue 270: FirstIterationConvergenceCriterion never triggers coupling
  iterations + wrong post-processing time
- Issue 268: Parallel execution bug in behaviour integration
- Issue 264: Installation from the doc fails on a fresh setup
- Issue 262: Refactor creation of submeshes and associated finite
  element spaces to handle it through MeshDiscretization and
  FiniteElementSpacesManager enhancement
- Issue 260: Refactor commented examples documentation
- Issue 257: Allow NewtonSolver to discard linear solver failures
- Issue 254: Missing parameter option for GMRESSolver
- Issue 253: Incomplete linear solver convergence checks in
  NonLinearEvolutionProblemImplementation.cxx and NewtonSolver.cxx
- Issue 248: Improve `PartialQuadratureFunction` interface
- Issue 245: Check the size of the unknown when defining bricks
- Issue 240: Small bug in `LinearSolverFactory.cxx`
- Issue 237: [cmake] Add a build-tests target
- Issue 218: [performance] synchronize success of the setup methods at a
  higher level to minimize collective communications enhancement
- Issue 213: Add a simple way to resolve dependencies (material
  properties, external state variables) of a :cxx:`NonLinearEvolutionProblem`
  using the gradients, thermodynamic forces and internal state variables
  of another one.
- Issue 211: Add coupling schemes and models
- Issue 209: Introduce the :cxx:`Simulation` class
- Issue 206: Line-search-like handling of behaviour integration failures
- Issue 200: automatically assign materials and boundaries's names from
  |mfem|'s attributes
- Issue 198: Add the ability to define multiple behaviour integrators on
  the same material
- Issue 193: Add support for other types of search operators for
  computing a prediction of the solution at the end of the time step
  enhancement
- Issue 192: Take external forces into account when computing the
  prediction of the solution
- Issue 188: retrieve information about a quadrature space
- Issue 149: work on the prediction of the solution
