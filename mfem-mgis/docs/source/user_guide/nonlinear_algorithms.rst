.. _mfem_mgis_nonlinear_algorithms:

==============================
Available nonlinear algorithms
==============================

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Prediction of the solution
==========================

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

The members of the :cxx:`PredictionPolicy` structure are:

- :cxx:`strategy`: the prediction strategy
  (:cxx:`PredictionStrategy::DEFAULT_PREDICTION` by default);
- :cxx:`prediction_operator`: the operator used by the
  :cxx:`BEGINNING_OF_TIME_STEP_PREDICTION` strategy
  (:cxx:`PredictionOperator::ELASTIC` by default);
- :cxx:`integration_operator`: the operator used by the
  :cxx:`CONSTANT_GRADIENTS_INTEGRATION_PREDICTION` strategy
  (:cxx:`IntegrationOperator::CONSISTENT_TANGENT` by default);
- :cxx:`null_time_increment`: if true, the behaviours are integrated with a
  null time increment by the :cxx:`CONSTANT_GRADIENTS_INTEGRATION_PREDICTION`
  strategy, which neglects the viscoplastic relaxation (:cxx:`false` by
  default).

Available strategies
--------------------

Default prediction (:cxx:`PredictionStrategy::DEFAULT_PREDICTION`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, a nonlinear evolution problem uses the solution at the
beginning of the time step, modified by applying Dirichlet boundary
conditions, as the initial guess of the solution at the end of the time
step.

.. warning::

   In mechanics, this may lead to very high increments of the deformation
   gradients or the strain in the neighboring elements of boundaries where
   evolving displacements are imposed.
   
Prediction for the state at the beginning of the time step (:cxx:`PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :cxx:`CONSTANT_GRADIENTS_INTEGRATION` strategy determines the
increment of the unknown :math:`\Delta\,{u}` by solving the following
linear system:
     
.. math::

   \tilde{\mathbb{K}}\,\cdot\,\Delta\,\mathbb{u} =
   \ets{\mathbb{F}_{e}}-\ets{\tilde{\mathbb{F}}_{i}}
     
where:
     
- :math:`\tilde{\mathbb{K}}` denotes the operator computed at the end of
  the behaviour integration, selected by the :cxx:`integration_operator`
  member of the prediction policy (see below).
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
  discarded by integrating the behaviour with a null time step
  (:cxx:`null_time_increment` member of the prediction policy set to
  :cxx:`true`).

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

.. note::

   Internally, the :cxx:`computePrediction` methods of the
   :cxx:`NonLinearEvolutionProblemImplementation` classes compute the
   opposite of the prediction.

Newton algorithm
================

The default nonlinear solver is the Newton's algorithm.

Let :math:`\mathbb{u}^{(n)}` be the estimate of the solution at the
:math:`n^{\text{th}}` iteration. The Newton's algorithm determines a new
estimate :math:`\mathbb{u}^{(n+1)}` of the solution by computing a
correction :math:`\delta\,\mathbb{u}^{(n)}` as follows:

.. math::

   \mathbb{K}\,\cdot\,\delta\,\mathbb{u}^{(n)} = \ets{\mathbb{F}_{e}}^{(n)}-\ets{\mathbb{F}_{i}}^{(n)}

where:

- :math:`\mathbb{K}` denotes one of the search operator.
- :math:`\ets{\mathbb{F}_{e}}^{(n)}` denotes current estimate of the
  external forces at the end of the time step.
- :math:`\ets{\mathbb{F}_{i}}^{(n)}` denotes current estimate of the inner
  forces at the end of the time step.

Currently, the search operator is given by the consistent tangent operator.

The new estimate :math:`\mathbb{u}^{(n+1)}` is calculated as follows:

.. math::

   \mathbb{u}^{(n+1)} = \mathbb{u}^{(n)}+\delta\,\mathbb{u}^{(n)}

Handling of integration failure
-------------------------------

When the new estimate leads to an integration failure, a new estimate of
the solution is calculated by keeping the same search direction but by
reducing the amplitude of the correction by a factor :math:`2` as
follows:

.. math::

   \mathbb{u}^{(n+1)} = \mathbb{u}^{(n)}+\frac{1}{2}\,\delta\,\mathbb{u}^{(n)}

This reduction of the amplitude is repeated until an estimation of the
solution that does not lead to an integration failure is found, or until
the maximum number of iterations is reached (each reduction counts as an
iteration).

Stagnation of the residual
--------------------------

The residual is computed with rounding errors: its norm cannot go below a
floor, which grows with the size of the problem and with the level of the
stresses, and which is much higher with ``float`` than with ``double``
(about :math:`1.3\,10^{-5}` to :math:`2\,10^{-5}` of the reference
residual for a beam of :math:`10^5` unknowns in small strain, :math:`3\,10^{-4}`
to :math:`8\,10^{-4}` for a cube of one element in finite strain). A
relative tolerance below this floor can't be reached.

The Newton solver therefore accepts as converged a residual which has not
decreased by more than 10 % during the last two iterations, when its norm,
relative to the reference residual (the residual of the prediction, if
any, or else the initial residual of the time step), is below the option
:cxx:`StagnationRelativeTolerance`. The best iterate of the resolution is
then kept: the residual at the floor is noisy, and the last iterate may be
worse than a previous one. A slow convergence, even with a factor 0.9 per
iteration, is not taken for a stagnation. The default value of
the option is :math:`10^4` times the machine epsilon of the numerical type:
about :math:`1.2\,10^{-3}` with ``float``, :math:`2.2\,10^{-12}` with
``double``. A null value disables this criterion. A message is printed
when the criterion stops the iterations. Finite strain computations with
``float`` on large meshes may have a floor close to this default value:
the H = F - I formulation, not available yet, would lower it. The
:cxx:`FixedPoint` solver does not have this criterion.

.. code-block:: c++

   problem.setSolverParameters(ctx, {{"RelativeTolerance", 1e-8},
                                     {"StagnationRelativeTolerance", 1e-4}}) |
       or_die;

.. _mfem_mgis_fixed_point_algorithm:

Fixed point algorithm
=====================

The nonlinear solver is chosen by the :cxx:`NonLinearSolver` parameter of
the nonlinear evolution problem: :cxx:`Newton` (default) or
:cxx:`FixedPoint`, given with its options, which can be changed later by
the :cxx:`setSolverParameters` method:

.. code-block:: c++

   {"NonLinearSolver",
    mfem_mgis::Parameters{{"FixedPoint", mfem_mgis::Parameters{...}}}}

The :cxx:`FixedPoint` solver performs the following iterations:

.. math::

   \mathbb{u}^{(n+1)} = \mathbb{u}^{(n)}-\mathbb{S}\,\cdot\,\mathbb{r}\left(\mathbb{u}^{(n)}\right)

where :math:`\mathbb{r}` denotes the residual. The iterations may be
accelerated (option :cxx:`Acceleration`, for instance :cxx:`Anderson` or
:cxx:`Cast3M`, none by default).

Methods
-------

The option :cxx:`Method` selects :math:`\mathbb{S}`:

- :cxx:`ElasticOperator` (default): an exact or approximate inverse of the
  elastic operator :math:`\mathbb{K}_{0}`, assembled from the elastic
  operators of the behaviours (option :cxx:`ElasticOperatorComputation`:
  :cxx:`Prediction`, by default, uses the prediction operator of the
  behaviours; :cxx:`Integration` integrates the behaviours and uses their
  elastic operator, for the behaviours which do not provide the prediction
  operator) at each resolution
  or once (option :cxx:`ElasticOperatorUpdate`: :cxx:`EachResolution`, by
  default, or :cxx:`Once`). The behaviours are integrated without tangent
  operator at each iteration.
- :cxx:`ModifiedNewton`: the inverse of the consistent tangent operator at
  the first iteration of each resolution, given to the linear solver of the
  problem.
- :cxx:`NewtonK0`: Newton's method, the consistent tangent operator being
  inverted by an iterative solver of :code:`MFEM` (option
  :cxx:`TangentLinearSolver`: :cxx:`CGSolver`, by default,
  :cxx:`GMRESSolver`, :cxx:`FGMRESSolver`, :cxx:`BiCGSTABSolver`,
  :cxx:`MINRESSolver` or :cxx:`SLISolver`, given with its options except
  :cxx:`Preconditioner`, with a relative tolerance of :math:`10^{-4}` and
  at most 1000 iterations if not given),
  preconditioned by the linear solver of the elastic operator (option
  :cxx:`ElasticOperatorLinearSolver`, required), whose setup is kept: for
  instance :cxx:`HypreBoomerAMG` set up on :math:`\mathbb{K}_{0}` instead
  of a new setup at each iteration. The tangent operator can be applied
  without assembly (see :ref:`mfem_mgis_batch_computations_matrix_free`).

The linear solver of the elastic operator (option
:cxx:`ElasticOperatorLinearSolver`) is one of the linear solvers of the
problems (see :doc:`glossary`), for instance :cxx:`HypreBoomerAMG` for one
cycle of BoomerAMG per iteration in a parallel computation. With the
:cxx:`ElasticOperator` method, if it is not given, the linear solver of the
problem is used, and set up with :math:`\mathbb{K}_{0}` at each
resolution. The value :cxx:`None` removes the linear solver of the elastic
operator (option :cxx:`ElasticOperatorLinearSolver`) or the acceleration
(option :cxx:`Acceleration`). The :cxx:`NewtonK0` method can't be
accelerated, and the :cxx:`ModifiedNewton` method refuses the option
:cxx:`ElasticOperatorLinearSolver`.

The option :cxx:`PredictionLinearSolver` gives the linear solver of the
elastic prediction of the problem: :cxx:`Problem` (default), the linear
solver of the problem, or an iterative solver of :code:`MFEM`
preconditioned by the linear solver of the elastic operator, given as for
:cxx:`TangentLinearSolver` (:cxx:`ElasticOperator` being a synonym of
:cxx:`CGSolver`), with a relative tolerance of :math:`10^{-8}` and at most
1000 iterations if not given. The latter is used by the elastic
predictions made at the beginning of the time steps which follow the first
resolution, with the :cxx:`ElasticOperator` and :cxx:`NewtonK0` methods and
a linear solver of the elastic operator.

The other options are the ones of the iterative solvers
(:cxx:`RelativeTolerance`, :cxx:`AbsoluteTolerance`,
:cxx:`MaximumNumberOfIterations`, :cxx:`VerbosityLevel`) and
:cxx:`DiscardLinearSolverFailure`.

The relative tolerance of the iterative linear solver of the Newton solver
and of the :cxx:`NewtonK0` method can be adapted to the decrease of the
residual (option :cxx:`AdaptiveLinearSolverTolerance`, see
:doc:`glossary`). This option is refused by the other methods of the
:cxx:`FixedPoint` solver.

When the integration of the behaviours fails at a new estimate, this
estimate is replaced by the middle of this estimate and of the previous
one, ten times at most. The iterations stop if the norm of the residual is
not finite.

The following options select the :cxx:`NewtonK0` method in a parallel
computation:

.. code-block:: c++

   const auto newton_k0 = mfem_mgis::Parameters{
       {"Method", "NewtonK0"},
       {"ElasticOperatorUpdate", "Once"},
       {"ElasticOperatorLinearSolver", "HypreBoomerAMG"},
       {"TangentLinearSolver",
        mfem_mgis::Parameters{
            {"CGSolver",
             mfem_mgis::Parameters{{"RelativeTolerance", 1e-4},
                                   {"MaximumNumberOfIterations", 1000}}}}}};
   // given to the problem as
   // {"NonLinearSolver", mfem_mgis::Parameters{{"FixedPoint", newton_k0}}}
