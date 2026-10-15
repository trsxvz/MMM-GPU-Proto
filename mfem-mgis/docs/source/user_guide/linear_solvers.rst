.. _linear_solvers:

=======================================
Linear solvers and preconditioners
=======================================

The linear solver of a problem is chosen by its name and its options,
given to the :cxx:`setLinearSolver` method of the problem:

.. code-block:: cpp

   problem.setLinearSolver(
       ctx, "CGSolver",
       {{"RelativeTolerance", 1e-10},
        {"MaximumNumberOfIterations", 500},
        {"Preconditioner", Parameters{{"Name", "HypreBoomerAMG"},
                                      {"Options", Parameters{}}}}}) |
       or_die;

The linear solvers are registered in the factory of linear solvers
(:cxx:`LinearSolverFactory`), the preconditioners in the factory of
preconditioners (:cxx:`PreconditionerFactory`). A preconditioner is
selected by the :cxx:`Preconditioner` option of a linear solver: its
name (:cxx:`Name`, :cxx:`None` for no preconditioner) and its options
(:cxx:`Options`).

Each linear solver and each preconditioner is registered with a
description (:cxx:`LinearSolverDescription`), returned by the
:cxx:`getDescription` method of the factories (the :cxx:`getNames`
method lists the registered names):

- the library which provides it;
- its kind: iterative solver (converging up to its tolerances), direct
  solver, approximate solver (a fixed number of cycles of a multigrid,
  used for instance by the fixed point solvers) or preconditioner;
- where it computes: on the host only, on the GPU only, where the device
  of MFEM computes, or where the memory of hypre is. When the execution
  settings place the linear solvers on the CPU or on the GPU (see
  :ref:`mfem_mgis_batch_computations`), a linear solver or a preconditioner which
  computes elsewhere is refused with a message;
- the preconditioners accepted by a linear solver: none (direct solvers,
  self-contained solvers), any preconditioner, or only the ones of its
  library (the Krylov solvers of hypre only accept the preconditioners of
  hypre). The other combinations are refused with a message;
- its options.

Linear solvers
==============

The solvers marked as parallel are only available in parallel.

.. list-table::
   :header-rows: 1
   :widths: 18 10 12 12 14 34

   * - Name
     - Library
     - Kind
     - Computes
     - Preconditioners
     - Options
   * - :cxx:`CGSolver`, :cxx:`GMRESSolver`, :cxx:`FGMRESSolver`,
       :cxx:`BiCGSTABSolver`, :cxx:`MINRESSolver`, :cxx:`SLISolver`
     - MFEM
     - iterative
     - as MFEM
     - any
     - :cxx:`VerbosityLevel`, :cxx:`RelativeTolerance`,
       :cxx:`AbsoluteTolerance`, :cxx:`MaximumNumberOfIterations`,
       :cxx:`Preconditioner`, :cxx:`KDim` (GMRES solvers)
   * - :cxx:`HyprePCG`, :cxx:`HypreGMRES`, :cxx:`HypreFGMRES` (parallel)
     - hypre
     - iterative
     - as hypre
     - hypre
     - :cxx:`VerbosityLevel`, :cxx:`Tolerance`, :cxx:`AbsoluteTolerance`,
       :cxx:`MaximumNumberOfIterations`, :cxx:`Preconditioner`,
       :cxx:`KDim` (GMRES solvers)
   * - :cxx:`HypreBoomerAMG` (parallel)
     - hypre
     - approximate
     - as hypre
     - none
     - the options of the :cxx:`HypreBoomerAMG` preconditioner
   * - :cxx:`MUMPSSolver` (parallel)
     - MUMPS
     - direct
     - host
     - none
     - :cxx:`Symmetric`, :cxx:`PositiveDefinite`, :cxx:`ReuseAnalysis`
   * - :cxx:`CuDSSSolver` (parallel)
     - cuDSS
     - direct
     - GPU
     - none
     - :cxx:`Symmetric`, :cxx:`PositiveDefinite`, :cxx:`ReuseAnalysis`
   * - :cxx:`UMFPackSolver` (sequential)
     - SuiteSparse
     - direct
     - host
     - none
     -
   * - :cxx:`AmgX`
     - AmgX
     - approximate
     - GPU
     - none
     - :cxx:`Configuration`, :cxx:`ConfigurationFile`, :cxx:`VerbosityLevel`,
       :cxx:`AmgXParameters`, :cxx:`Precision`
   * - :cxx:`AmgXPCG`
     - AmgX
     - iterative
     - GPU
     - none
     - :cxx:`VerbosityLevel`, :cxx:`RelativeTolerance`,
       :cxx:`AbsoluteTolerance`, :cxx:`MaximumNumberOfIterations`,
       :cxx:`AmgXParameters`, :cxx:`Precision`

Preconditioners
===============

.. list-table::
   :header-rows: 1
   :widths: 22 12 12 54

   * - Name
     - Library
     - Computes
     - Options
   * - :cxx:`OperatorJacobiSmoother`
     - MFEM
     - as MFEM
     - :cxx:`Damping`
   * - :cxx:`OperatorChebyshevSmoother`
     - MFEM
     - as MFEM
     - :cxx:`Order`, :cxx:`NumberOfPowerIterations`
   * - :cxx:`AssembledOperator`
     - mfem-mgis
     - as the wrapped preconditioner
     - :cxx:`Preconditioner` (the wrapped preconditioner, set up on the
       assembled version of an operator applied without assembly),
       :cxx:`UpdatePeriod`
   * - :cxx:`HypreBoomerAMG` (parallel)
     - hypre
     - as hypre
     - :cxx:`Strategy`, :cxx:`StrengthThreshold`, :cxx:`CoarseningType`,
       :cxx:`AggressiveCoarseningLevels`, :cxx:`InterpolationType`,
       :cxx:`RelaxationType`, :cxx:`RelaxationSweeps`, :cxx:`CycleType`,
       :cxx:`NumberOfCycles`, :cxx:`MaximumNumberOfLevels`,
       :cxx:`VerbosityLevel`
   * - :cxx:`HypreILU` (parallel)
     - hypre
     - as hypre
     - :cxx:`HypreILUType`, :cxx:`HypreILULevelOfFill`, :cxx:`VerbosityLevel`
   * - :cxx:`HypreDiagScale` (parallel)
     - hypre
     - as hypre
     - :cxx:`VerbosityLevel`
   * - :cxx:`HypreEuclid`, :cxx:`HypreParaSails` (parallel)
     - hypre
     - host
     - :cxx:`VerbosityLevel`
   * - :cxx:`AmgX`
     - AmgX
     - GPU
     - :cxx:`Configuration`, :cxx:`ConfigurationFile`, :cxx:`VerbosityLevel`,
       :cxx:`AmgXParameters`, :cxx:`Precision`

The libraries which are not available (MFEM built without them, or
mfem-mgis built without AmgX) register nothing.
