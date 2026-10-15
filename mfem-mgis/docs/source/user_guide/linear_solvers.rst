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
   * - :cxx:`PETScLinearSolver` (parallel)
     - PETSc
     - iterative
     - host
     - any
     - :cxx:`Options` (see :ref:`petsc_linear_solvers`),
       :cxx:`VerbosityLevel`, :cxx:`RelativeTolerance`,
       :cxx:`AbsoluteTolerance`, :cxx:`MaximumNumberOfIterations`,
       :cxx:`Preconditioner`

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
   * - :cxx:`PETScPreconditioner` (parallel)
     - PETSc
     - host
     - :cxx:`Options` (see :ref:`petsc_linear_solvers`)

The libraries which are not available (MFEM built without them, or
mfem-mgis built without AmgX) register nothing.

.. _petsc_linear_solvers:

Linear solvers and preconditioners of PETSc
===========================================

The :cxx:`PETScLinearSolver` linear solver is a Krylov solver of PETSc
(:cxx:`KSP`), and the :cxx:`PETScPreconditioner` preconditioner a
preconditioner of PETSc (:cxx:`PC`). Both are configured by the options
of PETSc, given by the :cxx:`Options` parameter without the leading dash:

.. code-block:: cpp

   problem.setLinearSolver(
       ctx, "PETScLinearSolver",
       {{"RelativeTolerance", 1e-10},
        {"Options", Parameters{{"ksp_type", "cg"},
                               {"pc_type", "gamg"},
                               {"pc_gamg_threshold", 0.02},
                               {"ksp_monitor", ""}}}}) |
       or_die;

The value of an option is a boolean, an integer, a floating-point value
or a string, an empty string giving an option without value (such as
:cxx:`ksp_monitor`). The options are written in the options database of
PETSc with a prefix which is unique to each solver, so that the options of
two solvers do not interfere, and removed from it when the solver is
destroyed. They are read by PETSc when the solver is first used: they
thus override the tolerances and the maximum number of iterations given by
the parameters of mfem-mgis, and an invalid value is only detected then.
PETSc ignores the options it does not know: the options which PETSc never
read (for instance a misspelt option) are reported when a solver which
has been used is destroyed.

Without options, the defaults of PETSc apply: the GMRES solver restarted
every 30 iterations, preconditioned by an incomplete LU factorization in
sequential, and by the block Jacobi method with incomplete LU
factorizations of the blocks in parallel. The convergence is checked on
the norm of the preconditioned residual (:cxx:`ksp_norm_type`).

The operators are converted to the AIJ format of PETSc each time they are
set, so that most preconditioners of PETSc can be built, and the
preconditioner of PETSc is then built again: the operators shall be
assembled (an operator applied without assembly is converted column by
column).

A preconditioner of mfem-mgis may be given by the :cxx:`Preconditioner`
parameter of :cxx:`PETScLinearSolver`:

- a preconditioner of another library is seen by PETSc as a shell
  preconditioner (:cxx:`pc_type shell`);
- the :cxx:`PETScPreconditioner` preconditioner becomes the preconditioner
  of the Krylov solver, its options being given by its own :cxx:`Options`
  option.

In both cases, the options of the preconditioner of PETSc (beginning with
:cxx:`pc_`) are refused in the options of the solver, since they would
replace the given preconditioner. The :cxx:`PETScPreconditioner`
preconditioner may also be used by the Krylov solvers of MFEM:

.. code-block:: cpp

   problem.setLinearSolver(
       ctx, "CGSolver",
       {{"RelativeTolerance", 1e-10},
        {"Preconditioner",
         Parameters{{"Name", "PETScPreconditioner"},
                    {"Options",
                     Parameters{{"Options",
                                 Parameters{{"pc_type", "gamg"}}}}}}}}) |
       or_die;

PETSc is initialized when the first linear solver or preconditioner of
PETSc is built, which shall thus be built by all the processes. It is
then finalized by the :cxx:`finalize` function, unless it was initialized
by the program itself. The command line is not given to PETSc: the
options of the solvers are only given by their :cxx:`Options` parameter.
The non linear solver of PETSc, selected by the :cxx:`--use-petsc`
command line option, is not related to these solvers.
