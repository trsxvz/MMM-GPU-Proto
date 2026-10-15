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
   * - :cxx:`HyprePCG`, :cxx:`HypreGMRES`, :cxx:`HypreFGMRES`,
       :cxx:`HypreBiCGSTAB`, :cxx:`HypreCOGMRES`, :cxx:`HypreLGMRES`
       (parallel)
     - hypre
     - iterative
     - as hypre (host for :cxx:`HypreCOGMRES`)
     - hypre (applied by hypre, i.e. not :cxx:`HypreSmoother`)
     - :cxx:`VerbosityLevel`, :cxx:`Tolerance`, :cxx:`AbsoluteTolerance`,
       :cxx:`MaximumNumberOfIterations`, :cxx:`Preconditioner`,
       :cxx:`KDim` (GMRES solvers, 50 by default)
   * - :cxx:`HypreHybrid` (parallel)
     - hypre
     - iterative
     - as hypre
     - none
     - :cxx:`VerbosityLevel`, :cxx:`Tolerance` (relative, 1e-6 by
       default), :cxx:`MaximumNumberOfIterations` (of each phase, 1000 by
       default), :cxx:`SolverType` (:cxx:`PCG`, by default, :cxx:`GMRES`,
       :cxx:`BiCGSTAB`), :cxx:`KDim` (50 by default),
       :cxx:`ConvergenceRate` (0.9 by default), options of BoomerAMG
       (:cxx:`StrengthThreshold`, :cxx:`CoarseningType`,
       :cxx:`InterpolationType`, :cxx:`RelaxationType`,
       :cxx:`AggressiveCoarseningLevels`, :cxx:`MaximumNumberOfLevels`); no
       absolute tolerance
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
   * - :cxx:`SuperLUSolver` (parallel)
     - SuperLU_DIST
     - direct
     - host
     - none
     - :cxx:`ReuseAnalysis`, :cxx:`ColumnPermutation` (:cxx:`Natural`,
       :cxx:`MMD_ATA`, :cxx:`MMD_AT_PLUS_A`, :cxx:`COLAMD`,
       :cxx:`METIS_AT_PLUS_A`, :cxx:`ParMETIS`), :cxx:`RowPermutation`
       (:cxx:`None`, :cxx:`LargeDiag_MC64`, :cxx:`LargeDiag_HWPM`),
       :cxx:`IterativeRefinement` (:cxx:`None`, :cxx:`Single`,
       :cxx:`Double`, :cxx:`Extra`), :cxx:`Equilibrate`,
       :cxx:`ReplaceTinyPivot`, :cxx:`SymmetricPattern`,
       :cxx:`ParallelSymbolicFactorization`, :cxx:`NumberOfLookAheads`,
       :cxx:`DeviceOffload`, :cxx:`VerbosityLevel`
   * - :cxx:`STRUMPACKSolver`, :cxx:`STRUMPACKMixedPrecisionSolver`
       (parallel)
     - STRUMPACK
     - direct
     - host
     - none
     - :cxx:`Options` (see :ref:`strumpack_linear_solvers`),
       :cxx:`ReuseAnalysis`, :cxx:`VerbosityLevel`,
       :cxx:`RelativeTolerance`, :cxx:`AbsoluteTolerance`,
       :cxx:`MaximumNumberOfIterations`
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
   * - :cxx:`GinkgoCG`, :cxx:`GinkgoGMRES`, ... (see
       :ref:`ginkgo_linear_solvers`)
     - Ginkgo
     - iterative, direct (:cxx:`GinkgoDirect`)
     - as MFEM
     - any (Ginkgo only for :cxx:`GinkgoBiCG` and :cxx:`GinkgoPipeCG`)
     - :cxx:`Configuration`, :cxx:`VerbosityLevel`,
       :cxx:`RelativeTolerance`, :cxx:`AbsoluteTolerance`,
       :cxx:`MaximumNumberOfIterations`, :cxx:`Preconditioner`, :cxx:`KDim`
       (GMRES, CB-GMRES and GCR solvers)
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
       :cxx:`TruncationFactor`, :cxx:`MaximumElementsPerRow`,
       :cxx:`AggressiveInterpolationType`, :cxx:`NumberOfPaths`,
       :cxx:`MaximumCoarseSize`, :cxx:`MinimumCoarseSize`,
       :cxx:`RelaxationWeight`, :cxx:`ChebyshevOrder`,
       :cxx:`KeepTranspose`, :cxx:`VerbosityLevel`
   * - :cxx:`HypreILU` (parallel)
     - hypre
     - as hypre
     - :cxx:`HypreILUType`, :cxx:`HypreILULevelOfFill`, :cxx:`VerbosityLevel`
   * - :cxx:`HypreDiagScale` (parallel)
     - hypre
     - as hypre
     - :cxx:`VerbosityLevel`
   * - :cxx:`HypreSmoother` (parallel)
     - hypre
     - as hypre
     - :cxx:`Type` (:cxx:`Jacobi`, :cxx:`l1Jacobi`, :cxx:`l1GS`,
       :cxx:`l1GStr`, :cxx:`lumpedJacobi`, :cxx:`GS`, :cxx:`OPFS`,
       :cxx:`Chebyshev`, :cxx:`Taubin`, :cxx:`FIR`),
       :cxx:`RelaxationSweeps`, :cxx:`RelaxationWeight`, :cxx:`Omega`,
       :cxx:`PolynomialOrder`, :cxx:`PolynomialFraction`,
       :cxx:`EigenvalueEstimationIterations` (10 by default, 0 selecting
       the estimate of hypre), :cxx:`TaubinLambda`, :cxx:`TaubinMu`,
       :cxx:`TaubinIterations`, :cxx:`Window`; applied by MFEM, it can't
       precondition the Krylov solvers of hypre; the lumped Jacobi smoother
       divides by the sums of the rows, which vanish for a stiffness
       matrix; the :cxx:`OPFS`, :cxx:`Taubin`, :cxx:`FIR` and
       :cxx:`lumpedJacobi` smoothers only compute on the host
   * - :cxx:`HypreFSAI` (parallel)
     - hypre
     - as hypre
     - :cxx:`AlgorithmType`, :cxx:`LocalSolveType`, :cxx:`MaximumSteps`,
       :cxx:`MaximumStepSize`, :cxx:`MaximumNonZerosPerRow`,
       :cxx:`NumberOfLevels`, :cxx:`Threshold`, :cxx:`KapTolerance`,
       :cxx:`EigenvalueEstimationIterations`, :cxx:`VerbosityLevel`
   * - :cxx:`HypreEuclid`, :cxx:`HypreParaSails` (parallel)
     - hypre
     - host
     - :cxx:`VerbosityLevel`
   * - :cxx:`AmgX`
     - AmgX
     - GPU
     - :cxx:`Configuration`, :cxx:`ConfigurationFile`, :cxx:`VerbosityLevel`,
       :cxx:`AmgXParameters`, :cxx:`Precision`
   * - :cxx:`GinkgoJacobi`, :cxx:`GinkgoILU`, :cxx:`GinkgoMultigrid`, ...
       (see :ref:`ginkgo_linear_solvers`)
     - Ginkgo
     - as MFEM
     - :cxx:`Configuration`, :cxx:`L1Smoother` (local preconditioners in
       parallel)
   * - :cxx:`PETScPreconditioner` (parallel)
     - PETSc
     - host
     - :cxx:`Options` (see :ref:`petsc_linear_solvers`)
   * - :cxx:`STRUMPACKPreconditioner` (parallel)
     - STRUMPACK
     - host
     - :cxx:`Options` (see :ref:`strumpack_linear_solvers`),
       :cxx:`ReuseAnalysis`, :cxx:`VerbosityLevel`

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

.. _strumpack_linear_solvers:

Direct solvers of SuperLU_DIST and STRUMPACK
============================================

The matrices of hypre are copied to the formats of SuperLU_DIST and
STRUMPACK, which factorize them on the host (SuperLU_DIST may offload
parts of the factorization to the GPU, see :cxx:`DeviceOffload`).

The :cxx:`ReuseAnalysis` option of :cxx:`SuperLUSolver` keeps the column
permutation of the previous factorization for the next operators, which
must have the same sparsity pattern (as the successive operators of a
Newton solver); the solver is built again when the size of the operator
changes. The parallel symbolic factorization of SuperLU_DIST, selected by
the :cxx:`ParallelSymbolicFactorization` option or by the :cxx:`ParMETIS`
column permutation (unless the :cxx:`ParallelSymbolicFactorization` option
is false), requires at least two processes and can't be combined with the
reuse of the analysis.

The options of STRUMPACK are given by the :cxx:`Options` parameter, by
their names on the command line of STRUMPACK without the leading dashes,
for instance:

.. code-block:: cpp

   problem.setLinearSolver(
       ctx, "STRUMPACKSolver",
       {{"RelativeTolerance", 1e-10},
        {"Options", Parameters{{"sp_compression", "blr"},
                               {"blr_rel_tol", 1e-2},
                               {"sp_Krylov_solver", "pgmres"}}}}) |
       or_die;

The :cxx:`sp_compression` option compresses the frontal matrices (block
low-rank, HSS, HODLR, lossy compression), the factorization being then
approximate and used as the preconditioner of a Krylov solver of STRUMPACK
(:cxx:`sp_Krylov_solver`, whose tolerances are given by the options of
mfem-mgis). :cxx:`STRUMPACKMixedPrecisionSolver` factorizes in single
precision and refines the solution in double precision (up to the relative
tolerance of STRUMPACK, 1e-6, unless :cxx:`RelativeTolerance` is given);
STRUMPACK refuses a Krylov solver without preconditioner as its outer
solver. The :cxx:`STRUMPACKPreconditioner` preconditioner applies an
approximate factorization once (block low-rank compression by default) and
may be used by the Krylov solvers of MFEM. The :cxx:`ReuseAnalysis` option
keeps the reordering of the matrix for the next operators of the same
sparsity pattern; the solver is built again when the size of the operator
changes.

An option without value (a flag) is given by an empty string or by the
boolean :cxx:`true`, the boolean :cxx:`false` giving nothing. The options
unknown to STRUMPACK, misspelled ones included, are ignored without
message. The :cxx:`VerbosityLevel`, :cxx:`RelativeTolerance`,
:cxx:`AbsoluteTolerance` and :cxx:`MaximumNumberOfIterations` options of
mfem-mgis are applied after the :cxx:`Options` parameter, and override the
:cxx:`sp_verbose`, :cxx:`sp_rel_tol`, :cxx:`sp_abs_tol` and
:cxx:`sp_maxit` options of STRUMPACK. STRUMPACK reports a success when its
Krylov solver or its iterative refinement stops at the maximum number of
iterations: mfem-mgis reports such a resolution as not converged.

.. _ginkgo_linear_solvers:

Linear solvers and preconditioners of Ginkgo
============================================

The linear solvers and the preconditioners of Ginkgo are described by
configurations of Ginkgo (see the documentation of
:cxx:`gko::config::parse`), built on the classes of MFEM wrapping Ginkgo,
which convert the matrices and the vectors of MFEM. Ginkgo computes where
the device of MFEM computes: on the GPU (CUDA or HIP executor of Ginkgo)
when MFEM computes on a GPU, with OpenMP when the device of MFEM uses
OpenMP and Ginkgo is built with OpenMP, and with the reference executor
of Ginkgo otherwise. The solvers and the preconditioners are available in
parallel when Ginkgo is built with MPI.

Each entry is a configuration of a given type, completed by the keys of
its :cxx:`Configuration` option, which gives access to all the parameters
of Ginkgo:

.. code-block:: cpp

   problem.setLinearSolver(
       ctx, "GinkgoGMRES",
       {{"RelativeTolerance", 1e-10},
        {"MaximumNumberOfIterations", 500},
        {"KDim", 50},
        {"Configuration", Parameters{{"flexible", true}}},
        {"Preconditioner",
         Parameters{{"Name", "GinkgoILU"},
                    {"Options",
                     Parameters{{"Configuration",
                                 Parameters{{"factorization",
                                             Parameters{{"type",
                                                         "factorization::ParIlu"}}}}}}}}}}) |
       or_die;

The keys and the values of a configuration are the ones of Ginkgo
(:cxx:`snake_case` keys, a nested configuration being given by a map of
parameters, and a list by a vector of parameters). A floating-point value
must be written as such (:cxx:`1.` and not :cxx:`1`). Ginkgo refuses the
keys it does not know. The keys defined by an entry can't be redefined:
its type, the Krylov dimension given by the :cxx:`KDim` option, the
triangular solvers of :cxx:`GinkgoILUISAI` and :cxx:`GinkgoICISAI`. The
stopping criteria (:cxx:`criteria` key) and the preconditioner
(:cxx:`preconditioner` key, :cxx:`solver` for the iterative refinement)
can't be given both by the configuration and by the options of mfem-mgis.
The factorization of :cxx:`GinkgoDirect` (LU by default, for instance
:cxx:`factorization::Cholesky`) and the levels (:cxx:`mg_level`) and the
number of cycles (:cxx:`criteria`) of :cxx:`GinkgoMultigrid` can be
redefined. Each key of a configuration is checked alone before the
configuration is used, since Ginkgo stops the program on some errors in a
configuration whose other keys are not read yet.

Linear solvers
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Name
     - Solver of Ginkgo
     - Availability
   * - :cxx:`GinkgoCG`, :cxx:`GinkgoFCG`
     - conjugate gradient, flexible variant
     - sequential, parallel
   * - :cxx:`GinkgoBiCGSTAB`, :cxx:`GinkgoCGS`
     - BiCGSTAB, conjugate gradient squared
     - sequential, parallel
   * - :cxx:`GinkgoGMRES`, :cxx:`GinkgoGCR`
     - GMRES, generalized conjugate residual (:cxx:`KDim` option)
     - sequential, parallel
   * - :cxx:`GinkgoMINRES`
     - MINRES
     - sequential, parallel
   * - :cxx:`GinkgoIR`
     - iterative refinement: stationary iterations with the
       preconditioner as inner solver
     - sequential, parallel
   * - :cxx:`GinkgoIDR`
     - induced dimension reduction
     - sequential
   * - :cxx:`GinkgoPipeCG`
     - pipelined conjugate gradient (preconditioners of Ginkgo and
       assembled operators only)
     - sequential
   * - :cxx:`GinkgoBiCG`
     - BiCG (preconditioners of Ginkgo which can be transposed:
       :cxx:`GinkgoJacobi`, :cxx:`GinkgoILU`, :cxx:`GinkgoILUISAI`,
       :cxx:`GinkgoIC`, :cxx:`GinkgoICISAI`; assembled operators only)
     - sequential
   * - :cxx:`GinkgoCBGMRES`
     - GMRES storing its Krylov basis in a reduced precision
       (:cxx:`KDim` option, :cxx:`storage_precision` key)
     - sequential
   * - :cxx:`GinkgoDirect`
     - sparse LU factorization (direct solver)
     - sequential
   * - :cxx:`GinkgoSolver`
     - any iterative solver of Ginkgo (:cxx:`solver::Cg`,
       :cxx:`solver::Chebyshev`, :cxx:`solver::Multigrid`, ...), described
       by the whole configuration given by the :cxx:`Configuration` option
       (with its :cxx:`type` key)
     - as the solver

The iterative solvers accept the :cxx:`RelativeTolerance` (relative to
the norm of the initial residual), :cxx:`AbsoluteTolerance`,
:cxx:`MaximumNumberOfIterations` (ten by default, as for the solvers of
MFEM) and :cxx:`VerbosityLevel` options, unless the configuration defines
its own stopping criteria (:cxx:`criteria` key). The norm of the residual
computed by the solver is used, except for :cxx:`GinkgoMINRES`, which uses
the norm of the preconditioned residual (MINRES does not compute the
residual). The relative tolerance of a solver whose stopping criteria are
given by the options of mfem-mgis can be changed (see
:cxx:`setRelativeTolerance`). The iterative solvers accept any
preconditioner, except the restrictions above:

- a preconditioner of Ginkgo is given to the solver by its factory, Ginkgo
  generating it for each operator;
- a preconditioner of another library is wrapped by MFEM, and receives the
  operators given to the solver.

Preconditioners
---------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Name
     - Preconditioner of Ginkgo
   * - :cxx:`GinkgoJacobi`
     - block Jacobi (blocks of at most 32 unknowns by default, the
       :cxx:`max_block_size` key being one for the scalar Jacobi method)
   * - :cxx:`GinkgoILU`, :cxx:`GinkgoIC`
     - incomplete LU and Cholesky factorizations (:cxx:`factorization`
       key for the parallel variants of Ginkgo)
   * - :cxx:`GinkgoILUISAI`, :cxx:`GinkgoICISAI`
     - the same, with triangular factors applied by incomplete sparse
       approximate inverses
   * - :cxx:`GinkgoSOR`, :cxx:`GinkgoGaussSeidel`
     - successive over-relaxation, Gauss-Seidel
   * - :cxx:`GinkgoDirect`
     - sparse LU factorization
   * - :cxx:`GinkgoMultigrid`
     - one cycle of the algebraic multigrid of Ginkgo (parallel graph
       match aggregation)
   * - :cxx:`GinkgoPreconditioner`
     - any preconditioner of Ginkgo, described by the whole configuration
       given by the :cxx:`Configuration` option

In parallel, all the preconditioners except :cxx:`GinkgoMultigrid` and
:cxx:`GinkgoPreconditioner` are local: they are the local solvers of the
Schwarz preconditioner of Ginkgo (one block per process). Its L1 smoother
is selected by the :cxx:`L1Smoother` option, which sorts the columns of
the diagonal block of the matrix of hypre in place (the matrix is then
reassembled before being used by hypre). The configuration of
:cxx:`GinkgoPreconditioner` must then describe a distributed
preconditioner (for instance :cxx:`preconditioner::Schwarz`), without the
L1 smoother. The preconditioners of Ginkgo require assembled operators.
