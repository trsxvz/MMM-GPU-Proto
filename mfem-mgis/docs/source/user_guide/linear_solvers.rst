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
     - :cxx:`Symmetric`, :cxx:`PositiveDefinite`, :cxx:`ReuseAnalysis`,
       :cxx:`ReorderingStrategy` (:cxx:`Automatic`, :cxx:`AMD`,
       :cxx:`AMF`, :cxx:`PORD`, :cxx:`METIS`, :cxx:`ParMETIS`,
       :cxx:`Scotch`, :cxx:`PTScotch`; MUMPS replaces a sequential
       strategy it was built without, while a missing parallel one,
       :cxx:`ParMETIS` or :cxx:`PTScotch`, may make MFEM stop the program),
       :cxx:`BLRTolerance` (block low-rank factorization of the fronts
       large enough, approximate: its result is not checked, small values
       such as 1e-10 keep a direct solver), :cxx:`VerbosityLevel` (print
       level of MUMPS, 0: nothing, 1: errors, by default, 2: warnings and
       statistics, up to 4; unlike the iterative solvers, 0 hides the
       errors)
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
   * - :cxx:`MKLPardisoSolver` (sequential),
       :cxx:`MKLClusterSparseSolver` (parallel)
     - oneMKL
     - direct
     - host
     - none
     - :cxx:`Symmetric`, :cxx:`PositiveDefinite`, :cxx:`ReuseAnalysis`,
       :cxx:`VerbosityLevel` (see :ref:`mkl_linear_solvers`)
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
   * - :cxx:`KLUSolver` (sequential)
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
       :cxx:`AmgXParameters`, :cxx:`Precision`, :cxx:`Multigrid`,
       :cxx:`NodalBlocks`
   * - :cxx:`AmgXPCG`
     - AmgX
     - iterative
     - GPU
     - none
     - :cxx:`VerbosityLevel`, :cxx:`RelativeTolerance`,
       :cxx:`AbsoluteTolerance`, :cxx:`MaximumNumberOfIterations`,
       :cxx:`AmgXParameters`, :cxx:`Precision`, :cxx:`Multigrid`,
       :cxx:`NodalBlocks`
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
   * - :cxx:`MAGMACG`, :cxx:`MAGMABiCGSTAB`, :cxx:`MAGMAGMRES`,
       :cxx:`MAGMAIDR`, :cxx:`MAGMACGS`, :cxx:`MAGMAQMR`,
       :cxx:`MAGMATFQMR` (sequential)
     - MAGMA
     - iterative
     - GPU
     - MAGMA
     - :cxx:`VerbosityLevel`, :cxx:`RelativeTolerance`,
       :cxx:`AbsoluteTolerance`, :cxx:`MaximumNumberOfIterations`,
       :cxx:`Preconditioner`, :cxx:`KDim` (GMRES),
       :cxx:`ShadowSpaceDimension` (IDR) (see :ref:`magma_linear_solvers`)
   * - :cxx:`PETScLinearSolver` (parallel)
     - PETSc
     - iterative
     - host
     - any
     - :cxx:`Options` (see :ref:`petsc_linear_solvers`),
       :cxx:`VerbosityLevel`, :cxx:`RelativeTolerance`,
       :cxx:`AbsoluteTolerance`, :cxx:`MaximumNumberOfIterations`,
       :cxx:`Preconditioner`, :cxx:`RigidBodyModes`
   * - :cxx:`BelosSolver` (parallel)
     - Trilinos
     - iterative
     - host, or GPU (Trilinos built with CUDA or HIP)
     - any
     - :cxx:`Solver`, :cxx:`Parameters` (see :ref:`trilinos_linear_solvers`),
       :cxx:`VerbosityLevel`, :cxx:`RelativeTolerance`,
       :cxx:`AbsoluteTolerance` (null only), :cxx:`MaximumNumberOfIterations`,
       :cxx:`Preconditioner`
   * - :cxx:`Amesos2Solver` (parallel)
     - Trilinos
     - direct
     - host, or GPU (Trilinos built with CUDA or HIP)
     - none
     - :cxx:`Solver`, :cxx:`Parameters` (see :ref:`trilinos_linear_solvers`)

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
   * - :cxx:`GSSmoother`
     - MFEM
     - host
     - :cxx:`Type` (:cxx:`Symmetric`, by default, :cxx:`Forward`,
       :cxx:`Backward`), :cxx:`NumberOfIterations` (1 by default)
   * - :cxx:`DSmoother`
     - MFEM
     - host
     - :cxx:`Type` (:cxx:`Jacobi`, by default, :cxx:`l1Jacobi`,
       :cxx:`LumpedJacobi`), :cxx:`Damping` (1 by default),
       :cxx:`NumberOfIterations` (1 by default); the lumped Jacobi smoother
       divides by the sums of the rows, meant for matrices whose rows have
       positive sums (mass matrices): the sums of the rows of a stiffness
       matrix vanish up to the rounding errors, and MFEM then stops the
       program or the preconditioner is meaningless
   * - :cxx:`BlockILU`
     - MFEM
     - host
     - :cxx:`BlockSize` (by default, the number of components of the
       unknowns if they are ordered by vector dimension, one otherwise),
       :cxx:`Reordering` (:cxx:`MinimumDiscardedFill`, by default,
       :cxx:`None`); without fill, the only level supported by MFEM
   * - :cxx:`UMFPack`, :cxx:`KLU`
     - SuiteSparse
     - host
     - none; exact solvers in sequential
   * - :cxx:`OrthoSolver`
     - MFEM
     - host
     - :cxx:`Preconditioner` (the wrapped preconditioner, applied to the
       vectors orthogonalized to the constant vector: operators whose
       kernel is the constant vector, such as the ones of scalar problems
       with pure Neumann boundary conditions, not the ones of mechanics,
       whose kernel holds the rigid body motions)
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
   * - :cxx:`HypreMGR` (parallel)
     - hypre
     - as hypre
     - :cxx:`CoarseComponents`, :cxx:`RelaxationType`,
       :cxx:`NumberOfRelaxationSweeps`, :cxx:`InterpolationType`,
       :cxx:`RestrictionType`, :cxx:`VerbosityLevel` (see
       :ref:`hypre_mgr`)
   * - :cxx:`NodalBlockJacobi`
     - MFEM
     - as MFEM
     - :cxx:`Backend`, :cxx:`BlockSize`, :cxx:`Factorization`,
       :cxx:`Damping` (see :ref:`nodal_block_jacobi`)
   * - :cxx:`HypreEuclid`, :cxx:`HypreParaSails` (parallel)
     - hypre
     - host
     - :cxx:`VerbosityLevel`
   * - :cxx:`AmgX`
     - AmgX
     - GPU
     - :cxx:`Configuration`, :cxx:`ConfigurationFile`, :cxx:`VerbosityLevel`,
       :cxx:`AmgXParameters`, :cxx:`Precision`, :cxx:`Multigrid`,
       :cxx:`NodalBlocks`
   * - :cxx:`GinkgoJacobi`, :cxx:`GinkgoILU`, :cxx:`GinkgoMultigrid`, ...
       (see :ref:`ginkgo_linear_solvers`)
     - Ginkgo
     - as MFEM
     - :cxx:`Configuration`, :cxx:`L1Smoother` (local preconditioners in
       parallel)
   * - :cxx:`MAGMAJacobi`, :cxx:`MAGMAILU`, :cxx:`MAGMAIC` (see
       :ref:`magma_linear_solvers`)
     - MAGMA
     - GPU
     - none
   * - :cxx:`MUMPSPreconditioner` (parallel)
     - MUMPS
     - host
     - the options of :cxx:`MUMPSSolver`: the factorization, approximate
       with :cxx:`BLRTolerance` when the fronts are large enough to be
       compressed, is applied once; with the conjugate gradient, the matrix
       shall be declared :cxx:`Symmetric` and :cxx:`PositiveDefinite`, the
       approximate factorization being otherwise a non symmetric LU
       factorization
   * - :cxx:`PETScPreconditioner` (parallel)
     - PETSc
     - host
     - :cxx:`Options` (see :ref:`petsc_linear_solvers`),
       :cxx:`RigidBodyModes`
   * - :cxx:`Ifpack2`, :cxx:`MueLu`, :cxx:`FROSch`,
       :cxx:`Amesos2Preconditioner` (parallel)
     - Trilinos
     - host, or GPU (Trilinos built with CUDA or HIP)
     - :cxx:`Type` (:cxx:`Ifpack2`), :cxx:`Solver` (:cxx:`Amesos2Preconditioner`),
       :cxx:`Parameters` (see :ref:`trilinos_linear_solvers`)
   * - :cxx:`STRUMPACKPreconditioner` (parallel)
     - STRUMPACK
     - host
     - :cxx:`Options` (see :ref:`strumpack_linear_solvers`),
       :cxx:`ReuseAnalysis`, :cxx:`VerbosityLevel`

The preconditioners :cxx:`GSSmoother`, :cxx:`DSmoother`, :cxx:`BlockILU`,
:cxx:`UMFPack` and :cxx:`KLU` are defined on the sparse matrices of MFEM:
in parallel, each process applies them to its diagonal block of the matrix
of hypre (block Jacobi method, one block per process), the number of
iterations then growing with the number of processes. They compute on the
host, the matrix being copied on the host and the vectors from and to the
device when MFEM computes on a GPU. With the execution settings
(:cxx:`--execution`), they are refused when the linear solvers are placed
on the GPU, and, when they are placed on the CPU in the :cxx:`mixed` mode,
the Krylov solvers of MFEM, which compute where MFEM computes, are refused:
they are then given to the Krylov solvers of PETSc.

The libraries which are not available (MFEM built without them, or
mfem-mgis built without AmgX, Trilinos or the sparse part of MAGMA)
register nothing.

.. _nodal_block_jacobi:

Jacobi preconditioner by nodal blocks
-------------------------------------

The :cxx:`NodalBlockJacobi` preconditioner inverts, node by node, the
block of the matrix which couples the components of the unknown at this
node (:math:`3 \times 3` blocks in three dimensional mechanics), which is
more robust than the Jacobi preconditioner when these components are
strongly coupled. The blocks are extracted from the matrix (from the
diagonal block of each process in parallel), or computed by kernels for the
tangent operator applied without assembly (see
:ref:`mfem_mgis_batch_computations_matrix_free`), then inverted or
factorized by the batched dense linear algebra of :code:`MFEM`, which
computes where :code:`MFEM` computes. Its options are:

- :cxx:`Backend`: :cxx:`Native` (kernels of :code:`MFEM`), :cxx:`GPUBLAS`
  (cuBLAS or hipBLAS) or :cxx:`MAGMA`. By default, the active backend of
  :code:`MFEM` (:code:`BatchedLinAlg::GetActiveBackend`, initially the first
  one available among MAGMA, cuBLAS or hipBLAS, and the kernels of
  :code:`MFEM`), or :cxx:`Native` if it is not available. The backends
  :cxx:`GPUBLAS` and :cxx:`MAGMA` require :code:`MFEM` computing on a GPU;
- :cxx:`BlockSize`: the size of the blocks, the number of components of
  the unknowns (by default) or one (the Jacobi preconditioner); other
  sizes are refused. With an operator applied without assembly, the size
  of the blocks must be the number of components of the unknowns: blocks
  of size one are refused when the operator is given, during the
  resolution (the :cxx:`OperatorJacobiSmoother` preconditioner is the
  Jacobi preconditioner of such an operator);
- :cxx:`Factorization`: :cxx:`Inverse` (the inverses of the blocks are
  applied by batched products, by default) or :cxx:`LU` (the LU factors are
  applied by batched solves, slightly more accurate in single precision);
- :cxx:`Damping`: a damping factor, one by default. It does not change the
  iterates of the Krylov solvers (the conjugate gradient, GMRES), only the
  ones of the stationary iterations (:cxx:`SLISolver`) and of the
  smoothers. The Krylov solvers of :code:`MFEM` stop on the preconditioned
  residual: with an :cxx:`AbsoluteTolerance`, the damping changes the
  number of iterations.

The unknowns may be ordered by nodes or by vector dimension: when they are
ordered by nodes, the vectors are gathered by nodes before the batched
products or solves, and scattered back. The backend :cxx:`MAGMA` computes
on its own queue: the device is synchronized before and after the
operations of the backends :cxx:`GPUBLAS` and :cxx:`MAGMA`. A singular
block (a node whose components are not coupled to the rest of the
problem, for instance) stops the computation in :code:`MFEM`. On a non
conforming mesh, the blocks of the operator applied without assembly are
approximated as its diagonal (see
:ref:`mfem_mgis_batch_computations_matrix_free`), and differ from the
blocks of the assembled operator.

.. _hypre_mgr:

Multigrid reduction of hypre
----------------------------

The :cxx:`HypreMGR` preconditioner applies one cycle of the multigrid
reduction of :code:`hypre` (MGR), from a null initial guess. The unknowns
are split by their components: at each level of reduction, the unknowns of
some components are kept (coarse points), and the other ones (fine points)
are relaxed (F-relaxation) and eliminated, through an interpolation and a
restriction; the coarsest level is solved by BoomerAMG. The components are
given to :code:`hypre` by markers, whatever the ordering of the unknowns.
The preconditioner of :code:`hypre` is created again, with the settings,
for each new operator. Its options are:

- :cxx:`CoarseComponents`: the components kept, a list of components for
  one level of reduction, or a list of lists of components, one per level,
  each level keeping a non empty strict subset of the components of the
  previous one. By default, one level keeps all the components but the last
  one. The unknowns must have at least two components;
- :cxx:`RelaxationType`: the F-relaxation (:code:`HYPRE_MGRSetRelaxType`):
  :cxx:`0` (Jacobi, by default), :cxx:`2` (algebraic multigrid on the fine
  points), :cxx:`1` (V-cycle), :cxx:`3` to :cxx:`8`, :cxx:`13`, :cxx:`14`
  and :cxx:`18` (relaxations of BoomerAMG), :cxx:`9`, :cxx:`99` and
  :cxx:`199` (dense direct solvers of the block of the fine points, for
  small problems only);
- :cxx:`NumberOfRelaxationSweeps`: the number of sweeps of the
  F-relaxation;
- :cxx:`InterpolationType`: :cxx:`0`, :cxx:`1`, :cxx:`2` (by default),
  :cxx:`3`, :cxx:`4`, :cxx:`12` or :cxx:`14`, see
  :code:`HYPRE_MGRSetInterpType` (the interpolation :cxx:`13` is refused:
  it gives wrong results);
- :cxx:`RestrictionType`: :cxx:`0` (injection, by default), :cxx:`1` to
  :cxx:`5` or :cxx:`12` to :cxx:`15`, see :code:`HYPRE_MGRSetRestrictType`;
- :cxx:`VerbosityLevel`: the print level of the multigrid reduction,
  :cxx:`0` to :cxx:`7` (information on the setup, the resolution and the
  parameters; the other bits of the print level of :code:`hypre` write
  files).

The default cycle (Jacobi interpolation, injection, F-relaxation before
the coarse level only) is not symmetric: the preconditioner is meant for
the GMRES solvers rather than for the conjugate gradients.

The block interpolations and restrictions (:cxx:`12` to :cxx:`15`, and the
block F-relaxation of the Jacobi relaxation with the interpolation
:cxx:`12`) group consecutive fine unknowns, which are the components
eliminated at a node only when the unknowns are ordered by vector
dimension: with the ordering by nodes, they are refused when a level
eliminates several components (with one component, they are pointwise).
With :code:`hypre` 3.2.0, the restriction :cxx:`12` stops the program when
a level eliminates several components, and the restrictions :cxx:`14` and
:cxx:`15` give invalid values when a level eliminates one component: these
cases are refused.

The multigrid reduction is meant for coupled problems whose fine points are
easy to relax. On the beam in linear elasticity of the test
:code:`HypreMGRTest`, the Jacobi F-relaxation of :code:`hypre` makes it
slower than the Jacobi preconditioner, while, with the F-relaxation by
algebraic multigrid (:cxx:`RelaxationType` set to :cxx:`2`), the GMRES
solver of :code:`MFEM` converges in about 35 iterations on the host and 50
on a GPU, against 150 to 300 with the Jacobi preconditioner (relative
tolerance of :math:`10^{-12}`).

When :code:`hypre` computes on a GPU, the F-relaxations :cxx:`0`,
:cxx:`2`, :cxx:`7`, :cxx:`9`, :cxx:`18`, :cxx:`99` and :cxx:`199`, the
interpolations :cxx:`0`, :cxx:`1`, :cxx:`2` and :cxx:`12`, and the
restrictions :cxx:`0`, :cxx:`1`, :cxx:`2` and :cxx:`12` are accepted; the
other values are refused: :code:`hypre` 3.2.0 applies the relaxations,
interpolations and restrictions of the host to the matrices on the GPU and
stops the program, and the V-cycle relaxation frees memory which is not
its own. The interpolation :cxx:`12` is also refused when a level
eliminates several components (invalid memory accesses). :code:`hypre`
replaces the Jacobi relaxation :cxx:`0` by its relaxation of the GPU,
:cxx:`7`, except with the interpolation :cxx:`12`, which it then replaces
by :cxx:`2` at the levels eliminating one component: mfem-mgis then gives
the relaxation :cxx:`7` to these levels, which would otherwise apply the
relaxation of the host.

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

The solvers and the preconditioners of PETSc compute on the host, even
when :code:`MFEM` computes on a GPU. The matrices of hypre are copied in
the AIJ format of the host by mfem-mgis: the conversion of :code:`MFEM`
goes through the :code:`MATHYPRE` format of PETSc, whose creation
initializes the GPU when PETSc and hypre are built for it (a context per
process, about 260 MiB of the GPU and 0.4 second each with CUDA 13), even
for a computation on the host. When :code:`MFEM` computes on a GPU, the
vectors are copied from and to the vectors of :code:`MFEM`, the
preconditioners of mfem-mgis given to the Krylov solvers of PETSc being
applied on the vectors of :code:`MFEM`. With :code:`MFEM`
4.10 and :code:`PETSc` 3.26 built with CUDA, the vectors of :code:`MFEM`
placed in the vectors of PETSc (whose internal offload mask :code:`MFEM`
sets) were seen as null by the Krylov solvers of PETSc and by the shell
preconditioners, and the vectors of PETSc of the GPU made
:code:`PetscFinalize` abort the program when freeing the memory of PETSc on
the GPU.

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

When the number of components of the unknowns is the space dimension (2 or
3) and the unknowns are ordered by vector dimension, the converted matrices
receive the number of components as block size and the rigid body modes
of the finite element space as near null space (computed from the
coordinates of its nodes): the algebraic multigrid of PETSc (GAMG) then
aggregates the nodes and builds its coarse spaces from these modes (15
conjugate gradient iterations instead of 51 on the beam of the tests). The
option :cxx:`RigidBodyModes` (true by default) of
:cxx:`PETScLinearSolver` and :cxx:`PETScPreconditioner` disables them;
when :cxx:`PETScPreconditioner` preconditions :cxx:`PETScLinearSolver`,
the two share the same object of PETSc and the option of the
preconditioner applies. They are not given with the unknowns ordered by
nodes (the default of mfem-mgis, see the parameter
:cxx:`FiniteElementSpaceOrdering`), the block size of PETSc grouping
consecutive rows, nor to a matrix of another size than the finite element
space. The block size is also seen by the other preconditioners of PETSc
(the point block Jacobi method uses blocks of this size, the field split
splits by default by components). The nonlinear solver :cxx:`SNES` gives
them likewise to the tangent operator converted for the linear solvers of
PETSc (option :cxx:`RigidBodyModes` of :cxx:`SNES`).

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

The external packages which PETSc is built with are available through
the same options, for instance the sparse direct solvers MUMPS, PaStiX,
SuperLU_DIST or STRUMPACK (:cxx:`ksp_type` :cxx:`preonly`,
:cxx:`pc_type` :cxx:`lu` and :cxx:`pc_factor_mat_solver_type`
:cxx:`pastix`), or the Krylov solvers of HPDDM, among which GCRODR
recycles a subspace between the resolutions (:cxx:`ksp_type` :cxx:`hpddm`,
:cxx:`ksp_hpddm_type` :cxx:`gcrodr`, :cxx:`ksp_hpddm_recycle`).

PETSc is initialized when the first linear solver or preconditioner of
PETSc is built, which shall thus be built by all the processes. It is
then finalized by the :cxx:`finalize` function, unless it was initialized
by the program itself. The command line is not given to PETSc: the
options of the solvers are only given by their :cxx:`Options` parameter.
The non linear solver of PETSc, selected by the :cxx:`--use-petsc`
command line option, is not related to these solvers.

.. _trilinos_linear_solvers:

Linear solvers and preconditioners of Trilinos
==============================================

When mfem-mgis is built with Trilinos (CMake option
:code:`enable-trilinos`, see :doc:`../installation_guide/installation_guide`),
the following entries are available in parallel, in double precision:

- :cxx:`BelosSolver`: a Krylov solver of Belos, selected by the
  :cxx:`Solver` option (:cxx:`GMRES` by default, :cxx:`CG`,
  :cxx:`BICGSTAB`, :cxx:`GCRODR`, :cxx:`MINRES`, :cxx:`TFQMR` and the other
  names of the factory of solvers of Belos), preconditioned on the right;
- :cxx:`Amesos2Solver`: a direct solver of Amesos2, selected by the
  :cxx:`Solver` option (:cxx:`KLU2` by default, and the other solvers
  which Trilinos is built with, such as :cxx:`Basker` or :cxx:`Tacho`);
- :cxx:`Ifpack2`: a preconditioner of Ifpack2, selected by the :cxx:`Type`
  option (:cxx:`RELAXATION` by default, a Jacobi relaxation, :cxx:`CHEBYSHEV`,
  :cxx:`RILUK`, :cxx:`ILUT`, :cxx:`SCHWARZ`...);
- :cxx:`MueLu`: one V-cycle of the smoothed aggregation algebraic multigrid
  of MueLu;
- :cxx:`FROSch`: the two level overlapping Schwarz preconditioner of FROSch
  (:code:`TwoLevelPreconditioner`: algebraic overlap, one layer by default,
  coarse space of the interface partition of unity, the subdomains and the
  coarse problem being solved by KLU2 by default), additive (the values of
  the overlaps are added, :code:`Combine Values in Overlap` being
  :code:`Full`, unless given in the sublist
  :code:`AlgebraicOverlappingOperator`), hence symmetric, unlike the
  restricted Schwarz method used by FROSch by default;
- :cxx:`Amesos2Preconditioner`: the factorization of a direct solver of
  Amesos2, applied once.

The :cxx:`Parameters` option of each entry is the list of parameters of
the solver of Trilinos (a :code:`Teuchos::ParameterList`): booleans,
integers, reals and strings, a list of parameters being a sublist. For
Amesos2, it is the content of the list :code:`Amesos2` of Amesos2 (the
parameters common to the solvers, and a sublist per solver, such as
:code:`KLU2`), which Amesos2 validates when the matrix is set. The names
of the first level of the parameters of Belos (and those of the outer
solver of the hybrid GMRES) are checked when the solver is created, which
Belos only does when they
are given again (a change of the tolerance): a misspelt parameter is
refused, except by LSQR and the native solvers of Tpetra
(:cxx:`TPETRA GMRES`...), which ignore the unknown parameters. The
parameters of Ifpack2 are given, when the solver is created, to a
preconditioner of the chosen type built without matrix: the relaxation
(:cxx:`RELAXATION`) refuses a misspelt parameter, the other preconditioners
ignore the unknown parameters at this stage. The block relaxations check
their parameters at the setup, and the additive Schwarz method the
parameters of its inner preconditioner; :cxx:`AMESOS2` reads its
parameters in a sublist :code:`Amesos2`, and :cxx:`HIPTMAIR`, which needs
the matrices of an auxiliary space, can't be used. The names of
the first level of the parameters of MueLu are checked when the solver is
created in its simplified syntax, the default one; the other syntaxes
(factories given in a sublist :code:`Hierarchy`, facade :code:`MueLu
preconditioner`, parameters of ML with :code:`parameterlist: syntax` set
to :code:`ml`) are checked by MueLu only (MueLu 17.2.1 translates the
parameters of ML to an option of the aggregation that its version on
Kokkos does not support). FROSch does not check its parameters:

.. code-block:: cpp

   problem.setLinearSolver(
       ctx, "BelosSolver",
       {{"Solver", "CG"},
        {"RelativeTolerance", 1e-10},
        {"Preconditioner",
         Parameters{{"Name", "MueLu"},
                    {"Options",
                     Parameters{{"Parameters",
                                 Parameters{{"smoother: type",
                                             "CHEBYSHEV"}}}}}}}}) |
       or_die;

The options of mfem-mgis set the tolerance (:code:`Convergence
Tolerance`), the maximum number of iterations (:code:`Maximum Iterations`)
and the verbosity (:code:`Verbosity`, and :code:`Output Frequency` for the
details of each iteration) of Belos, when the solver accepts these
parameters, the parameters of Belos overriding them; a relative tolerance
given to a solver without :code:`Convergence Tolerance` (LSQR, which has
its own tolerances) is refused. The solvers of Belos only have a relative
tolerance: by default, the norm of the residual is compared to the one of
the initial residual. The relative tolerance can be changed after the
creation of the solver (it is then given to the adaptive tolerance of the
Newton solver).
The GMRES solvers keep 50 Krylov vectors (:code:`Num Blocks`, 300 by
default in most GMRES solvers of Belos), and are restarted as many times
as the maximum number of iterations allows (:code:`Maximum Restarts`),
unless these parameters are given. The hybrid GMRES of Belos
(:cxx:`HYBRID BLOCK GMRES`, or its aliases :cxx:`GMRESPOLY` and
:cxx:`SEED GMRES`, a GMRES polynomial used as preconditioner) iterates
with an outer solver (:code:`Outer Solver`,
:code:`Block GMRES` by default, among the solvers of the generic factory
of Belos), whose parameters (the sublist :code:`Outer Solver Params`)
receive the tolerance, the maximum number of iterations, the verbosity
and the parameters of the restarted GMRES. With an empty name of outer
solver, the polynomial is applied alone, without tolerance nor maximum
number of iterations; with a null degree (:code:`Maximum Degree`), the
outer solver is applied alone, and is built by the factory of Belos for
Tpetra. The solvers of Belos are built again when the matrix changes
(the hybrid GMRES then builds its polynomial again, and PCPG, which keeps
a seed space only valid for the same matrix, forgets it), except GCRODR,
which keeps its recycled subspace for the matrices of the same rows. A
failure of a solver of
Belos (for instance the loss of the positive definiteness detected by the
conjugate gradient) is reported as a resolution which did not converge,
its message being printed.

The matrices of hypre are converted to the sparse matrices of Tpetra each
time they are set (the preconditioners and the factorizations being then
built again): the operators shall be assembled. When the unknowns have
several components and are ordered by nodes, the unknowns of each process
are seen by Trilinos ordered by vector dimension, the components of a
node being contiguous, as MueLu expects: the rows and the columns of each
process are permuted, and the vectors accordingly. MueLu aggregates the
nodes (its parameter :code:`number of equations` is the number of
components, unless it is given) and, when the number of components is the
space dimension (2 or 3), its near null space is made of the rigid body
modes of the finite element space, computed from the coordinates of its
nodes; the rigid body modes are not used when another number of equations
is given. The parameter :code:`number of equations` is only set in the
simplified syntax of MueLu; the rigid body modes and the coordinates of
the nodes (user data :code:`Nullspace` and :code:`Coordinates`) are given
in all the syntaxes.

MueLu solves directly, without multigrid, the systems smaller than its
parameter :code:`coarse: max size` (2000 unknowns by default). The
coordinates of the nodes are also given to MueLu (user data
:code:`Coordinates`), for its options based on them, unless another number
of equations is given. Unless they are given, mfem-mgis sets, in the
simplified syntax of MueLu, the parameters :code:`number of equations`
(see above), :code:`aggregation: deterministic` (reproducible aggregates,
MueLu aggregating the nodes in a random order by default),
:code:`verbosity` (:code:`none`, a silent setup) and, on CUDA, in the
sublist :code:`matrixmatrix: kernel params`, :code:`cuda: algorithm`
(:code:`SPGEMM_KK`: the products of sparse matrices of Kokkos Kernels,
those of cuSPARSE, the default of Tpetra on CUDA, running out of memory on
a GPU of 16 GB from a few hundred thousand unknowns). On a GPU, the
default smoother of MueLu (symmetric Gauss-Seidel) is slow: on a beam of
421,443 unknowns, a resolution took 3.8 seconds, against 0.071 seconds
with two sweeps of a Jacobi relaxation (:code:`smoother: type` set to
:code:`RELAXATION`, :code:`smoother: params` to :code:`{relaxation: type =
Jacobi, relaxation: sweeps = 2, relaxation: damping factor = 0.6}`), for
22 to 23 iterations instead of 20; the setup of MueLu (0.84 seconds) then
dominates.

The hierarchy of MueLu is built again for each operator, unless its
parameter :code:`reuse: type` is given (:code:`RP`: the prolongators and
the restrictions of the first operator are kept, the coarse operators
being computed again; :code:`tP`: the tentative prolongators are
kept...): the hierarchy is then updated by
:code:`MueLu::ReuseTpetraPreconditioner` when the distribution of the rows
of the new operator is the same, and built again otherwise. The smoothers
are kept and set up again for the new operator, except the estimate of
the largest eigenvalue of the Chebyshev smoother, kept from the first
operator: with a stiffening operator, the Jacobi smoother is safer. On
the beam of 421,443 unknowns, on the GPU, with the Jacobi smoother above,
:code:`RP` sets up a second operator in 0.17 seconds instead of 0.82.
Keeping the whole hierarchy (:code:`full`) keeps the coarse operators of
the first operator, which suits only operators that do not change.

FROSch is given the number of components as number of degrees of freedom
per node, and the space dimension (2 or 3). Its coarse space is built from
the null space of the linear elasticity (computed by FROSch from the
coordinates of the nodes) when the number of components is the space
dimension, from the constants otherwise, unless its parameter
:code:`Null Space Type` is given. The rows of the Dirichlet boundary
conditions are not given to FROSch (an empty list is given, FROSch
building its own list otherwise): with FROSch 17.2.1, the list of the rows
whose entries are null except the diagonal one (built by
:code:`FindOneEntryOnlyRowsGlobal` with a tolerance, the matrices of
:code:`MFEM` keeping the eliminated entries) made its initialization hang
on several processes with three unknowns per node, and did not change the
number of iterations on the Laplace equation; its parameter
:code:`DirichletBoundary Numerical Tolerance` is refused. The setup of
FROSch is less verbose (:code:`Verbosity` set to :code:`None` in the
sublists of the overlapping operator and of the coarse operator, unless
given), FROSch still printing some tables on the first process. With
large subdomains, its default solver (KLU2) may be very slow: a solver of
Amesos2 suited to the finite element matrices (Tacho, for instance) can be
given in its sublists. FROSch does not support empty subdomains (it
crashed when a process owned no unknown): the setup fails, on all the
processes, when a process owns no unknown.

A preconditioner of Trilinos given to :cxx:`BelosSolver` shares the matrix
of the solver and is applied to the vectors of Tpetra; the preconditioners
of the other libraries are applied to vectors of :code:`MFEM`, the vectors
of Tpetra being copied, on the GPU when :code:`MFEM` and Tpetra compute
there. The preconditioners of Trilinos may also be used by the Krylov
solvers of :code:`MFEM`.

The solvers of Trilinos compute where :code:`MFEM` computes. When
:code:`MFEM` computes on the GPU and Trilinos is built for it, Tpetra
computes with its default node, on the GPU. Otherwise, Tpetra computes on
the host with its serial node (each process of :code:`MPI` computing on
one core); when Tpetra is not instantiated for its serial node, the
solvers compute on
the GPU even when :code:`MFEM` computes on the host. FROSch and the
factorization of Tacho always compute with the default node of Tpetra, on
the GPU when Trilinos is built for one (FROSch 17.2.1 runs some loops in
the default execution space of :code:`Kokkos`, as Tacho, whatever the
node). :cxx:`BelosSolver` computes at the same place as its
preconditioner of Trilinos. Some options of MueLu may not be available on
the host in this case (the distance Laplacian with weights of MueLu 17.2.1
allocates its weights in the default memory space). The matrices are
converted where Tpetra computes, by kernels of :code:`Kokkos`, from the
arrays of hypre on the GPU when hypre computes there. The vectors of
:code:`MFEM` are copied on the GPU when :code:`MFEM` and Tpetra compute
on the GPU, and from and to the host otherwise. Trilinos requires an
initialized :code:`Kokkos`,
which is initialized, if needed, when the first entry of Trilinos is
built, on the GPU of the device of :code:`MFEM`, and finalized by
:cxx:`finalize` before :code:`MPI` (see
:doc:`../installation_guide/installation_guide`): the solvers of Trilinos,
whose data are allocated by :code:`Kokkos`, shall be destroyed before.
When Trilinos is built without GPU and its default node is the OpenMP
node, the computations use the number of threads given by
:code:`OMP_NUM_THREADS`, which, with several processes per node, shall be
set accordingly (one thread per process if the processes occupy all the
cores).

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

.. _mkl_linear_solvers:

Sparse direct solvers of oneMKL
===============================

When :code:`MFEM` is built with them (:code:`MFEM_USE_MKL_PARDISO`,
:code:`MFEM_USE_MKL_CPARDISO`), the PARDISO solver of oneMKL
(:cxx:`MKLPardisoSolver`, sequential) and its Cluster Sparse Solver
(:cxx:`MKLClusterSparseSolver`, parallel, the rows of the matrix being
distributed as the ones of the matrix of hypre) are available. They are
called directly by mfem-mgis, not through the classes of :code:`MFEM`,
whose resolutions on the GPU and checks of the errors are not reliable.
The matrix is copied on the host in the compressed row format of oneMKL
(the upper triangle for the symmetric types, :cxx:`Symmetric`, with
:cxx:`PositiveDefinite` for a positive definite matrix), the vectors are
copied on the host. The analysis (reordering and symbolic factorization)
is reused for the next operators of the same pattern with
:cxx:`ReuseAnalysis` (with the scaling and the matching computed from the
values of the first matrix, which suit the next ones as long as their
values are close), and computed again otherwise. For the symmetric types,
only the upper triangle of the matrix is read. A failure of the analysis or
of the factorization stops the resolution on all the processes. oneMKL
perturbs the pivots which are too small, and factorizes a singular matrix
without error: as for cuDSS, a resolution is reported as converged only if
its residual is lower than :math:`10^{-3}` times the right hand side.
:cxx:`VerbosityLevel` prints the statistics of oneMKL.

.. _magma_linear_solvers:

Sparse solvers of MAGMA
=======================

When :code:`MFEM` is built with MAGMA (:code:`MFEM_USE_MAGMA`, used by
:code:`MFEM` for its batched dense algebra) in double precision, and the
library of the sparse part of MAGMA (:code:`magma_sparse`) and its header
(:code:`magmasparse.h`) are found next to the library of MAGMA used by
:code:`MFEM` (or given by the CMake variables
:code:`MFEM_MGIS_MAGMA_SPARSE_LIBRARY` and
:code:`MFEM_MGIS_MAGMA_SPARSE_INCLUDE_DIR`), mfem-mgis registers the Krylov
solvers and the preconditioners of MAGMA, which compute on one GPU (an
error is reported when generating them without GPU):

- the Krylov solvers :cxx:`MAGMACG`, :cxx:`MAGMABiCGSTAB`,
  :cxx:`MAGMAGMRES` (restarted flexible GMRES, :cxx:`KDim` giving the
  dimension of the Krylov space, 50 by default), :cxx:`MAGMAIDR`
  (:cxx:`ShadowSpaceDimension` giving the dimension :math:`s` of the shadow
  space of IDR(s), 1 by default, bounded by MAGMA by the size of the
  system; MAGMA reads 50 as its default dimension, and this value is
  refused), :cxx:`MAGMACGS`, :cxx:`MAGMAQMR` and :cxx:`MAGMATFQMR`, in
  sequential. Their tolerances are relative to the norm of the right hand
  side (MAGMA's defaults: :math:`10^{-10}` relative, :math:`10^{-16}`
  absolute, 1000 iterations). A resolution is reported as converged if
  MAGMA reports a success and if the residual recomputed by MAGMA at the
  end is lower than ten times the tolerance (the preconditioned TFQMR of
  MAGMA 2.10 reports its breakdowns as successes). They only accept the
  preconditioners of MAGMA, which they apply themselves. GMRES counts the
  products by the matrix, one more at each restart. QMR, without
  look-ahead, may break down (reported as not converged), and refuses the
  incomplete Cholesky factorization, whose transpose MAGMA 2.10 applies
  wrongly: for a symmetric matrix, :cxx:`MAGMAILU` gives the same
  factorization;
- the preconditioners :cxx:`MAGMAJacobi`, :cxx:`MAGMAILU` and
  :cxx:`MAGMAIC` (incomplete LU and Cholesky factorizations without fill-in
  of cuSPARSE or hipSPARSE, :cxx:`MAGMAIC` reading the lower triangle of
  the matrix), without options, given to the Krylov solvers of MAGMA or to
  any other solver, for instance the Krylov solvers of MFEM. In parallel,
  each process applies them to its diagonal block of the matrix of hypre,
  as :cxx:`BlockILU`. A null or missing diagonal entry is refused (on all
  the processes in parallel).

The matrix is copied on the host, its columns being sorted, then on the
GPU for each operator, and the preconditioner set up at the same time. The
vectors of MFEM are given to MAGMA when MFEM computes on the GPU and they
use the memory of its device (:cxx:`UseDevice`); otherwise they are copied
at each resolution or application of a preconditioner.

.. code-block:: cpp

   problem.setLinearSolver(
       ctx, "MAGMACG",
       {{"RelativeTolerance", 1e-10},
        {"Preconditioner", Parameters{{"Name", "MAGMAIC"}}}}) |
       or_die;

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
