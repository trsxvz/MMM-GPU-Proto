======================
Glossary Of Parameters
======================

``MFEM-MGIS`` offers a simplified API based on encapsulating
``Parameter`` objects in a ``Parameters`` class, which is then given to
various functions to implement the problem to be solved. It is important
to note that if a Parameter is not listed in the glossary, the code will
stop and an error message will specify the available options. In this
section we describe the options available for each part of the code,
both for setting the parameters of the linear solver and for defining
the problem to be solved.

.. contents::
    :depth: 3
    :local:

.. note:: 

  The list of default parameters is adapted according to compilation
  flags. For example, if ``MFEM_USE_MUMPS`` is not active, ``MUMPS`` will
  not appear in the list. So if one of the options listed below doesn't
  appear, ``MFEM-MGIS`` may not have been compiled with all the options.

Setup Your Linear Solver
========================

In this section, we'll distinguish between two types of solver: direct
solvers, which are not always listed in the parallel version of the
problems, and Krylov solvers, which are available in both serial and
parallel versions (except for ``hypre`` solvers).

To choose your solver, pass a Parameters object to the `setLinearSolver`
function.

.. note::

  With the option ``--execution``, the linear solvers and the
  preconditioners are checked against the execution settings: a solver
  computing elsewhere than the settings place the linear solvers is
  refused, with a message giving the options which would accept it (see
  :ref:`mfem_mgis_batch_computations_execution`).

Direct Solver
-------------

MUMPS Solver
^^^^^^^^^^^^

- Solver name: ``MUMPSSolver``
- Dependencies:
 
  - ``MFEM_USE_MUMPS``
  - parallel

- Parameters:

+--------------------------+------------------------------------------------------+
| Key                      | Description                                          |
+==========================+======================================================+
| ``Symmetric``            | Specify if the matrix is symmetric. [bool]           |
+--------------------------+------------------------------------------------------+
| ``PositiveDefinite``     | Specify if the matrix is positive definite. [bool]   |
+--------------------------+------------------------------------------------------+
| ``ReuseAnalysis``        | Reuse the analysis (ordering and symbolic            |
|                          | factorization) of the previous matrix when the new   |
|                          | one has the same pattern on all processes. [bool,    |
|                          | false by default]                                    |
+--------------------------+------------------------------------------------------+
| ``ReorderingStrategy``   | Reordering of the matrix: ``Automatic`` (by          |
|                          | default), ``AMD``, ``AMF``, ``PORD``, ``METIS``,     |
|                          | ``ParMETIS``, ``Scotch``, ``PTScotch`` (MUMPS must   |
|                          | be built with the corresponding library). [string]   |
+--------------------------+------------------------------------------------------+
| ``BLRTolerance``         | Tolerance of the block low-rank factorization, which |
|                          | is approximate and whose result is not checked       |
|                          | (zero, exact factorization, by default). [double]    |
+--------------------------+------------------------------------------------------+
| ``VerbosityLevel``       | Print level of MUMPS (``ICNTL(4)``): 0 nothing, 1    |
|                          | errors (by default), 2 warnings and statistics, up   |
|                          | to 4. [int]                                          |
+--------------------------+------------------------------------------------------+

The matrices of the successive iterations of a non linear solver usually
have the same pattern: with ``ReuseAnalysis``, only the first one is
analysed. The pattern of each new matrix is compared to the one of the
last analysed matrix, and a new analysis is made if it changed on one of
the processes. The matrix being distributed, MUMPS does not use the
maximum transversal permutation and computes the scaling at each
factorization: its analysis only depends on the pattern of the matrix, and
reusing it does not change the factorization. The function
``getNumberOfAnalyses`` gives the number of analyses made.

**Example:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "MUMPSSolver", {{"Symmetric", true}}) | or_die;

cuDSS Solver
^^^^^^^^^^^^

Direct solver of NVIDIA on the GPU, through the ``CuDSSSolver`` class of
MFEM.

- Solver name: ``CuDSSSolver``
- Dependencies:

  - ``MFEM_USE_CUDSS``
  - parallel
  - a device on a CUDA GPU (option ``--execution gpu``, or
    ``--device cuda``, of ``mfem_mgis::initialize``): the solver is refused
    otherwise; with ``--execution mixed``, the linear solvers must be
    placed on the GPU
  - when hypre computes on the host (option ``--hypre-on-host``), a GPU
    reading the pageable memory of the host (heterogeneous memory
    management): the solver is refused otherwise
  - on several processes, a CUDA-aware MPI and the communication layer of
    cuDSS (variable ``CUDSS_COMM_LIB``, or the library given to MFEM at its
    configuration)

- Parameters:

+----------------------+----------------------------------------------------+
| Key                  | Description                                        |
+======================+====================================================+
| ``Symmetric``        | Specify if the matrix is symmetric: only its lower |
|                      | triangle is then read. [bool]                      |
+----------------------+----------------------------------------------------+
| ``PositiveDefinite`` | Specify if the matrix is positive definite. [bool] |
+----------------------+----------------------------------------------------+
| ``ReuseAnalysis``    | Reuse the analysis (ordering and symbolic          |
|                      | factorization) of the previous matrix when the new |
|                      | one has the same pattern on all processes. [bool,  |
|                      | false by default]                                  |
+----------------------+----------------------------------------------------+

As for MUMPS, the pattern of each new matrix is compared to the one of the
last analysed matrix, and the function ``getNumberOfAnalyses`` gives the
number of analyses made. The operators must be ``HypreParMatrix`` objects,
which exist as long as the solver is used.

cuDSS does not report the failures of the factorization (null pivots,
matrix declared positive definite which is not): each solution is checked
by its relative residual, which must be lower than ``1e-3``, and a failure
is reported by ``hasConverged``, so that the non linear solvers reduce the
time step. The other errors of cuDSS stop the program.

With cuDSS 0.7.1, the distributed analysis (two processes or more) fails
for some matrices: when some rows only have a diagonal entry, for
tridiagonal matrices coupling the processes, or for some uneven
distributions of the rows (observed with several processes on one GPU).
The tangent operators of ``mfem-mgis`` keep the entries of the rows of the
imposed degrees of freedom; on one process, no failure was observed.

**Example:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "CuDSSSolver",
                          {{"Symmetric", true}, {"ReuseAnalysis", true}})
      | or_die;

UMFPack Solver
^^^^^^^^^^^^^^

- Solver name: ``UMFPackSolver``
- Dependencies: 

  - ``MFEM_USE_SUITESPARSE``
  - not parallel

- Parameters: None

**Example:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "UMFPackSolver", {}) | or_die;

KLU Solver
^^^^^^^^^^

- Solver name: ``KLUSolver``
- Dependencies:

  - ``MFEM_USE_SUITESPARSE``
  - not parallel

- Parameters: None

A singular matrix makes the resolution fail.

**Example:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "KLUSolver", {}) | or_die;

.. note::

  Direct solvers do not need any tolerance parameter or preconditioner parameter.

.. note::

  We advise you not to use the direct solvers with the parallel option for large systems,
  this advice is mainly due to the fact that these solvers don't scale
  in memory on large problems (e.g. meshes with more than a hundred
  thousand points). Nevertheless, these solvers are maintained within the
  framework of small computations achievable on a standard laptop.

Krylov Solver
-------------

In this section, we'll start with the solvers present in ``mfem`` and
then describe the ``hypre`` solvers available only in the parallel
version.

Conjugate Gradient
^^^^^^^^^^^^^^^^^^

- Solver name: ``CGSolver``
- Dependencies: None
- Parameters:


+-------------------------------+--------------------------------------------------+
| Key                           | Description                                      |
+===============================+==================================================+
| ``Preconditioner``            | Define your preconditioner. [Parameters]         |
+-------------------------------+--------------------------------------------------+
| ``VerbosityLevel``            | Define the verbosity of the solver. [int: 0,1]   |
+-------------------------------+--------------------------------------------------+
| ``RelativeTolerance``         | Define the relative tolerance, >0. [double]      |
+-------------------------------+--------------------------------------------------+
| ``AbsoluteTolerance``         | Define the absolute tolerance, >= 0. [double]    |
+-------------------------------+--------------------------------------------------+
| ``MaximumNumberOfIterations`` | Maximum number of Krylov iterations, >= 1. [int] |
+-------------------------------+--------------------------------------------------+

**Example:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "CGSolver",
                          {{"VerbosityLevel", 1},
                           {"AbsoluteTolerance", 1e-12},
                           {"RelativeTolerance", 1e-12},
                           {"MaximumNumberOfIterations", 5000}}) |
      or_die;


Generalized Minimal Residual (GMRES)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``GMRESSolver``
- Dependencies: None
- Parameters:

+-------------------------------+--------------------------------------------------+
| Key                           | Description                                      |
+===============================+==================================================+
| ``Preconditioner``            | Define your preconditioner. [Parameters]         |
+-------------------------------+--------------------------------------------------+
| ``VerbosityLevel``            | Define the verbosity of the solver. [int: 0,1]   |
+-------------------------------+--------------------------------------------------+
| ``RelativeTolerance``         | Define the relative tolerance, >0. [double]      |
+-------------------------------+--------------------------------------------------+
| ``AbsoluteTolerance``         | Define the absolute tolerance, >= 0. [double]    |
+-------------------------------+--------------------------------------------------+
| ``MaximumNumberOfIterations`` | Maximum number of Krylov iterations, >= 1. [int] |
+-------------------------------+--------------------------------------------------+
| ``KDim``                      | Define the Krylov dimension, >= 1. [int]         |
+-------------------------------+--------------------------------------------------+

**Example:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "GMRESSolver",
                          {{"VerbosityLevel", 1},
                           {"AbsoluteTolerance", 1e-12},
                           {"RelativeTolerance", 1e-12},
                           {"MaximumNumberOfIterations", 100000}}) |
      or_die;


Flexible Generalized Minimal Residual (FGMRES)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``FGMRESSolver``
- Dependencies: None
- Parameters: the ones of ``GMRESSolver``. Unlike GMRES, FGMRES allows a
  preconditioner which changes from one iteration to the next, for instance
  with the ``NewtonK0`` method of the fixed point solver, when the linear
  solver of the elastic operator is an iterative solver with a tolerance.

**Example:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "FGMRESSolver",
                          {{"VerbosityLevel", 1},
                           {"RelativeTolerance", 1e-12},
                           {"MaximumNumberOfIterations", 1000},
                           {"KDim", 100}}) |
      or_die;


Biconjugate Gradient Stabilized (BiCGSTAB)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``BiCGSTABSolver``
- Dependencies: None
- Parameters:

+-------------------------------+--------------------------------------------------+
| Key                           | Description                                      |
+===============================+==================================================+
| ``Preconditioner``            | Define your preconditioner. [Parameters]         |
+-------------------------------+--------------------------------------------------+
| ``VerbosityLevel``            | Define the verbosity of the solver. [int: 0,1]   |
+-------------------------------+--------------------------------------------------+
| ``RelativeTolerance``         | Define the relative tolerance, >0. [double]      |
+-------------------------------+--------------------------------------------------+
| ``AbsoluteTolerance``         | Define the absolute tolerance, >= 0. [double]    |
+-------------------------------+--------------------------------------------------+
| ``MaximumNumberOfIterations`` | Maximum number of Krylov iterations, >= 1. [int] |
+-------------------------------+--------------------------------------------------+

**Example:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "BiCGSTABSolver",
                          {{"VerbosityLevel", 1},
                           {"AbsoluteTolerance", 1e-12},
                           {"RelativeTolerance", 1e-12},
                           {"MaximumNumberOfIterations", 1000}}) |
      or_die;


Minimal Residual (MINRES)
^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``MINRESSolver``
- Dependencies: None
- Parameters:

+-------------------------------+--------------------------------------------------+
| Key                           | Description                                      |
+===============================+==================================================+
| ``Preconditioner``            | Define your preconditioner. [Parameters]         |
+-------------------------------+--------------------------------------------------+
| ``VerbosityLevel``            | Define the verbosity of the solver. [int: 0,1]   |
+-------------------------------+--------------------------------------------------+
| ``RelativeTolerance``         | Define the relative tolerance, >0. [double]      |
+-------------------------------+--------------------------------------------------+
| ``AbsoluteTolerance``         | Define the absolute tolerance, >= 0. [double]    |
+-------------------------------+--------------------------------------------------+
| ``MaximumNumberOfIterations`` | Maximum number of Krylov iterations, >= 1. [int] |
+-------------------------------+--------------------------------------------------+

**Example:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "MINRESSolver",
                          {{"VerbosityLevel", 1},
                           {"AbsoluteTolerance", 1e-12},
                           {"RelativeTolerance", 1e-12},
                           {"MaximumNumberOfIterations", 1000}}) |
      or_die;

Stationary Linear Iteration (SLI)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``SLISolver``
- Dependencies: None
- Parameters:

+-------------------------------+--------------------------------------------------+
| Key                           | Description                                      |
+===============================+==================================================+
| ``Preconditioner``            | Define your preconditioner. [Parameters]         |
+-------------------------------+--------------------------------------------------+
| ``VerbosityLevel``            | Define the verbosity of the solver. [int: 0,1]   |
+-------------------------------+--------------------------------------------------+
| ``RelativeTolerance``         | Define the relative tolerance, >0. [double]      |
+-------------------------------+--------------------------------------------------+
| ``AbsoluteTolerance``         | Define the absolute tolerance, >= 0. [double]    |
+-------------------------------+--------------------------------------------------+
| ``MaximumNumberOfIterations`` | Maximum number of Krylov iterations, >= 1. [int] |
+-------------------------------+--------------------------------------------------+

**Example:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "SLISolver",
                          {{"VerbosityLevel", 1},
                           {"AbsoluteTolerance", 1e-12},
                           {"RelativeTolerance", 1e-12},
                           {"MaximumNumberOfIterations", 1000}}) |
      or_die;

Preconditioned Conjugate Gradient (HyprePCG)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``HyprePCG``
- Dependencies:

  - ``hypre``
  - parallel

- Parameters:

+-------------------------------+--------------------------------------------------+
| Key                           | Description                                      |
+===============================+==================================================+
| ``Preconditioner``            | Define your preconditioner. [Parameters]         |
+-------------------------------+--------------------------------------------------+
| ``VerbosityLevel``            | Define the verbosity of the solver. [int: 0,1]   |
+-------------------------------+--------------------------------------------------+
| ``Tolerance``                 | Define the tolerance, >0. [double]               |
+-------------------------------+--------------------------------------------------+
| ``AbsoluteTolerance``         | Define the absolute tolerance, >= 0. [double]    |
+-------------------------------+--------------------------------------------------+
| ``MaximumNumberOfIterations`` | Maximum number of Krylov iterations, >= 1. [int] |
+-------------------------------+--------------------------------------------------+


**Example with a Jacobi preconditioner:**

.. code-block:: cpp

    problem.setLinearSolver(ctx, "HyprePCG",
                            {{"VerbosityLevel", 1},
                             {"Tolerance", 1e-12},
                             {"Preconditioner", diagscale},
                             {"MaximumNumberOfIterations", 5000}}) |
        or_die;

.. note::
  
  Only ``hypre`` preconditioners are allowed with the HyprePCG linear solver.

Generalized Minimal Residual (HypreGMRES)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``HypreGMRES``
- Dependencies:

  - ``hypre``
  - parallel

- Parameters:

+-------------------------------+--------------------------------------------------+
| Key                           | Description                                      |
+===============================+==================================================+
| ``KDim``                      | Define the Krylov dimension, >= 1. [int]         |
+-------------------------------+--------------------------------------------------+
| ``Preconditioner``            | Define your preconditioner. [Parameters]         |
+-------------------------------+--------------------------------------------------+
| ``VerbosityLevel``            | Define the verbosity of the solver. [int: 0,1]   |
+-------------------------------+--------------------------------------------------+
| ``Tolerance``                 | Define the tolerance, >0. [double]               |
+-------------------------------+--------------------------------------------------+
| ``AbsoluteTolerance``         | Define the absolute tolerance, >= 0. [double]    |
+-------------------------------+--------------------------------------------------+
| ``MaximumNumberOfIterations`` | Maximum number of Krylov iterations, >= 1. [int] |
+-------------------------------+--------------------------------------------------+

**Example with a  Parasail preconditioner:**

.. code-block:: cpp

   problem.setLinearSolver(ctx, "HypreGMRES",
                           {{"VerbosityLevel", 1},
                            {"Tolerance", 1e-12},
                            {"Preconditioner", parasail},
                            {"MaximumNumberOfIterations", 5000}}) |
       or_die;

.. note::
  
  Only ``hypre`` preconditioners are allowed with the HypreGMRES linear solver.


Flexible GMRES (HypreFGMRES)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``HypreFGMRES``
- Dependencies:

  - ``hypre``
  - parallel

- Parameters:

+-------------------------------+--------------------------------------------------+
| Key                           | Description                                      |
+===============================+==================================================+
| ``KDim``                      | Define the Krylov dimension, >= 1. [int]         |
+-------------------------------+--------------------------------------------------+
| ``Preconditioner``            | Define your preconditioner. [Parameters]         |
+-------------------------------+--------------------------------------------------+
| ``VerbosityLevel``            | Define the verbosity of the solver. [int: 0,1]   |
+-------------------------------+--------------------------------------------------+
| ``Tolerance``                 | Define the tolerance, >0. [double]               |
+-------------------------------+--------------------------------------------------+
| ``AbsoluteTolerance``         | Define the absolute tolerance, >= 0. [double]    |
+-------------------------------+--------------------------------------------------+
| ``MaximumNumberOfIterations`` | Maximum number of Krylov iterations, >= 1. [int] |
+-------------------------------+--------------------------------------------------+

**Example with a ILU preconditioner:**

.. code-block:: cpp

  problem.setLinearSolver(ctx, "HypreFGMRES",
                          {{"VerbosityLevel", 1},
                           {"Tolerance", 1e-12},
                           {"Preconditioner", ilu},
                           {"MaximumNumberOfIterations", 5000}}) |
      or_die;


.. note::
  
  Only ``hypre`` preconditioners are allowed with the HypreFGMRES linear solver.

Setup Your Preconditioner
=========================

The preconditioners are the ones of ``hypre`` (in parallel), the Jacobi
and Chebyshev smoothers of ``MFEM``, which only use products by the operator
and its diagonal, ``AssembledOperator``, which sets up a preconditioner on
the assembled version of an operator applied without assembly, and
``AmgX``, if ``mfem-mgis`` is built with it.

Algebraic MultiGrid (AMG) 
-------------------------

- Preconditioner name: ``HypreBoomerAMG``, also available as a linear
  solver under the same name (cycles of BoomerAMG, the options being then
  given directly to the solver).
- Options:

+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| Key                            | Description                                                                                                   |
+================================+===============================================================================================================+
| ``Strategy``                   | Propose strategies that can lead to improve convergence and scalability, [Elasticity, System, None]. [string] |
|                                | ``System`` by default for vector unknowns, ``None`` otherwise.                                                |
|                                | ``Elasticity`` requires the ``FiniteElementSpaceOrdering`` parameter set to ``byVDIM`` and is replaced by     |
|                                | ``System``, with a warning, otherwise or when ``hypre`` computes on a GPU.                                    |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| ``VerbosityLevel``             | Define the verbosity of the solver. [int: 0,1]                                                                |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| ``StrengthThreshold``          | Strength threshold of ``hypre``. [double]                                                                     |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| ``CoarseningType``             | Coarsening type of ``hypre``. [int]                                                                           |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| ``AggressiveCoarseningLevels`` | Number of levels of aggressive coarsening of ``hypre``. [int]                                                 |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| ``InterpolationType``          | Interpolation type of ``hypre``. [int]                                                                        |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| ``RelaxationType``             | Relaxation type of ``hypre``. [int]                                                                           |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| ``RelaxationSweeps``           | Number of sweeps before and after the coarse levels. [int, strictly positive]                                 |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| ``CycleType``                  | Type of cycle (1: V-cycle, 2: W-cycle). [int, strictly positive]                                              |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| ``NumberOfCycles``             | Number of cycles per application (1 by default). [int, strictly positive]                                     |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+
| ``MaximumNumberOfLevels``      | Maximum number of levels. [int, strictly positive]                                                            |
+--------------------------------+---------------------------------------------------------------------------------------------------------------+

The values of the options of ``hypre`` are the ones of ``hypre`` (see its
documentation).

- website: https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html

**Example:**

.. code-block:: cpp

    auto options = mfem_mgis::Parameters{{"VerbosityLevel", 0}, {"Strategy", "System"}};
    auto amg = mfem_mgis::Parameters{{"Name","HypreBoomerAMG"}, {"Options",options}};

**Reference:**

.. code-block:: text

  J. W. Ruge and K. Stüben. Algebraic multigrid (AMG). In S. F. McCormick, editor, Multigrid Methods, volume 3 of Frontiers in Applied Mathematics, pages 73–130. SIAM, Philadelphia, PA, 1987.

Euclid (HypreEuclid)
--------------------
- Preconditioner name: ``HypreEuclid``
- Options:

+--------------------+------------------------------------------------+
| Key                | Description                                    |
+====================+================================================+
| ``VerbosityLevel`` | Define the verbosity of the solver. [int: 0,1] |
+--------------------+------------------------------------------------+

- website: https://hypre.readthedocs.io/en/latest/solvers-euclid.html 

**Example:**

.. code-block:: cpp

    auto options = mfem_mgis::Parameters{{"VerbosityLevel", 0}};
    auto euclid = mfem_mgis::Parameters{{"Name","HypreEuclid"}, {"Options",options}};

**Reference:**

.. code-block:: text

  D. Hysom and A. Pothen. Efficient parallel computation of ILU(k) preconditioners. In Proceedings of Supercomputing ‘99. ACM, November 1999. Published on CDROM, ISBN #1-58113-091-0, ACM Order #415990, IEEE Computer Society Press Order # RS00197.

Incomplete LU factorization (HypreILU)
--------------------------------------

- Preconditioner name: ``HypreILU``
- Dependencies:

  - ``MFEM_HYPRE_VERSION`` >= 21900

- Options:

+-------------------------+------------------------------------------------+
| Key                     | Description                                    |
+=========================+================================================+
| ``HypreILUType``        | [int] Type of the factorisation of hypre:      |
|                         | ILU(k) (0, the default, 10, 20, 30, 40), ILUT  |
|                         | (1, 11, 21, 31, 41) or modified ILU(0) (50),   |
|                         | with a block Jacobi (0, 1), GMRES on the Schur |
|                         | complement (10, 11, 40, 41), Newton-Schulz-    |
|                         | Hotelling (20, 21), restricted additive        |
|                         | Schwarz (30, 31) or RAP (50) global method.    |
+-------------------------+------------------------------------------------+
| ``HypreILULevelOfFill`` | [int] Level of fill of the ILU(k)              |
|                         | factorisations (not allowed for the ILUT       |
|                         | ones). One by default (default of MFEM), zero  |
|                         | when hypre computes on a GPU, where a greater  |
|                         | level is refused (failure of hypre 3.2.0).     |
+-------------------------+------------------------------------------------+
| ``VerbosityLevel``      | Define the verbosity of the solver. [int: 0,1] |
+-------------------------+------------------------------------------------+

- website: https://hypre.readthedocs.io/en/latest/solvers-ilu.html

**Example:**

.. code-block:: cpp

    auto ilu = mfem_mgis::Parameters{{"Name","HypreILU"}, {"Options", mfem_mgis::Parameters{{"HypreILULevelOfFill", 1}}}};


Sparse Approximate Inverse (HypreParaSails)
-------------------------------------------

- Description: ParaSails is a parallel implementation of a sparse approximate inverse preconditioner, using a priori sparsity patterns and least-squares (Frobenius norm) minimization.
- Preconditioner name: ``HypreParaSails``
- Not supported when hypre computes on a GPU (ParaSails only computes on
  the host): the preconditioner is then refused.
- Options:

+--------------------+------------------------------------------------+
| Key                | Description                                    |
+====================+================================================+
| ``VerbosityLevel`` | Define the verbosity of the solver. [int: 0,1] |
+--------------------+------------------------------------------------+

- website: https://hypre.readthedocs.io/en/latest/solvers-parasails.html

**Example:**

.. code-block:: cpp

    auto options = mfem_mgis::Parameters{{"VerbosityLevel", 0}};
    auto parasail = mfem_mgis::Parameters{{"Name","HypreParaSails"}, {"Options",options}};

**Reference:**

.. code-block:: text

  E. Chow. A priori sparsity patterns for parallel sparse approximate inverse preconditioners. SIAM J. Sci. Comput., 21:1804–1822, 2000.

.. warning::

  ParaSails is not actively supported by the hypre development team. 

Jacobi (HypreDiagScale)
-----------------------

- Preconditioner name: ``HypreDiagScale``
- Options:

+--------------------+------------------------------------------------+
| Key                | Description                                    |
+====================+================================================+
| ``VerbosityLevel`` | Define the verbosity of the solver. [int: 0,1] |
+--------------------+------------------------------------------------+

**Example:**

.. code-block:: cpp

    auto options = mfem_mgis::Parameters{{"VerbosityLevel", 0}};
    auto diagscale = mfem_mgis::Parameters{{"Name","HypreDiagScale"}, {"Options",options}};

Jacobi smoother (OperatorJacobiSmoother)
----------------------------------------

Jacobi smoother of ``MFEM``, whose diagonal is computed by the operator: it
only uses products by the operator and its diagonal, and is thus usable with
an operator applied without assembly (``TangentOperatorComputation`` set to
``MatrixFree``), in sequential and in parallel.

- Preconditioner name: ``OperatorJacobiSmoother``
- Options:

+-------------+-----------------------------------------+
| Key         | Description                             |
+=============+=========================================+
| ``Damping`` | Damping factor, 1 by default. [double]  |
+-------------+-----------------------------------------+

Chebyshev smoother (OperatorChebyshevSmoother)
----------------------------------------------

Chebyshev smoother of ``MFEM``: polynomial of the operator scaled by the
inverse of its diagonal, the largest eigenvalue of the scaled operator being
estimated by the power method each time the operator changes (ten products
by the operator by default). As the Jacobi
smoother, it only uses products by the operator and its diagonal. The
essential degrees of freedom shall be eliminated with a diagonal set to one.

- Preconditioner name: ``OperatorChebyshevSmoother``
- Options:

+-----------------------------+------------------------------------------------+
| Key                         | Description                                    |
+=============================+================================================+
| ``Order``                   | Order of the polynomial, from 1 to 5, 2 by     |
|                             | default. [int]                                 |
+-----------------------------+------------------------------------------------+
| ``NumberOfPowerIterations`` | Iterations of the power method, 10 by default. |
|                             | [int]                                          |
+-----------------------------+------------------------------------------------+

**Example:**

.. code-block:: cpp

    auto chebyshev = mfem_mgis::Parameters{
        {"Name", "OperatorChebyshevSmoother"},
        {"Options", mfem_mgis::Parameters{{"Order", 3}}}};

Preconditioner on the assembled operator (AssembledOperator)
------------------------------------------------------------

Sets up another preconditioner on the assembled version of a tangent
operator applied without assembly (``TangentOperatorComputation`` set to
``MatrixFree``, see :ref:`mfem_mgis_batch_computations_matrix_free`).

- Preconditioner name: ``AssembledOperator``
- Options:

+--------------------+--------------------------------------------------+
| Key                | Description                                      |
+====================+==================================================+
| ``Preconditioner`` | Preconditioner set up on the assembled operator, |
|                    | given by its ``Name`` and ``Options``.           |
|                    | Required.                                        |
+--------------------+--------------------------------------------------+
| ``UpdatePeriod``   | Number of calls between two setups, 1 by         |
|                    | default. [int]                                   |
+--------------------+--------------------------------------------------+

**Example:**

.. code-block:: cpp

    auto assembled = mfem_mgis::Parameters{
        {"Name", "AssembledOperator"},
        {"Options",
         mfem_mgis::Parameters{
             {"UpdatePeriod", 2},
             {"Preconditioner",
              mfem_mgis::Parameters{{"Name", "HypreBoomerAMG"}}}}}};

Algebraic multigrid of AmgX (AmgX)
----------------------------------

- Preconditioner name: ``AmgX``, also available as a linear solver under
  the same name. The linear solver ``AmgXPCG`` is the conjugate gradient
  of ``AmgX``.
- Dependencies: ``mfem-mgis`` built with the CMake option ``enable-amgx``.
- Options: see :ref:`mfem_mgis_batch_computations_amgx`.

Inexact Newton method (AdaptiveLinearSolverTolerance)
-----------------------------------------------------

The relative tolerance of the iterative linear solver of the Newton solver,
and of the ``NewtonK0`` method of the ``FixedPoint`` solver, may be adapted
to the decrease of the residual (Eisenstat and Walker, 1996), as
``mfem::NewtonSolver::SetAdaptiveLinRtol``: loose at the first iterations,
tighter near the solution. The linear solver shall be an iterative solver of
``MFEM`` or a Krylov solver of ``hypre``; its relative tolerance is restored
at the end of each resolution (the predictions keep it).

- Option name: ``AdaptiveLinearSolverTolerance``, with the value ``None``
  (default), ``EisenstatWalker1`` or ``EisenstatWalker2``, possibly with
  options.
- The first choice uses the norm of the residual of the previous linear
  solve (one more product by the tangent operator per iteration), as in the
  paper of Eisenstat and Walker and in PETSc; MFEM multiplies it by
  ``Gamma`` and uses the safeguard of the second choice.
- The tolerance does not decrease too fast (safeguard), is bounded by the
  maximum tolerance, and is not smaller than half the ratio of the
  convergence criterion of the non linear solver to the norm of the residual
  (no oversolving of the last linear system, Kelley): it may be lower than
  the tolerance given to the linear solver. The iterative solvers stop on
  their own norm (preconditioned norm for the conjugate gradient of MFEM).
- Options:

+----------------------+---------------------------------------------------+
| Key                  | Description                                       |
+======================+===================================================+
| ``InitialTolerance`` | Tolerance of the first linear solve, 0.5 by       |
|                      | default (or the maximum tolerance if lower).      |
+----------------------+---------------------------------------------------+
| ``MaximumTolerance`` | Maximum tolerance, 0.9 by default.                |
+----------------------+---------------------------------------------------+
| ``Gamma``            | Factor of the second choice, 1 by default.        |
+----------------------+---------------------------------------------------+
| ``Alpha``            | Exponent of the second choice, the golden ratio   |
|                      | by default.                                       |
+----------------------+---------------------------------------------------+
| ``Threshold``        | Threshold of the safeguard, 0.1 by default.       |
+----------------------+---------------------------------------------------+

**Example:**

.. code-block:: cpp

  problem.setSolverParameters(
      ctx, {{"RelativeTolerance", 1e-8},
            {"AdaptiveLinearSolverTolerance",
             mfem_mgis::Parameters{
                 {"EisenstatWalker1",
                  mfem_mgis::Parameters{{"InitialTolerance", 0.1}}}}}}) |
      or_die;

Stagnation of the residual (StagnationRelativeTolerance)
--------------------------------------------------------

Option of the Newton solver: a residual which has not decreased by more than
10 % during the last two iterations, at a level of the order of the rounding
errors, is accepted as converged (see :doc:`nonlinear_algorithms`).

- Option name: ``StagnationRelativeTolerance``, a real in :math:`[0, 1[`,
  relative to the reference residual; ``0`` disables the criterion.
- Default value: :math:`10^4` times the machine epsilon of the numerical
  type (about :math:`1.2\,10^{-3}` with ``float``, :math:`2.2\,10^{-12}`
  with ``double``).

**Example:**

.. code-block:: cpp

  problem.setSolverParameters(
      ctx, {{"RelativeTolerance", 1e-8},
            {"StagnationRelativeTolerance", 1e-4}}) |
      or_die;

Setup Your NonLinearProblem
===========================

.. warning::

  This section is under construction.

FiniteElementDiscretization
---------------------------

These parameters define the mesh and the finite element space. They are
given to the constructor of the ``FiniteElementDiscretization`` class, and
to the ones of the non linear evolution problems.

+----------------------------+----------------------------------------------------------+
| Key                        | Description                                              |
+============================+==========================================================+
| FiniteElementFamily        | Finite Element Family, ex: H1 [string]                   |
+----------------------------+----------------------------------------------------------+
| MeshFileName               | Path to the mesh file [string]                           |
+----------------------------+----------------------------------------------------------+
| MeshReadMode               | FromScratch (default) or Restart (each process reads its |
|                            | part of a mesh split beforehand) [string]                |
+----------------------------+----------------------------------------------------------+
| FiniteElementOrder         | Finite Element Order, >= 1 [int]                         |
+----------------------------+----------------------------------------------------------+
| FiniteElementSpaceOrdering | Ordering of the dofs: byNODES (default), byVDIM [string] |
+----------------------------+----------------------------------------------------------+
| IntegrationRuleOrder       | Order of the integration rules of the behaviour          |
|                            | integrators, >= 0, optional [int]                        |
+----------------------------+----------------------------------------------------------+
| UnknownsSize               | Number of unknowns, >=1 [int]                            |
+----------------------------+----------------------------------------------------------+
| NumberOfUniformRefinements | Number of times the mesh is uniformly refined, >= 0 [int]|
+----------------------------+----------------------------------------------------------+
| Parallel                   | Run parallel execution [bool]                            |
+----------------------------+----------------------------------------------------------+
| Materials                  | Names associated with the material identifiers           |
|                            | [Parameters]                                             |
+----------------------------+----------------------------------------------------------+
| Boundaries                 | Names associated with the boundary identifiers           |
|                            | [Parameters]                                             |
+----------------------------+----------------------------------------------------------+
| Points                     | Points declared by the user, see :doc:`mesh/curves`      |
+----------------------------+----------------------------------------------------------+
| PointsSets                 | Sets of points declared by the user, see                 |
|                            | :doc:`mesh/curves`                                       |
+----------------------------+----------------------------------------------------------+
| GeneralVerbosityLevel      | Verbosity level [int]                                    |
+----------------------------+----------------------------------------------------------+

.. note::

  In parallel, the uniform refinements are made after splitting the mesh across the MPI
  processes, except the first ones if the mesh has less than 8 elements per process.

NonLinearEvolutionProblem
-------------------------

The constructor of the ``NonLinearEvolutionProblem`` class takes the
parameters of the ``FiniteElementDiscretization`` class and the following
ones:

+-------------------------------------+-------------------------------------------------------------+
| Key                                 | Description                                                 |
+=====================================+=============================================================+
| Hypothesis                          | Modelling hypothesis, ex: Tridimensional, required [string] |
+-------------------------------------+-------------------------------------------------------------+
| UseMultiMaterialNonLinearIntegrator | Use the multi-material non linear integrator [bool]         |
+-------------------------------------+-------------------------------------------------------------+
| NonLinearSolver                     | Non linear solver, Newton (default) or FixedPoint, see      |
|                                     | :doc:`nonlinear_algorithms`                                 |
+-------------------------------------+-------------------------------------------------------------+
| BehaviourIntegration                | PointByPoint or Batch [string]                              |
+-------------------------------------+-------------------------------------------------------------+
| GradientsComputation                | ElementByElement or Kernels [string]                        |
+-------------------------------------+-------------------------------------------------------------+
| InnerForcesComputation              | ElementByElement or Kernels [string]                        |
+-------------------------------------+-------------------------------------------------------------+
| TangentOperatorComputation          | ElementByElement, Kernels or MatrixFree [string]            |
+-------------------------------------+-------------------------------------------------------------+
| TangentOperatorAssembly             | Atomic (default) or Ordered [string]                        |
+-------------------------------------+-------------------------------------------------------------+

The last five parameters and their default values are described in
:ref:`mfem_mgis_batch_computations`.

PeriodicNonLinearEvolutionProblem
---------------------------------

The ``PeriodicNonLinearEvolutionProblem`` class is built from a finite
element discretization, either two corners of the computation domain or
the type of the condition fixing the rigid body motions
(``mfem_mgis::FIX_XMIN`` by default), and optional parameters of the non
linear evolution problem other than the ones of the finite element
discretization (for instance ``NonLinearSolver``).

CPP example:

.. code-block:: cpp

  auto ctx = mfem_mgis::Context{};
  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      ctx, mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
                                 {"FiniteElementFamily", "H1"},
                                 {"FiniteElementOrder", 2},
                                 {"UnknownsSize", 3},
                                 {"NumberOfUniformRefinements", 2},
                                 {"Parallel", true}});
  mfem_mgis::PeriodicNonLinearEvolutionProblem problem(ctx, fed);


Material parameters
-------------------

Material Identifier
^^^^^^^^^^^^^^^^^^^

- Functions: ``getMaterialsIdentifiers``
- Key: ``Material`` (one material) or ``Materials`` (a list of
  materials), which can't be given together. A material is given by its
  identifier (an integer) or by a name, interpreted as a regular
  expression.

**Example:**

.. code-block:: cpp

  problem.addPostProcessing(
      ctx, "StoredEnergy",
      {{"OutputFileName", "energy.txt"},
       {"Materials", std::vector<mfem_mgis::Parameter>{"Attr1", "Attr2"}}}) |
      or_die;

Boundary Condition parameters
-----------------------------

Boundary Identifier
^^^^^^^^^^^^^^^^^^^

- Functions: ``getBoundariesIdentifiers``
- Key: ``Boundary`` (one boundary) or ``Boundaries`` (a list of
  boundaries), which can't be given together, given as the materials.

**Example:**

.. code-block:: cpp

  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "LowerBoundary"}, {"Component", 1}}) |
      or_die;

Dirichlet Condition
-------------------

This is how to apply a dirichlet boundary condition to several boundaries.

**Example:**

.. code-block:: cpp

  // boundary conditions
  for (const auto boundary : {"left", "right"}) {
    for (const auto dof : {0, 1}) {
      problem.addBoundaryCondition(
          ctx, std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
                   problem.getFiniteElementDiscretizationPointer(), boundary,
                   dof)) |
          or_die;
    }
  }

