.. _linear_solver_modules:

=====================================================
Adding the linear solvers of a library
=====================================================

The linear solvers and the preconditioners of each library are registered
by a module (see :file:`src/LinearSolverModules.hxx`). A module is a pair
of functions, called by the constructors of the factories when the
library is available:

- :cxx:`declare<Library>LinearSolvers(LinearSolverFactory<parallel>&)`;
- :cxx:`declare<Library>Preconditioners(PreconditionerFactory<parallel>&)`.

The existing modules are the ones of MFEM and mfem-mgis
(:file:`src/MFEMLinearSolvers.cxx`), hypre
(:file:`src/HypreLinearSolvers.cxx`), MUMPS
(:file:`src/MUMPSSolvers.cxx`), cuDSS (:file:`src/CuDSSSolvers.cxx`),
AmgX (:file:`src/AmgXSolvers.cxx`) and PETSc
(:file:`src/PETScLinearSolvers.cxx`).

Adding a library consists of:

1. writing its module in :file:`src/<Library>LinearSolvers.cxx` (or
   :file:`src/<Library>Solvers.cxx`), compiled when the library is
   available (:cxx:`MFEM_USE_<LIBRARY>` for the libraries wrapped by MFEM,
   an option of mfem-mgis for the others, as :cxx:`MFEM_MGIS_HAVE_AMGX`);
2. declaring its functions in :file:`src/LinearSolverModules.hxx` and
   calling them in :cxx:`declareLinearSolvers` and
   :cxx:`declarePreconditioners` (:file:`src/LinearSolverModules.cxx`);
3. adding its entries to :cxx:`getUnavailableLibraryHint`
   (:file:`src/LinearSolverModules.cxx`), so that the error message of an
   unknown entry explains that the library is not available.

Each entry is registered by the :cxx:`declare` function with its
description (:cxx:`LinearSolverDescription`): library, kind, placement,
accepted preconditioners, options and a short text. The factories check
the placement against the execution settings, refuse a preconditioner
given to a linear solver which accepts none, and the
:cxx:`internals::makePreconditioner` function, called by the generators of
the linear solvers which accept preconditioners, refuses the
preconditioners of another library when the solver only accepts the ones
of its library.

A linear solver which is neither an iterative solver of MFEM nor a Krylov
solver of hypre implements the :cxx:`LinearSolverConvergenceReport`
interface, read by the :cxx:`hasConverged`,
:cxx:`getNumberOfIterationsAtConvergence`, :cxx:`getRelativeTolerance` and
:cxx:`setRelativeTolerance` functions (the last two are used by the
adaptive tolerance of the linear solver of the Newton solver).

The generic test :file:`tests/LinearSolverFactoriesTest.cxx` solves the
Laplace equation with each registered linear solver and each registered
preconditioner (given to the GMRES solver of MFEM), in sequential and in
parallel, and checks the descriptions: each new entry is thus solved
once with its default options, its specific options and behaviours being
checked by a dedicated test.
The table of the linear solvers of the user guide
(:ref:`linear_solvers`) is updated with the new entries.
