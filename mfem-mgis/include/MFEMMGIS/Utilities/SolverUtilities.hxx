/*!
 * \file   MFEMMGIS/Utilities/SolverUtilities.hxx
 * \brief  This file declares functions setting the parameters of solvers and
 * checking their convergence
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#ifndef LIB_MFEM_MGIS_UTILITIES_SOLVERUTILITIES_HXX
#define LIB_MFEM_MGIS_UTILITIES_SOLVERUTILITIES_HXX

#include <optional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \return the list of parameters for an iterative solver
   * \note the following parameters are available:
   * - `AbstractNonLinearEvolutionProblem::SolverVerbosityLevel`, aka
   *   `"VerbosityLevel"`,
   * - `AbstractNonLinearEvolutionProblem::SolverRelativeTolerance`, aka
   *   `"RelativeTolerance"`,
   * - `AbstractNonLinearEvolutionProblem::SolverAbsoluteTolerance`, aka
   *   `"AbsoluteTolerance"`,
   * - `AbstractNonLinearEvolutionProblem::SolverMaximumNumberOfIterations`, aka
   *   `"MaximumNumberOfIterations"`,
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<std::string>
  getIterativeSolverParametersList();

  /*!
   * \brief set the parameters of an iterative solver
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] s: iterative solver
   * \param[in] params: parameters
   * \return true on success
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool setSolverParameters(
      Context& ctx, IterativeSolver& s, const Parameters& params) noexcept;

#ifdef MFEM_USE_PETSC
  /*!
   * \brief set the parameters of a PETSc solver
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] s: solver
   * \param[in] params: parameters
   * \return true on success
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool setSolverParameters(
      Context& ctx,
      mfem::PetscNonlinearSolver& s,
      const Parameters& params) noexcept;
#endif /* MFEM_USE_PETSC */

  /*!
   * \brief set the parameters of an iterative solver
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in, out] s: iterative solver
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[deprecated]] void setSolverParameters(
      attributes::Throwing throwing,
      IterativeSolver& s,
      const Parameters& params);

#ifdef MFEM_USE_PETSC
  /*!
   * \brief set the parameters of a PETSc solver
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in, out] s: solver
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[deprecated]] void setSolverParameters(
      attributes::Throwing throwing,
      mfem::PetscNonlinearSolver& s,
      const Parameters& params);
#endif /* MFEM_USE_PETSC */

  /*!
   * \brief clear the errors raised by hypre during the last linear solve
   *
   * The Krylov solvers of hypre raise `HYPRE_ERROR_CONV` when they reach their
   * maximum number of iterations and MFEM does not clear it: the next hypre
   * solver driven by MFEM would then abort. This function does nothing
   * without MPI.
   */
  MFEM_MGIS_EXPORT void clearHypreErrors() noexcept;
  /*!
   * \brief check if the linear solver has converged
   * \param[in] ls: linear solver
   * \return if the linear solver has converged: convergence reported by the
   * iterative solvers of MFEM, the Krylov solvers of hypre (number of
   * iterations below the maximum) and the solvers of AmgX; for cuDSS,
   * relative residuals of the last solutions below 1e-3
   * (`CuDSSDirectSolver::hasSucceeded`); true for the other solvers
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool hasConverged(
      const LinearSolver& ls) noexcept;  // end of hasConverged
  /*!
   * \brief get the number of iterations of an iterative linear solver
   *
   * Supported solvers: iterative solvers of MFEM, Krylov solvers of hypre
   * (`HyprePCG`, `HypreGMRES`, `HypreFGMRES`) and iterative solvers of AmgX
   * (`AmgXPCG`).
   *
   * \param[in] ls: linear solver
   * \return the number of iterations, empty for the other solvers
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<int>
  getNumberOfIterationsAtConvergence(const LinearSolver& ls) noexcept;
  /*!
   * \return the number of analyses (ordering and symbolic factorization)
   * made by a direct solver which reuses its analysis for the matrices of
   * the same pattern (parameter `ReuseAnalysis` of the `MUMPSSolver` and
   * `CuDSSSolver` solvers of the factory of linear solvers), empty for the
   * other solvers
   * \param[in] ls: linear solver
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<size_type> getNumberOfAnalyses(
      const LinearSolver& ls) noexcept;
  /*!
   * \return the relative tolerance of the iterative solvers of MFEM, of the
   * Krylov solvers of hypre (`HyprePCG`, `HypreGMRES`, `HypreFGMRES`) and of
   * the solvers of the other libraries which give access to it (see
   * `LinearSolverConvergenceReport`), empty for the other solvers (for
   * instance the solvers of AmgX)
   * \param[in] ls: linear solver
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<real> getRelativeTolerance(
      const LinearSolver& ls) noexcept;
  /*!
   * \brief sets the relative tolerance of an iterative linear solver
   * (iterative solvers of MFEM, Krylov solvers of hypre, solvers of the other
   * libraries which give access to it, see `LinearSolverConvergenceReport`)
   * \return false if the relative tolerance of the solver can't be set (for
   * instance for a solver of AmgX)
   * \param[in, out] ls: linear solver
   * \param[in] tol: relative tolerance
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool setRelativeTolerance(
      LinearSolver& ls, const real tol) noexcept;
  /*!
   * \return a relative tolerance reachable with the numerical type of
   * mfem-mgis: the given tolerance with double, bounded below by 100 times
   * the machine epsilon with float (1.2e-5): the relative residuals of the
   * nonlinear and linear problems stagnate at about ten epsilons in float
   * (1.2e-6 for the tests of mfem-mgis).
   * The relative tolerances given to the solvers (`RelativeTolerance` and
   * `Tolerance` parameters, `setRelativeTolerance`) are bounded this way;
   * a null tolerance, which fixes the number of iterations, and the
   * absolute tolerances, which have a unit, are not changed.
   * \param[in] tol: relative tolerance
   */
  MFEM_MGIS_EXPORT [[nodiscard]] real boundRelativeTolerance(
      const double tol) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_UTILITIES_SOLVERUTILITIES_HXX */
