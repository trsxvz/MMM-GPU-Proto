/*!
 * \file   MFEMMGIS/NonLinearSolvers/NewtonSolver.hxx
 * \brief  This file declares the `NewtonSolver` class
 * \author Thomas Helfer
 * \date   29/03/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEARSOLVERS_NEWTONSOLVER_HXX
#define LIB_MFEM_MGIS_NONLINEARSOLVERS_NEWTONSOLVER_HXX

#include <limits>
#include <vector>
#include <optional>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverBase.hxx"
#include "MFEMMGIS/NonLinearSolvers/AdaptiveLinearSolverTolerance.hxx"

namespace mfem_mgis {

  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

  /*!
   * \brief custom implementation of the Newton Solver
   *
   * Options: the ones of the iterative solvers (`RelativeTolerance`,
   * `AbsoluteTolerance`, `MaximumNumberOfIterations`, `VerbosityLevel`),
   * `DiscardLinearSolverFailure`, `AdaptiveLinearSolverTolerance`
   * (relative tolerance of the iterative linear solver adapted to the
   * decrease of the residual, see `AdaptiveLinearSolverTolerance`) and
   * `StagnationRelativeTolerance` (see `stagnation_relative_tolerance`).
   */
  struct MFEM_MGIS_EXPORT NewtonSolver : public NonLinearSolverBase {
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    NewtonSolver(NonLinearEvolutionProblemImplementation<true>& p);
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    NewtonSolver(NonLinearEvolutionProblemImplementation<false>& p);
    /*!
     * \brief solve the non linear problem
     * \param[in] b: right hand side, unused
     * \param[in, out] x: initial guess, then solution
     */
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override;
    [[nodiscard]] bool setSolverParameters(Context&,
                                           const Parameters&) noexcept override;
    //! \brief destructor
    ~NewtonSolver() override;

   protected:
    /*!
     * \brief compute the correction associated with the given residual
     * \param[out] c: opposite of the Newton correction
     * \param[in] r: residual
     * \param[in] u: current estimate of the unknowns
     * \param[in, out] t: adaptive tolerance of the linear solver, if any,
     * given the norm of the residual of the linear solve if required
     * \return true on success
     */
    virtual bool computeNewtonCorrection(
        mfem::Vector& c,
        const mfem::Vector& r,
        const mfem::Vector& u,
        AdaptiveLinearSolverToleranceUpdater* const t) const noexcept;
    /*!
     * \brief settings of the relative tolerance of the linear solver adapted
     * to the decrease of the residual (option
     * `AdaptiveLinearSolverTolerance`), disabled by default
     */
    AdaptiveLinearSolverTolerance adaptive_linear_solver_tolerance;
    /*!
     * \brief relative residual below which a stagnation of the residual is
     * accepted as convergence (option `StagnationRelativeTolerance`, zero
     * disables it). The residual is computed with rounding errors: it
     * cannot go below a floor, which grows with the size of the problem and
     * with the level of the stresses, and which is larger in finite strain
     * with float (the deformation gradient, close to the identity, is
     * rounded). Below this tolerance, a residual whose smallest value has not
     * been halved during the last two iterations is considered at this
     * floor. By default, 10^4 times the machine epsilon of the numerical
     * type: about 1.2e-3 with float, 2.2e-12 with double.
     */
    real stagnation_relative_tolerance =
        10000 * std::numeric_limits<real>::epsilon();
  };  // end of struct NewtonSolver

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEARSOLVERS_NEWTONSOLVER_HXX */
