/*!
 * \file   src/KINSOLSolver.hxx
 * \brief  This file declares the `KINSOLSolver` class, the non linear
 * solver based on KINSOL (SUNDIALS)
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_KINSOLSOLVER_HXX
#define LIB_MFEM_MGIS_KINSOLSOLVER_HXX

#include "mfem/config/config.hpp"

#ifdef MFEM_MGIS_HAVE_SUNDIALS

#include <limits>
#include <memory>
#include <optional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverBase.hxx"

namespace mfem_mgis::internals {

  // forward declaration
  struct KINSOLSolverImplementation;

  /*!
   * \brief non linear solver based on KINSOL (SUNDIALS), through the
   * `mfem::KINSolver` class of MFEM: Newton's method, with or without line
   * search, or Picard's method, optionally accelerated by Anderson's method.
   *
   * The residual is the one of the problem: KINSOL evaluates it after each
   * new estimate of the unknowns. When the integration of the behaviours
   * fails, KINSOL halves the step of Newton's method, at most five times,
   * at the first trial point of an iteration; a failure elsewhere (later in
   * the line search, or with Picard's method) stops the resolution. The
   * tangent operator is inverted by the linear solver of the problem, whose
   * failure is reported to KINSOL unless the option
   * `DiscardLinearSolverFailure` is true: KINSOL then computes the tangent
   * operator again if it is not current (Newton's method), or stops.
   *
   * KINSOL stops when the maximum norm of the residual is lower than
   * `max(RelativeTolerance * r_ref, AbsoluteTolerance)`, where `r_ref` is
   * the maximum norm of the initial residual, multiplied by the ratio of
   * the reference norm of the residual, if given (prediction), to the norm
   * of the initial residual: unlike the other non linear solvers, the
   * absolute tolerance applies to the maximum norm. A null tolerance gives
   * a fixed number of iterations. The resolution is stopped when the
   * residual stagnates (three evaluations without a decrease of more than
   * 10 %) below `StagnationRelativeTolerance` times the reference norm, and
   * accepted, as a resolution which KINSOL stops without reaching the
   * tolerance with a residual below this value.
   *
   * The behaviours are integrated again at the solution returned by KINSOL
   * if their last integration was made elsewhere (a trial point of the line
   * search), so that the state of the materials matches the solution.
   *
   * Options: the ones of the iterative solvers (`RelativeTolerance`,
   * `AbsoluteTolerance`, `MaximumNumberOfIterations`, strictly positive,
   * `VerbosityLevel`), `DiscardLinearSolverFailure`, `Strategy`
   * (`LineSearch`, by default, `Newton` or `Picard`), `MaximumSetupCalls`
   * (number of iterations between two updates of the tangent operator, one
   * by default, the tangent operator of the first iteration being kept by
   * default with Picard's method; KINSOL also updates it after a failure of
   * the linear solver or a large step with Newton's method; with a tangent
   * operator applied without assembly, only the preconditioner is kept),
   * `StagnationRelativeTolerance` (10^4 times the machine epsilon by
   * default, as for the Newton solver), for Newton's method
   * `ScaledStepTolerance` (disabled by default, the unknowns being not
   * scaled) and `MaximumNewtonStep` (maximum euclidean norm of a step, no
   * limit by default, while KINSOL limits it to 1000 times the norm of the
   * initial guess), and, for Picard's method, `Damping`,
   * `AndersonAccelerationSubspaceSize` (zero, no acceleration, by default;
   * KINSOL bounds it by the maximum number of iterations minus one),
   * `AndersonAccelerationDelay` and `AndersonAccelerationDamping`.
   *
   * KINSOL prints its error messages on the error stream.
   */
  struct KINSOLSolver final : public NonLinearSolverBase {
    //! \brief strategy of KINSOL
    enum struct Strategy { NEWTON, LINE_SEARCH, PICARD };
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    KINSOLSolver(NonLinearEvolutionProblemImplementation<true>& p);
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    KINSOLSolver(NonLinearEvolutionProblemImplementation<false>& p);
    /*!
     * \brief solve the non linear problem
     * \param[in] b: right hand side, unused
     * \param[in, out] x: initial guess, then solution
     */
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override;
    [[nodiscard]] bool setSolverParameters(Context&,
                                           const Parameters&) noexcept override;
    //! \brief destructor
    ~KINSOLSolver() override;

   private:
    friend struct KINSOLSolverImplementation;
    //! \brief strategy
    Strategy strategy = Strategy::LINE_SEARCH;
    /*!
     * \brief number of iterations between two updates of the tangent, if
     * given (one by default, except for Picard's method)
     */
    std::optional<int> maximum_setup_calls;
    //! \brief scaled step tolerance, if given (disabled by default)
    std::optional<real> scaled_step_tolerance;
    //! \brief maximum length of a Newton step, if given (no limit by default)
    std::optional<real> maximum_newton_step;
    //! \brief damping of Picard's method
    real damping = 1;
    //! \brief size of the subspace of Anderson's acceleration
    int anderson_subspace_size = 0;
    //! \brief delay of Anderson's acceleration
    int anderson_delay = 0;
    //! \brief damping of Anderson's acceleration
    real anderson_damping = 1;
    /*!
     * \brief relative residual below which a resolution stopped by KINSOL
     * without reaching the tolerance is accepted (see the option
     * `StagnationRelativeTolerance` of `NewtonSolver`)
     */
    real stagnation_relative_tolerance =
        10000 * std::numeric_limits<real>::epsilon();
  };  // end of struct KINSOLSolver

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_MGIS_HAVE_SUNDIALS */

#endif /* LIB_MFEM_MGIS_KINSOLSOLVER_HXX */
