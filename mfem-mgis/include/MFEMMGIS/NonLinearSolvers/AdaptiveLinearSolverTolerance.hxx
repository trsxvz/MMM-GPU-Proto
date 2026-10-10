/*!
 * \file   MFEMMGIS/NonLinearSolvers/AdaptiveLinearSolverTolerance.hxx
 * \brief  This file declares the `AdaptiveLinearSolverTolerance` class,
 * the relative tolerance of the iterative linear solvers of the Newton
 * methods adapted to the decrease of the residual (inexact Newton methods)
 * \date   10/10/2026
 */

#ifndef LIB_MFEMMGIS_NONLINEARSOLVERS_ADAPTIVELINEARSOLVERTOLERANCE_HXX
#define LIB_MFEMMGIS_NONLINEARSOLVERS_ADAPTIVELINEARSOLVERTOLERANCE_HXX

#include <optional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \brief relative tolerance of the iterative linear solver of a Newton
   * method adapted to the decrease of the residual (inexact Newton method,
   * Eisenstat and Walker, SIAM J. Sci. Comput. 17(1), 1996), as
   * `mfem::NewtonSolver::SetAdaptiveLinRtol` and the `-snes_ksp_ew` option
   * of PETSc. The tolerance `η_i` of the linear solve of the iteration `i`
   * is:
   *
   * - `η_0 = InitialTolerance` at the first iteration;
   * - choice 1 (`EisenstatWalker1`):
   *   `η_i = | ||F_i|| - ||F_{i-1} + J_{i-1} s_{i-1}|| | / ||F_{i-1}||`,
   *   which requires the norm of the residual of the previous linear solve
   *   (one product by the tangent operator per iteration; after a reduction
   *   of the correction, the residual of the complete linear solve is
   *   used);
   * - choice 2 (`EisenstatWalker2`): `η_i = γ (||F_i|| / ||F_{i-1}||)^α`.
   *
   * Safeguards, as Eisenstat and Walker and PETSc: when the following value
   * is above `Threshold`, the tolerance is not smaller than
   * `η_{i-1}^((1+√5)/2)` (choice 1) or `γ η_{i-1}^α` (choice 2), so that it
   * does not decrease too fast (MFEM uses the latter for both choices,
   * which only differ for values of `Gamma` and `Alpha` other than the
   * default ones); the tolerance is bounded by `MaximumTolerance`, and is
   * not smaller than `0.5 τ / ||F_i||`, `τ` being the convergence criterion
   * on the residual, which avoids oversolving the last linear system
   * (Kelley, Iterative methods for linear and nonlinear equations, SIAM,
   * 1995). The tolerance may thus be lower than the one given to the linear
   * solver.
   *
   * Option `AdaptiveLinearSolverTolerance` of the `Newton` solver and of the
   * `NewtonK0` method of the `FixedPoint` solver: `None` (default, the
   * tolerance of the linear solver is kept), `EisenstatWalker1` or
   * `EisenstatWalker2`, given by a name, by a dictionary associating a name
   * with its options, or by a dictionary with the entries `Name` and
   * `Options`. Options: `InitialTolerance` (0.5 by default, or the maximum
   * tolerance if lower), `MaximumTolerance` (0.9), `Gamma` (1, choice 2),
   * `Alpha` (the golden ratio, choice 2) and `Threshold` (0.1), as in MFEM.
   * The linear solver shall be an iterative
   * solver of MFEM or a Krylov solver of hypre; its relative tolerance is
   * restored at the end of each resolution.
   */
  struct MFEM_MGIS_EXPORT AdaptiveLinearSolverTolerance {
    //! \brief methods
    enum struct Method { NONE, EISENSTAT_WALKER_1, EISENSTAT_WALKER_2 };
    /*!
     * \return the adaptive tolerance given by the option
     * `AdaptiveLinearSolverTolerance`, or an invalid value on error
     * \param[in, out] ctx: execution context
     * \param[in] params: parameters, which shall contain the option
     */
    [[nodiscard]] static std::optional<AdaptiveLinearSolverTolerance> read(
        Context&, const Parameters&) noexcept;
    //! \return if the tolerance is adapted
    [[nodiscard]] bool isEnabled() const noexcept {
      return this->method != Method::NONE;
    }
    //! \brief method
    Method method = Method::NONE;
    //! \brief tolerance of the first linear solve
    real initial_tolerance = real{0.5};
    //! \brief maximum tolerance
    real maximum_tolerance = real{0.9};
    //! \brief factor `γ`
    real gamma = real{1};
    //! \brief exponent `α`
    real alpha = real{1.618033988749895};
    //! \brief threshold of the safeguard
    real threshold = real{0.1};
  };  // end of struct AdaptiveLinearSolverTolerance

  /*!
   * \brief sets the relative tolerance of an iterative linear solver during
   * a resolution, and restores its initial value on destruction
   */
  struct MFEM_MGIS_EXPORT AdaptiveLinearSolverToleranceUpdater {
    /*!
     * \param[in] t: adaptive tolerance
     * \param[in] s: linear solver
     * \throws CollectiveFailure if the linear solver is not an iterative
     * solver of MFEM or a Krylov solver of hypre (the type of the linear
     * solver being the same on all the processes)
     */
    AdaptiveLinearSolverToleranceUpdater(const AdaptiveLinearSolverTolerance&,
                                         LinearSolver&);
    AdaptiveLinearSolverToleranceUpdater(
        const AdaptiveLinearSolverToleranceUpdater&) = delete;
    AdaptiveLinearSolverToleranceUpdater& operator=(
        const AdaptiveLinearSolverToleranceUpdater&) = delete;
    /*!
     * \brief sets the tolerance of the linear solve of an iteration
     * \return the tolerance
     * \param[in] i: iteration number, starting at 0
     * \param[in] norm: norm of the residual `||F_i||`
     * \param[in] goal: convergence criterion on the norm of the residual
     */
    real update(const size_type, const real, const real) noexcept;
    /*!
     * \return if the norm of the residual of the linear solve,
     * `||F_i + J_i s_i||`, is required (choice 1)
     */
    [[nodiscard]] bool requiresLinearResidualNorm() const noexcept;
    /*!
     * \brief sets the norm of the residual of the last linear solve
     * \param[in] n: norm
     */
    void setLinearResidualNorm(const real) noexcept;
    /*!
     * \return the norm of the residual of the last linear solve, if set
     * since the last update of the tolerance
     */
    [[nodiscard]] std::optional<real> getLinearResidualNorm() const noexcept;
    //! \brief destructor, restoring the initial tolerance
    ~AdaptiveLinearSolverToleranceUpdater() noexcept;

   private:
    //! \brief adaptive tolerance
    const AdaptiveLinearSolverTolerance options;
    //! \brief linear solver
    LinearSolver& linear_solver;
    //! \brief initial relative tolerance of the linear solver
    const real initial_tolerance;
    //! \brief norm of the residual at the previous iteration
    real previous_norm = 0;
    //! \brief norm of the residual of the previous linear solve, if set
    std::optional<real> linear_residual_norm;
    //! \brief previous tolerance
    real previous_tolerance = 0;
  };  // end of struct AdaptiveLinearSolverToleranceUpdater

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_NONLINEARSOLVERS_ADAPTIVELINEARSOLVERTOLERANCE_HXX */
