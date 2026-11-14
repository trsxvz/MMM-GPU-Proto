/*!
 * \file   MFEMMGIS/NonLinearSolvers/LineSearch.hxx
 * \brief  This file declares the `LineSearch` class, the line search of the
 * Newton methods by backtracking on the norm of the residual
 * \date   14/11/2026
 */

#ifndef LIB_MFEMMGIS_NONLINEARSOLVERS_LINESEARCH_HXX
#define LIB_MFEMMGIS_NONLINEARSOLVERS_LINESEARCH_HXX

#include <optional>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \brief line search of a Newton method by backtracking on the norm of
   * the residual (inexact Newton backtracking method, Eisenstat and Walker,
   * SIAM J. Optim. 4(2), 1994, and SIAM J. Sci. Comput. 17(1), 1996, as the
   * `bt` line search of PETSc and as NITSOL). `c` being the opposite of the
   * correction of an iteration, an approximation of `J⁻¹ F(x)`, the new
   * estimate is `x - θ c`, where the step length `θ` is the first value of
   * `1, 1/2, 1/4...` such that:
   *
   * `||F(x - θ c)|| <= [1 - t θ (1 - η)] ||F(x)||`
   *
   * where `t` is the option `SufficientDecrease` and `η` the relative
   * tolerance of the linear solve which gave `c`, if adapted (see
   * `AdaptiveLinearSolverTolerance`), zero otherwise. An estimate at which
   * the integration of the behaviours fails, or whose residual is not
   * finite, is refused. The search fails if the step length reaches
   * `2^-n`, `n` being the option `MaximumNumberOfStepReductions`, without
   * being accepted: the resolution then fails.
   *
   * The behaviours are integrated with the tangent operator at each trial:
   * an accepted step, the full one in most iterations, requires no other
   * integration nor residual than the ones of the next iteration.
   *
   * Option `LineSearch` of the `Newton` solver and of the `NewtonK0` method
   * of the `FixedPoint` solver: `None` (default, the correction is always
   * applied) or `Backtracking`, given by a name, by a dictionary
   * associating a name with its options, or by a dictionary with the
   * entries `Name` and `Options`. Options: `SufficientDecrease` (1e-4, in
   * ]0, 1[) and `MaximumNumberOfStepReductions` (10, at least 1).
   */
  struct MFEM_MGIS_EXPORT LineSearch {
    //! \brief methods
    enum struct Method { NONE, BACKTRACKING };
    //! \brief step retained by a search
    struct Result {
      //! \brief if true, an acceptable step was found
      bool success = false;
      //! \brief step length `θ` of the last trial
      real step_length = 1;
      //! \brief norm of the residual at the accepted estimate, on success
      real norm = 0;
      //! \brief number of reductions of the step
      size_type reductions = 0;
    };
    /*!
     * \brief function integrating the behaviours at an estimate of the
     * unknowns and computing the residual: returns the norm of the
     * residual, or nothing if the integration failed
     */
    using Evaluation = std::function<std::optional<real>(const mfem::Vector&)>;
    /*!
     * \return the line search given by the option `LineSearch`, or an
     * invalid value on error
     * \param[in, out] ctx: execution context
     * \param[in] params: parameters, which shall contain the option
     */
    [[nodiscard]] static std::optional<LineSearch> read(
        Context& ctx, const Parameters& params) noexcept;
    //! \return if a line search is made
    [[nodiscard]] bool isEnabled() const noexcept {
      return this->method != Method::NONE;
    }
    /*!
     * \brief searches a step length satisfying the decrease criterion. The
     * trials are computed from the current estimate, which is not modified.
     * On success, the behaviours are integrated and the residual is
     * computed at the new estimate; on failure, at the last trial.
     * \return the step retained
     * \param[out] xn: new estimate `x - θ c` (the last trial on failure)
     * \param[in] x: current estimate
     * \param[in] c: opposite of the correction
     * \param[in] norm: norm of the residual at the current estimate
     * \param[in] eta: relative tolerance of the linear solve which gave the
     * correction, if adapted, zero otherwise
     * \param[in] evaluate: integration of the behaviours and computation of
     * the residual at a trial
     * \pre the line search is enabled
     */
    [[nodiscard]] Result search(mfem::Vector& xn,
                                const mfem::Vector& x,
                                const mfem::Vector& c,
                                const real norm,
                                const real eta,
                                const Evaluation& evaluate) const;
    //! \brief method
    Method method = Method::NONE;
    //! \brief coefficient `t` of the decrease criterion
    real sufficient_decrease = real{1e-4};
    //! \brief maximum number of reductions of the step of an iteration
    size_type maximum_number_of_step_reductions = 10;
  };  // end of struct LineSearch

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_NONLINEARSOLVERS_LINESEARCH_HXX */
