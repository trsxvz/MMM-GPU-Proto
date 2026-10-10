/*!
 * \file   src/AdaptiveLinearSolverTolerance.cxx
 * \brief  This file implements the `AdaptiveLinearSolverTolerance` class
 * \date   10/10/2026
 */

#include <cmath>
#include <tuple>
#include <string>
#include <vector>
#include <algorithm>
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/NonLinearSolvers/AdaptiveLinearSolverTolerance.hxx"

namespace mfem_mgis {

  std::optional<AdaptiveLinearSolverTolerance>
  AdaptiveLinearSolverTolerance::read(Context& ctx,
                                      const Parameters& params) noexcept {
    const auto* const option = "AdaptiveLinearSolverTolerance";
    const auto oa = getFactoryArgument(ctx, params, option);
    if (isInvalid(oa)) {
      return {};
    }
    const auto& [n, opts] = *oa;
    auto t = AdaptiveLinearSolverTolerance{};
    if (n == "None") {
      if (!opts.empty()) {
        return ctx.registerErrorMessage(
            "AdaptiveLinearSolverTolerance: no option expected with 'None'");
      }
      return t;
    }
    if ((n != "EisenstatWalker1") && (n != "EisenstatWalker2")) {
      return ctx.registerErrorMessage(
          "AdaptiveLinearSolverTolerance: invalid method '" + n +
          "' (expected 'None', 'EisenstatWalker1' or 'EisenstatWalker2')");
    }
    t.method = (n == "EisenstatWalker1") ? Method::EISENSTAT_WALKER_1
                                         : Method::EISENSTAT_WALKER_2;
    if (!checkParameters(
            ctx, opts,
            std::vector<std::string>{"InitialTolerance", "MaximumTolerance",
                                     "Gamma", "Alpha", "Threshold"})) {
      return {};
    }
    const auto o0 =
        get_if<double>(ctx, opts, "InitialTolerance", t.initial_tolerance);
    const auto omax =
        get_if<double>(ctx, opts, "MaximumTolerance", t.maximum_tolerance);
    const auto ogamma = get_if<double>(ctx, opts, "Gamma", t.gamma);
    const auto oalpha = get_if<double>(ctx, opts, "Alpha", t.alpha);
    const auto othreshold = get_if<double>(ctx, opts, "Threshold", t.threshold);
    if (isInvalid(o0) || isInvalid(omax) || isInvalid(ogamma) ||
        isInvalid(oalpha) || isInvalid(othreshold)) {
      return {};
    }
    // ranges given by Eisenstat and Walker
    if ((*omax <= 0) || (*omax >= 1)) {
      return ctx.registerErrorMessage(
          "AdaptiveLinearSolverTolerance: the maximum tolerance must be in "
          "]0, 1[");
    }
    if ((*o0 <= 0) || (*o0 > *omax)) {
      return ctx.registerErrorMessage(
          "AdaptiveLinearSolverTolerance: the initial tolerance must be in "
          "]0, MaximumTolerance]");
    }
    if ((*ogamma <= 0) || (*ogamma > 1)) {
      return ctx.registerErrorMessage(
          "AdaptiveLinearSolverTolerance: gamma must be in ]0, 1]");
    }
    if ((*oalpha <= 1) || (*oalpha > 2)) {
      return ctx.registerErrorMessage(
          "AdaptiveLinearSolverTolerance: alpha must be in ]1, 2]");
    }
    if ((*othreshold < 0) || (*othreshold >= 1)) {
      return ctx.registerErrorMessage(
          "AdaptiveLinearSolverTolerance: the threshold must be in [0, 1[");
    }
    t.initial_tolerance = *o0;
    t.maximum_tolerance = *omax;
    t.gamma = *ogamma;
    t.alpha = *oalpha;
    t.threshold = *othreshold;
    return t;
  }  // end of read

  /*!
   * \return the initial relative tolerance of a linear solver
   * \param[in] s: linear solver
   * \throws CollectiveFailure if the linear solver is not iterative
   */
  [[nodiscard]] static real getInitialRelativeTolerance(const LinearSolver& s) {
    const auto otol = getRelativeTolerance(s);
    if (!otol.has_value()) {
      throw CollectiveFailure(
          "AdaptiveLinearSolverTolerance: the linear solver is not an "
          "iterative solver of MFEM or a Krylov solver of hypre");
    }
    return *otol;
  }  // end of getInitialRelativeTolerance

  AdaptiveLinearSolverToleranceUpdater::AdaptiveLinearSolverToleranceUpdater(
      const AdaptiveLinearSolverTolerance& t, LinearSolver& s)
      : options(t),
        linear_solver(s),
        initial_tolerance(getInitialRelativeTolerance(s)) {
  }  // end of AdaptiveLinearSolverToleranceUpdater

  real AdaptiveLinearSolverToleranceUpdater::update(const size_type i,
                                                    const real norm,
                                                    const real goal) noexcept {
    const auto& o = this->options;
    auto eta = o.initial_tolerance;
    if ((i != 0) && (this->previous_norm > 0)) {
      if (o.method ==
          AdaptiveLinearSolverTolerance::Method::EISENSTAT_WALKER_1) {
        eta = std::abs(norm - this->linear_residual_norm.value_or(real{0})) /
              this->previous_norm;
      } else {
        eta = o.gamma * std::pow(norm / this->previous_norm, o.alpha);
      }
      // the tolerance does not decrease too fast
      const auto safeguard =
          o.gamma * std::pow(this->previous_tolerance, o.alpha);
      if (safeguard > o.threshold) {
        eta = std::max(eta, safeguard);
      }
    }
    // no oversolving of the last linear system
    if (norm > 0) {
      eta = std::max(eta, goal / (2 * norm));
    }
    eta = std::min(eta, o.maximum_tolerance);
    // the linear solver is iterative (see the constructor)
    std::ignore = setRelativeTolerance(this->linear_solver, eta);
    this->previous_norm = norm;
    this->previous_tolerance = eta;
    this->linear_residual_norm.reset();
    return eta;
  }  // end of update

  bool AdaptiveLinearSolverToleranceUpdater::requiresLinearResidualNorm()
      const noexcept {
    return this->options.method ==
           AdaptiveLinearSolverTolerance::Method::EISENSTAT_WALKER_1;
  }  // end of requiresLinearResidualNorm

  void AdaptiveLinearSolverToleranceUpdater::setLinearResidualNorm(
      const real n) noexcept {
    this->linear_residual_norm = n;
  }  // end of setLinearResidualNorm

  std::optional<real>
  AdaptiveLinearSolverToleranceUpdater::getLinearResidualNorm() const noexcept {
    return this->linear_residual_norm;
  }  // end of getLinearResidualNorm

  AdaptiveLinearSolverToleranceUpdater::
      ~AdaptiveLinearSolverToleranceUpdater() noexcept {
    std::ignore =
        setRelativeTolerance(this->linear_solver, this->initial_tolerance);
  }  // end of ~AdaptiveLinearSolverToleranceUpdater

}  // end of namespace mfem_mgis
