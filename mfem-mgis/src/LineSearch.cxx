/*!
 * \file   src/LineSearch.cxx
 * \brief  This file implements the `LineSearch` class
 * \date   14/11/2026
 */

#include <cmath>
#include <string>
#include <vector>
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NonLinearSolvers/LineSearch.hxx"

namespace mfem_mgis {

  std::optional<LineSearch> LineSearch::read(
      Context& ctx, const Parameters& params) noexcept {
    const auto* const option = "LineSearch";
    const auto oa = getFactoryArgument(ctx, params, option);
    if (isInvalid(oa)) {
      return {};
    }
    const auto& [n, opts] = *oa;
    auto l = LineSearch{};
    if (n == "None") {
      if (!opts.empty()) {
        return ctx.registerErrorMessage(
            "LineSearch: no option expected with 'None'");
      }
      return l;
    }
    if (n != "Backtracking") {
      return ctx.registerErrorMessage("LineSearch: invalid method '" + n +
                                      "' (expected 'None' or 'Backtracking')");
    }
    l.method = Method::BACKTRACKING;
    if (!checkParameters(
            ctx, opts,
            std::vector<std::string>{"SufficientDecrease",
                                     "MaximumNumberOfStepReductions"})) {
      return {};
    }
    const auto odecrease =
        get_if<double>(ctx, opts, "SufficientDecrease", l.sufficient_decrease);
    const auto oreductions =
        get_if<size_type>(ctx, opts, "MaximumNumberOfStepReductions",
                          l.maximum_number_of_step_reductions);
    if (isInvalid(odecrease) || isInvalid(oreductions)) {
      return {};
    }
    if ((*odecrease <= 0) || (*odecrease >= 1)) {
      return ctx.registerErrorMessage(
          "LineSearch: the sufficient decrease must be in ]0, 1[");
    }
    if (*oreductions < 1) {
      return ctx.registerErrorMessage(
          "LineSearch: the maximum number of step reductions must be at "
          "least 1");
    }
    l.sufficient_decrease = static_cast<real>(*odecrease);
    l.maximum_number_of_step_reductions = *oreductions;
    return l;
  }  // end of read

  LineSearch::Result LineSearch::search(mfem::Vector& xn,
                                        const mfem::Vector& x,
                                        const mfem::Vector& c,
                                        const real norm,
                                        const real eta,
                                        const Evaluation& evaluate) const {
    auto r = Result{};
    subtract(x, c, xn);
    while (true) {
      const auto v = evaluate(xn);
      if (v.has_value() && std::isfinite(*v) &&
          (*v <= (1 - this->sufficient_decrease * r.step_length * (1 - eta)) *
                     norm)) {
        r.success = true;
        r.norm = *v;
        return r;
      }
      if (r.reductions >= this->maximum_number_of_step_reductions) {
        return r;
      }
      // the trial is computed from the current estimate: the rounding errors
      // of the previous trials are not accumulated
      r.step_length /= 2;
      ++(r.reductions);
      add(x, -r.step_length, c, xn);
    }
  }  // end of search

}  // end of namespace mfem_mgis
