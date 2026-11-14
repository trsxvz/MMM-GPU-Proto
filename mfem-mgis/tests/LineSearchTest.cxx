/*!
 * \file   LineSearchTest.cxx
 * \brief  Test of the line search of the Newton methods by backtracking on
 * the norm of the residual (`LineSearch`), without mesh: options read and
 * refused, full step accepted with one evaluation, step lengths compared to
 * a direct evaluation of the decrease criterion on a scalar equation for
 * which Newton's method diverges, trials refused on an integration failure
 * or on a norm which is not finite, failure of the search, tolerance of the
 * linear solver, sufficient decrease and maximum number of reductions taken
 * into account.
 * \date   14/11/2026
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NonLinearSolvers/LineSearch.hxx"
#include "Tolerance.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

//! \return if two values are equal up to the rounding errors
static bool areClose(const double a, const double b) {
  return std::abs(a - b) <= mfem_mgis::unit_tests::tolerance(1e-14) *
                                std::max(std::abs(a), std::abs(b));
}  // end of areClose

//! \return the line search given by a value of the option `LineSearch`
static mfem_mgis::LineSearch makeLineSearch(const mfem_mgis::Parameter& value) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  const auto ol = LineSearch::read(ctx, Parameters{{"LineSearch", value}});
  if (!ol.has_value()) {
    std::cerr << ctx.getErrorMessage() << '\n';
    std::exit(EXIT_FAILURE);
  }
  return *ol;
}  // end of makeLineSearch

//! \brief options read, default values, invalid options refused
static void checkOptions() {
  using namespace mfem_mgis;
  const auto none = makeLineSearch("None");
  check(!none.isEnabled(), "'None' is enabled");
  check(!LineSearch{}.isEnabled(), "the default line search is enabled");
  const auto bt = makeLineSearch("Backtracking");
  check(bt.isEnabled() && (bt.method == LineSearch::Method::BACKTRACKING),
        "'Backtracking' is not enabled");
  check(areClose(bt.sufficient_decrease, 1e-4),
        "default sufficient decrease");
  check(bt.maximum_number_of_step_reductions == 10,
        "default maximum number of step reductions");
  const auto with_options = makeLineSearch(
      Parameters{{"Backtracking",
                  Parameters{{"SufficientDecrease", 0.25},
                             {"MaximumNumberOfStepReductions", 4}}}});
  check(areClose(with_options.sufficient_decrease, 0.25) &&
            (with_options.maximum_number_of_step_reductions == 4),
        "options given by a dictionary associating a name with its options");
  const auto with_name = makeLineSearch(
      Parameters{{"Name", "Backtracking"},
                 {"Options", Parameters{{"SufficientDecrease", 0.5}}}});
  check(areClose(with_name.sufficient_decrease, 0.5) &&
            (with_name.maximum_number_of_step_reductions == 10),
        "options given by the entries 'Name' and 'Options'");
  // refusal of a value, with part of the expected error message
  auto check_refusal = [](const Parameter& value, const std::string& name,
                          const char* const message) {
    auto ctx = Context{};
    const auto ol = LineSearch::read(ctx, Parameters{{"LineSearch", value}});
    check(!ol.has_value(), name + " is accepted");
    // reading the error message clears it
    const auto e = ctx.getErrorMessage();
    check(e.find(message) != std::string::npos,
          name + ": unexpected error message (" + e + ")");
  };
  auto bt_option = [](const char* const k, const Parameter& v) {
    return Parameters{{"Backtracking", Parameters{{k, v}}}};
  };
  check_refusal("Unknown", "an unknown method", "invalid method");
  check_refusal(
      Parameters{{"None", Parameters{{"SufficientDecrease", 1e-4}}}},
      "an option of 'None'", "no option expected");
  check_refusal(bt_option("Unknown", 1), "an unknown option", "Unknown");
  check_refusal(bt_option("SufficientDecrease", 0.),
                "a null sufficient decrease", "sufficient decrease");
  check_refusal(bt_option("SufficientDecrease", 1.),
                "a sufficient decrease of 1", "sufficient decrease");
  check_refusal(bt_option("MaximumNumberOfStepReductions", 0),
                "a null maximum number of step reductions",
                "maximum number of step reductions");
}  // end of checkOptions

/*!
 * \brief a full step which decreases the residual is accepted with one
 * evaluation
 */
static void checkFullStep() {
  using namespace mfem_mgis;
  const auto ls = makeLineSearch("Backtracking");
  auto x = mfem::Vector(3);
  auto c = mfem::Vector(3);
  auto xn = mfem::Vector(3);
  for (int i = 0; i != 3; ++i) {
    x(i) = real(i + 1);
    c(i) = real(2 * i + 1) / 4;
  }
  auto evaluations = 0;
  const auto r = ls.search(xn, x, c, 2, 0,
                           [&evaluations](const mfem::Vector&) noexcept {
                             ++evaluations;
                             return std::optional<real>{real{1}};
                           });
  check(r.success && (r.step_length == 1) && (r.reductions == 0) &&
            (r.norm == 1),
        "full step: step retained");
  check(evaluations == 1,
        "full step: " + std::to_string(evaluations) + " evaluations");
  auto same = true;
  for (int i = 0; i != 3; ++i) {
    same = same && (xn(i) == x(i) - c(i)) && (x(i) == real(i + 1));
  }
  check(same, "full step: new estimate");
}  // end of checkFullStep

/*!
 * \brief Newton's method on `atan(x) = 0` from 1.5, where it diverges
 * without line search: the step lengths are the ones of a direct evaluation
 * of the decrease criterion, and the iterations converge
 */
static void checkScalarEquation() {
  using namespace mfem_mgis;
  const auto ls = makeLineSearch("Backtracking");
  // without line search, the iterates move away from the solution
  {
    auto x = real{1.5};
    for (int i = 0; i != 5; ++i) {
      x -= std::atan(x) * (1 + x * x);
    }
    check(std::abs(x) > 100, "scalar equation: Newton's method converges");
  }
  auto x = mfem::Vector(1);
  auto c = mfem::Vector(1);
  auto xn = mfem::Vector(1);
  x(0) = real{1.5};
  auto reduced = 0;
  auto converged = false;
  constexpr auto goal = mfem_mgis::unit_tests::tolerance(1e-12, 1e-5);
  for (int i = 0; i != 20; ++i) {
    const auto norm = std::abs(std::atan(x(0)));
    if (norm < goal) {
      converged = true;
      break;
    }
    c(0) = std::atan(x(0)) * (1 + x(0) * x(0));
    // direct evaluation of the criterion
    auto theta = real{1};
    auto reductions = size_type{};
    while (!(std::abs(std::atan(x(0) - theta * c(0))) <=
             (1 - ls.sufficient_decrease * theta) * norm)) {
      theta /= 2;
      ++reductions;
    }
    const auto r =
        ls.search(xn, x, c, norm, 0, [](const mfem::Vector& u) {
          return std::optional<real>{std::abs(std::atan(u(0)))};
        });
    check(r.success && (r.step_length == theta) &&
              (r.reductions == reductions),
          "scalar equation: step length of iteration " + std::to_string(i) +
              " (" + std::to_string(r.step_length) + " instead of " +
              std::to_string(theta) + ")");
    check(r.success && (r.norm < norm) &&
              areClose(r.norm, std::abs(std::atan(xn(0)))),
          "scalar equation: norm of iteration " + std::to_string(i));
    if (!r.success) {
      break;
    }
    if (r.reductions != 0) {
      ++reduced;
    }
    x = xn;
  }
  check(converged, "scalar equation: no convergence");
  check(reduced > 0, "scalar equation: no step reduced");
}  // end of checkScalarEquation

/*!
 * \brief trials refused: failures of the evaluation (integration failures),
 * norms which are not finite, no acceptable step
 */
static void checkRefusedTrials() {
  using namespace mfem_mgis;
  const auto ls = makeLineSearch("Backtracking");
  auto x = mfem::Vector(2);
  auto c = mfem::Vector(2);
  auto xn = mfem::Vector(2);
  x(0) = real{1};
  x(1) = real{-2};
  c(0) = real{3};
  c(1) = real{0.5};
  // the first three evaluations fail
  {
    auto evaluations = 0;
    const auto r = ls.search(
        xn, x, c, 1, 0, [&evaluations](const mfem::Vector&) noexcept {
          ++evaluations;
          return evaluations <= 3 ? std::optional<real>{}
                                  : std::optional<real>{real{0.5}};
        });
    check(r.success && (r.reductions == 3) && (r.step_length == real{0.125}) &&
              (evaluations == 4),
          "failures of the evaluation: " + std::to_string(r.reductions) +
              " reductions");
    check((xn(0) == x(0) - c(0) / 8) && (xn(1) == x(1) - c(1) / 8),
          "failures of the evaluation: new estimate");
  }
  // norms which are not finite
  {
    auto evaluations = 0;
    const auto r = ls.search(
        xn, x, c, 1, 0, [&evaluations](const mfem::Vector&) noexcept {
          ++evaluations;
          if (evaluations == 1) {
            return std::optional<real>{
                std::numeric_limits<real>::quiet_NaN()};
          }
          if (evaluations == 2) {
            return std::optional<real>{std::numeric_limits<real>::infinity()};
          }
          return std::optional<real>{real{0.5}};
        });
    check(r.success && (r.reductions == 2) && (r.step_length == real{0.25}),
          "norms which are not finite: " + std::to_string(r.reductions) +
              " reductions");
  }
  // no acceptable step: failure after the maximum number of reductions, the
  // current estimate being unchanged
  for (const auto n : {10, 3}) {
    auto l = ls;
    if (n != 10) {
      l = makeLineSearch(
          Parameters{{"Backtracking",
                      Parameters{{"MaximumNumberOfStepReductions", n}}}});
    }
    auto evaluations = 0;
    const auto r =
        l.search(xn, x, c, 1, 0, [&evaluations](const mfem::Vector&) noexcept {
          ++evaluations;
          return std::optional<real>{real{2}};
        });
    const auto theta = std::ldexp(real{1}, -n);
    check((!r.success) && (r.reductions == n) && (r.step_length == theta) &&
              (evaluations == n + 1),
          "no acceptable step: " + std::to_string(r.reductions) +
              " reductions, " + std::to_string(evaluations) +
              " evaluations (" + std::to_string(n) +
              " reductions expected)");
    check((x(0) == real{1}) && (x(1) == real{-2}),
          "no acceptable step: current estimate modified");
    check((xn(0) == x(0) - theta * c(0)) && (xn(1) == x(1) - theta * c(1)),
          "no acceptable step: last trial");
  }
}  // end of checkRefusedTrials

/*!
 * \brief the tolerance of the linear solver and the sufficient decrease
 * are taken into account by the decrease criterion
 */
static void checkCriterion() {
  using namespace mfem_mgis;
  auto x = mfem::Vector(1);
  auto c = mfem::Vector(1);
  auto xn = mfem::Vector(1);
  x(0) = real{1};
  c(0) = real{1};
  // norm of the residual at the full step, half of it at the reduced ones
  auto search = [&](const LineSearch& ls, const real full, const real eta) {
    auto evaluations = 0;
    return ls.search(xn, x, c, 1, eta,
                     [&evaluations, full](const mfem::Vector&) noexcept {
                       ++evaluations;
                       return std::optional<real>{
                           evaluations == 1 ? full : full / 2};
                     });
  };
  // sufficient decrease of 0.1: the full step is accepted below 0.9 with
  // an exact linear solve, below 0.99 with a tolerance of 0.9
  const auto ls = makeLineSearch(
      Parameters{{"Backtracking", Parameters{{"SufficientDecrease", 0.1}}}});
  check(search(ls, real{0.95}, 0).reductions == 1,
        "decrease of 5 % accepted with an exact linear solve");
  check(search(ls, real{0.95}, real{0.9}).reductions == 0,
        "decrease of 5 % refused with a tolerance of the linear solve of 0.9");
  check(search(ls, real{0.85}, 0).reductions == 0,
        "decrease of 15 % refused with a sufficient decrease of 0.1");
  // default sufficient decrease: a decrease of 5 % is accepted
  check(search(makeLineSearch("Backtracking"), real{0.95}, 0).reductions == 0,
        "decrease of 5 % refused with the default sufficient decrease");
  // an increase is always refused
  check(search(makeLineSearch("Backtracking"), real{1.5}, real{0.9})
                .reductions == 1,
        "increase of the residual accepted");
}  // end of checkCriterion

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  checkOptions();
  checkFullStep();
  checkScalarEquation();
  checkRefusedTrials();
  checkCriterion();
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
