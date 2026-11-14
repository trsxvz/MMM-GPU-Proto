/*!
 * \file   AdaptiveLinearSolverToleranceTest.cxx
 * \brief  Test of the relative tolerance of the iterative linear solvers
 * adapted to the decrease of the residual (`AdaptiveLinearSolverTolerance`):
 * relative tolerances of the iterative solvers of MFEM and of the Krylov
 * solvers of hypre, tolerances given by the two choices of Eisenstat and
 * Walker and by the safeguards, compared to a direct evaluation of their
 * formulas, restoration of the initial tolerance, refusal of a direct
 * solver and of invalid options, norm of the linear model at the step
 * retained by a line search.
 * \date   10/10/2026
 */

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/NonLinearSolvers/AdaptiveLinearSolverTolerance.hxx"
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

/*!
 * \brief relative tolerance set and read back. In float, the relative
 * tolerances are bounded below by 100 times the machine epsilon (see
 * `mfem_mgis::boundRelativeTolerance`)
 */
static constexpr auto relative_tolerance =
    mfem_mgis::unit_tests::tolerance(1e-7, 1e-3);

/*!
 * \brief initial relative tolerance of the linear solvers, restored by the
 * adaptive tolerance, through the same bound
 */
static constexpr auto initial_tolerance =
    mfem_mgis::unit_tests::tolerance(1e-9, 1e-4);

/*!
 * \brief checks the relative tolerance of an iterative solver
 * \param[in] s: linear solver
 * \param[in] name: name of the solver
 */
static void checkRelativeTolerance(mfem_mgis::LinearSolver& s,
                                   const std::string& name) {
  using namespace mfem_mgis;
  check(setRelativeTolerance(s, relative_tolerance),
        name + ": tolerance not set");
  const auto otol = getRelativeTolerance(s);
  check(otol.has_value() && areClose(*otol, relative_tolerance),
        name + ": invalid relative tolerance");
}  // end of checkRelativeTolerance

//! \brief adaptive tolerance of the given method, with the default options
static mfem_mgis::AdaptiveLinearSolverTolerance makeAdaptiveTolerance(
    const char* const method) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  const auto ot = AdaptiveLinearSolverTolerance::read(
      ctx, Parameters{{"AdaptiveLinearSolverTolerance", method}});
  if (!ot.has_value()) {
    std::cerr << ctx.getErrorMessage() << '\n';
    std::exit(EXIT_FAILURE);
  }
  return *ot;
}  // end of makeAdaptiveTolerance

//! \brief tolerances of the second choice of Eisenstat and Walker
static void checkSecondChoice() {
  using namespace mfem_mgis;
  const auto t = makeAdaptiveTolerance("EisenstatWalker2");
  auto s = mfem::CGSolver{};
  s.SetRelTol(initial_tolerance);
  {
    auto u = AdaptiveLinearSolverToleranceUpdater{t, s};
    check(!u.requiresLinearResidualNorm(),
          "second choice: linear residual required");
    // norms of the residual and convergence criterion
    const auto goal = 1e-12;
    const auto norms = std::vector<double>{1, 1e-1, 1e-3, 1e-6};
    // direct evaluation of the formulas
    auto previous = 0.;
    for (std::size_t i = 0; i != norms.size(); ++i) {
      auto eta = static_cast<double>(t.initial_tolerance);
      if (i != 0) {
        eta = t.gamma * std::pow(norms[i] / norms[i - 1], t.alpha);
        const auto sg = t.gamma * std::pow(previous, t.alpha);
        if (sg > t.threshold) {
          eta = std::max(eta, sg);
        }
      }
      eta = std::min(std::max(eta, goal / (2 * norms[i])),
                     static_cast<double>(t.maximum_tolerance));
      const auto v = u.update(i, norms[i], goal);
      check(areClose(v, eta), "second choice: tolerance of iteration " +
                                  std::to_string(i) + " (" + std::to_string(v) +
                                  " instead of " + std::to_string(eta) + ")");
      check(areClose(*getRelativeTolerance(s), eta),
            "second choice: tolerance of the linear solver at iteration " +
                std::to_string(i));
      previous = eta;
    }
    // the safeguard, active at the second and third iterations, makes the
    // tolerance decrease slowly; the last one is the formula
    check(areClose(previous, std::pow(1e-3, t.alpha)),
          "second choice: tolerance of the last iteration");
    // no oversolving: the tolerance is not smaller than half the ratio of
    // the convergence criterion to the norm of the residual
    const auto v = u.update(4, 1e-10, 1e-8);
    check(areClose(v, t.maximum_tolerance),
          "second choice: oversolving of the last linear system");
  }
  check(areClose(*getRelativeTolerance(s), initial_tolerance),
        "second choice: initial tolerance not restored");
}  // end of checkSecondChoice

//! \brief tolerances of the first choice of Eisenstat and Walker
static void checkFirstChoice() {
  using namespace mfem_mgis;
  auto ctx = Context{};
  // the formula of the first choice does not depend on gamma (MFEM
  // multiplies it by gamma)
  const auto ot = AdaptiveLinearSolverTolerance::read(
      ctx, Parameters{{"AdaptiveLinearSolverTolerance",
                       Parameters{{"EisenstatWalker1",
                                   Parameters{{"InitialTolerance", 0.05},
                                              {"Gamma", 0.5}}}}}});
  check(ot.has_value(), "first choice: options refused");
  if (!ot.has_value()) {
    return;
  }
  auto s = mfem::GMRESSolver{};
  s.SetRelTol(initial_tolerance);
  {
    auto u = AdaptiveLinearSolverToleranceUpdater{*ot, s};
    check(u.requiresLinearResidualNorm(),
          "first choice: linear residual not required");
    check(areClose(u.update(0, 1, 1e-12), 0.05),
          "first choice: initial tolerance");
    // |0.1 - 0.06| / 1, the safeguard 0.05^α being below the threshold
    u.setLinearResidualNorm(0.06);
    check(areClose(u.update(1, 0.1, 1e-12), 0.04),
          "first choice: tolerance of the second iteration");
    // without the norm of the residual of the linear solve, formula of the
    // second choice, γ (0.01 / 0.1)^α, the safeguard 0.04^((1+√5)/2) being
    // below the threshold
    check(areClose(u.update(2, 0.01, 1e-12),
                   ot->gamma * std::pow(0.1, ot->alpha)),
          "first choice: tolerance without the residual of the linear solve");
  }
  check(areClose(*getRelativeTolerance(s), initial_tolerance),
        "first choice: initial tolerance not restored");
  // three iterations with the safeguard of the first choice, η^((1+√5)/2)
  // (Eisenstat and Walker), which does not depend on gamma
  const auto ot2 = AdaptiveLinearSolverTolerance::read(
      ctx, Parameters{{"AdaptiveLinearSolverTolerance",
                       Parameters{{"EisenstatWalker1",
                                   Parameters{{"InitialTolerance", 0.5},
                                              {"Gamma", 0.5}}}}}});
  check(ot2.has_value(), "first choice: options refused");
  if (!ot2.has_value()) {
    return;
  }
  const auto golden_ratio = (1 + std::sqrt(5.)) / 2;
  auto u = AdaptiveLinearSolverToleranceUpdater{*ot2, s};
  check(areClose(u.update(0, 1, 1e-12), 0.5),
        "first choice: initial tolerance");
  // |0.1 - 0.06| / 1 = 0.04, below the safeguard 0.5^φ
  u.setLinearResidualNorm(0.06);
  const auto eta1 = std::pow(0.5, golden_ratio);
  check(areClose(u.update(1, 0.1, 1e-12), eta1),
        "first choice: safeguard of the second iteration");
  // |0.01 - 0.02| / 0.1 = 0.1, below the safeguard eta1^φ
  u.setLinearResidualNorm(0.02);
  check(areClose(u.update(2, 0.01, 1e-12), std::pow(eta1, golden_ratio)),
        "first choice: safeguard of the third iteration");
}  // end of checkFirstChoice

/*!
 * \brief norm of the linear model at the step retained by a line search:
 * given the norm of the residual of the linear solve, `lr = J c - r`, and
 * its scalar product with the residual `r`, the norm of
 * `(1 - θ) r - θ lr` is the one computed on the vectors, and gives the
 * next tolerance of the first choice of Eisenstat and Walker
 */
static void checkStepLength() {
  using namespace mfem_mgis;
  auto ctx = Context{};
  // the safeguard of an initial tolerance of 0.05 is below its threshold:
  // the next tolerance is the one of the formula
  const auto ot = AdaptiveLinearSolverTolerance::read(
      ctx, Parameters{{"AdaptiveLinearSolverTolerance",
                       Parameters{{"EisenstatWalker1",
                                   Parameters{{"InitialTolerance", 0.05}}}}}});
  check(ot.has_value(), "step length: options refused");
  if (!ot.has_value()) {
    return;
  }
  // in float, the norm is computed with the rounding errors of this type
  auto isClose = [](const double a, const double b, const double tolerance) {
    return std::abs(a - b) <= tolerance * std::max(std::abs(a), std::abs(b));
  };
  constexpr auto norm_tolerance = mfem_mgis::unit_tests::tolerance(1e-13, 1e-5);
  constexpr auto eta_tolerance = mfem_mgis::unit_tests::tolerance(1e-12, 1e-4);
  auto s = mfem::CGSolver{};
  s.SetRelTol(initial_tolerance);
  const auto r = std::vector<double>{3, -1, 2, 0.5};
  const auto lr = std::vector<double>{0.4, 0.3, -0.2, 0.1};
  auto dot = [](const std::vector<double>& a, const std::vector<double>& b) {
    auto v = 0.;
    for (std::size_t i = 0; i != a.size(); ++i) {
      v += a[i] * b[i];
    }
    return v;
  };
  const auto nr = std::sqrt(dot(r, r));
  const auto nlr = std::sqrt(dot(lr, lr));
  // norm of the residual at the next iteration
  const auto f1 = 0.6 * nr;
  auto previous = 0.;
  for (const auto theta : {1., 0.5, 0.25}) {
    auto u = AdaptiveLinearSolverToleranceUpdater{*ot, s};
    static_cast<void>(u.update(0, nr, 1e-12));
    u.setLinearResidual(nlr, dot(lr, r));
    u.setStepLength(theta);
    auto m = std::vector<double>(r.size());
    for (std::size_t i = 0; i != r.size(); ++i) {
      m[i] = (1 - theta) * r[i] - theta * lr[i];
    }
    const auto expected = std::sqrt(dot(m, m));
    const auto name = "step length " + std::to_string(theta);
    const auto lm = u.getLinearResidualNorm();
    check(lm.has_value() && isClose(*lm, expected, norm_tolerance),
          name + ": norm of the linear model (" +
              std::to_string(lm.value_or(-1)) + " instead of " +
              std::to_string(expected) + ")");
    // the norms of the three step lengths differ
    check(!isClose(expected, previous, 1e-2),
          name + ": norm of the linear model not discriminating");
    previous = expected;
    // a second declaration of the step length changes nothing
    u.setStepLength(theta);
    check(u.getLinearResidualNorm() == lm,
          name + ": step length declared twice");
    const auto eta = std::abs(f1 - expected) / nr;
    const auto v = u.update(1, f1, 1e-12);
    check(isClose(v, eta, eta_tolerance),
          name + ": next tolerance (" + std::to_string(v) + " instead of " +
              std::to_string(eta) + ")");
  }
  // without the scalar product, the norm of the residual of the complete
  // linear solve is kept
  auto u = AdaptiveLinearSolverToleranceUpdater{*ot, s};
  static_cast<void>(u.update(0, nr, 1e-12));
  u.setLinearResidualNorm(nlr);
  u.setStepLength(0.5);
  const auto n = u.getLinearResidualNorm();
  check(n.has_value() && isClose(*n, nlr, norm_tolerance),
        "step length without scalar product: norm changed");
  // nothing is done without the norm of the residual of the linear solve
  static_cast<void>(u.update(1, f1, 1e-12));
  u.setStepLength(0.5);
  check(!u.getLinearResidualNorm().has_value(),
        "step length without residual of the linear solve: norm set");
}  // end of checkStepLength

//! \brief a direct solver is refused
static void checkDirectSolver() {
  using namespace mfem_mgis;
#ifdef MFEM_MGIS_HAVE_SUITESPARSE
  auto s = mfem::UMFPackSolver{};
  check(!getRelativeTolerance(s).has_value(),
        "direct solver: relative tolerance");
  check(!setRelativeTolerance(s, 1e-7), "direct solver: tolerance set");
  const auto t = makeAdaptiveTolerance("EisenstatWalker2");
  auto refused = false;
  try {
    auto u = AdaptiveLinearSolverToleranceUpdater{t, s};
  } catch (CollectiveFailure&) {
    refused = true;
  }
  check(refused, "direct solver: accepted");
#endif /* MFEM_MGIS_HAVE_SUITESPARSE */
}  // end of checkDirectSolver

//! \brief invalid options are refused
static void checkRefusals() {
  using namespace mfem_mgis;
  // refusal of a value, with part of the expected error message
  auto check_refusal = [](const Parameter& value, const std::string& name,
                          const char* const message) {
    auto ctx = Context{};
    const auto ot = AdaptiveLinearSolverTolerance::read(
        ctx, Parameters{{"AdaptiveLinearSolverTolerance", value}});
    check(!ot.has_value(), name + " is accepted");
    // reading the error message clears it
    const auto e = ctx.getErrorMessage();
    check(e.find(message) != std::string::npos,
          name + ": unexpected error message (" + e + ")");
  };
  check_refusal("Unknown", "an unknown method", "invalid method");
  check_refusal(Parameters{{"None", Parameters{{"Gamma", 0.5}}}},
                "an option of 'None'", "no option expected");
  check_refusal(Parameters{{"RelativeTolerance", 1e-4},
                           {"MaximumNumberOfIterations", 10}},
                "a dictionary of options without name",
                "invalid value of the parameter");
  auto ew2 = [](const char* const k, const double v) {
    return Parameters{{"EisenstatWalker2", Parameters{{k, v}}}};
  };
  check_refusal(ew2("MaximumTolerance", 1), "a maximum tolerance of 1",
                "maximum tolerance");
  check_refusal(ew2("InitialTolerance", 0.95),
                "an initial tolerance above the maximum one",
                "initial tolerance");
  check_refusal(ew2("Gamma", 0), "a null gamma", "gamma");
  check_refusal(ew2("Alpha", 1), "an alpha of 1", "alpha");
  check_refusal(ew2("Alpha", 2.5), "an alpha above 2", "alpha");
  check_refusal(ew2("Threshold", 1), "a threshold of 1", "threshold");
  check_refusal(ew2("Unknown", 1), "an unknown option", "Unknown");
  auto ctx = Context{};
  // the default initial tolerance is bounded by the maximum one
  const auto omax = AdaptiveLinearSolverTolerance::read(
      ctx, Parameters{{"AdaptiveLinearSolverTolerance",
                       ew2("MaximumTolerance", 0.3)}});
  check(omax.has_value() && areClose(omax->initial_tolerance, 0.3),
        "default initial tolerance above a maximum tolerance of 0.3");
  const auto ot = AdaptiveLinearSolverTolerance::read(
      ctx, Parameters{{"AdaptiveLinearSolverTolerance", "None"}});
  check(ot.has_value() && (!ot->isEnabled()), "'None' refused or enabled");
}  // end of checkRefusals

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  {
    auto cg = mfem::CGSolver{};
    checkRelativeTolerance(cg, "CGSolver");
    auto gmres = mfem::GMRESSolver{};
    checkRelativeTolerance(gmres, "GMRESSolver");
#ifdef MFEM_USE_MPI
    auto pcg = mfem::HyprePCG{MPI_COMM_WORLD};
    checkRelativeTolerance(pcg, "HyprePCG");
    auto hgmres = mfem::HypreGMRES{MPI_COMM_WORLD};
    checkRelativeTolerance(hgmres, "HypreGMRES");
    auto hfgmres = mfem::HypreFGMRES{MPI_COMM_WORLD};
    checkRelativeTolerance(hfgmres, "HypreFGMRES");
#endif /* MFEM_USE_MPI */
  }
  checkSecondChoice();
  checkFirstChoice();
  checkStepLength();
  checkDirectSolver();
  checkRefusals();
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
