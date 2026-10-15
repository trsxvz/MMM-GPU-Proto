/*!
 * \file   KINSOLSolverTest.cxx
 * \brief  Test of the non linear solver based on KINSOL (`KINSOL` non
 * linear solver): uniaxial tension of a cube made of a plastic material
 * with a linear hardening, in two time steps:
 *
 * - each strategy of KINSOL (Newton's method with and without line search,
 *   Picard's method, the tangent operator of the first iteration being
 *   kept, damped, or updated every five iterations) gives the solution of
 *   the Newton solver of mfem-mgis;
 * - Newton's method of KINSOL makes as many iterations as the Newton
 *   solver, within one iteration per time step (the criteria differ:
 *   maximum norm for KINSOL, euclidean norm for the Newton solver);
 * - the tangent operator updated every three iterations (modified Newton
 *   method) and Picard's method make more iterations, and Anderson's
 *   acceleration reduces the iterations of Picard's method;
 * - with a relative tolerance below the floor of the rounding errors, the
 *   resolution is stopped by the stagnation of the residual and accepted;
 * - a failure of the linear solver stops the resolution at the first
 *   iteration;
 * - the invalid options are refused.
 *
 * With the `--parallel` option, the computations are parallel (MUMPS).
 * \date   15/10/2026
 */

#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "SequentialLinearSolver.hxx"
#include "Tolerance.hxx"

//! \brief number of failed checks
static int failures = 0;
//! \brief if true, the computations are parallel
static bool parallel = false;
/*!
 * \brief tolerance on the relative difference between the solutions of
 * KINSOL and the one of the Newton solver
 */
static constexpr auto solution_tolerance =
    mfem_mgis::unit_tests::tolerance(1e-8, 3e-4);

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

//! \brief results of a computation
struct Results {
  //! \brief status
  bool status = false;
  //! \brief unknowns at the end of the computation
  std::vector<mfem_mgis::real> u;
  //! \brief total number of iterations
  mfem_mgis::size_type iterations = 0;
  //! \brief number of iterations of the time step which failed, if any
  mfem_mgis::size_type failed_iterations = 0;
};

/*!
 * \return the results of the computation
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] solver: value of the `NonLinearSolver` parameter, the Newton
 * solver being used if empty
 * \param[in] linear_solver: linear solver, MUMPS in parallel and the linear
 * solver of the sequential tests if null
 * \param[in] linear_solver_options: options of the linear solver
 * \param[in] relative_tolerance: relative tolerance of the non linear solver
 */
static Results compute(const char* const mesh,
                       const char* const library,
                       const mfem_mgis::Parameters& solver,
                       const char* const linear_solver = nullptr,
                       const mfem_mgis::Parameters& linear_solver_options = {},
                       const mfem_mgis::real relative_tolerance = 1e-12) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  auto params = Parameters{
      {"MeshFileName", mesh},
      {"FiniteElementFamily", "H1"},
      {"FiniteElementOrder", 1},
      {"UnknownsSize", 3},
      {"NumberOfUniformRefinements", 1},
      {"Materials", Parameters{{"Cube", 1}}},
      {"Boundaries",
       Parameters{{"Ymin", 1}, {"Zmin", 2}, {"Xmax", 3}, {"Xmin", 5}}},
      {"Hypothesis", "Tridimensional"},
      {"Parallel", parallel}};
  if (!solver.empty()) {
    params.insert(throwing, {{"NonLinearSolver", solver}});
  }
  auto problem = construct<NonLinearEvolutionProblem>(ctx, params) | or_die;
  problem.addBehaviourIntegrator(ctx, "Mechanics", "Cube", library,
                                 "IsotropicLinearHardeningPlasticity") |
      or_die;
  auto& m = problem.getMaterial(ctx, "Cube", 0) | or_die;
  for (auto* s : {&m.s0, &m.s1}) {
    mgis::behaviour::setExternalStateVariable(ctx, *s, "Temperature", 293.15) |
        or_die;
  }
  for (const auto& [boundary, component] :
       {std::pair{"Xmin", 0}, std::pair{"Ymin", 1}, std::pair{"Zmin", 2}}) {
    problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", boundary}, {"Component", component}}) |
        or_die;
  }
  problem.addUniformDirichletBoundaryCondition(
      ctx,
      {{"Boundary", "Xmax"},
       {"Component", 0},
       {"LoadingEvolution", [](const real t) noexcept { return 1e-2 * t; }}}) |
      or_die;
  const auto* const ls =
      (linear_solver != nullptr)
          ? linear_solver
          : (parallel ? "MUMPSSolver" : unit_tests::sequential_direct_solver);
  const auto ls_options =
      ((linear_solver == nullptr) && (!parallel))
          ? unit_tests::getSequentialDirectSolverParameters()
          : linear_solver_options;
  problem.setLinearSolver(ctx, ls, ls_options) | or_die;
  problem.setSolverParameters(
      ctx, Parameters{{"VerbosityLevel", 0},
                      {"RelativeTolerance", relative_tolerance},
                      {"AbsoluteTolerance", 0.},
                      {"MaximumNumberOfIterations", 100}}) |
      or_die;
  auto r = Results{};
  constexpr auto nsteps = 2;
  const auto dt = real{1} / nsteps;
  auto t = real{0};
  for (int i = 0; i != nsteps; ++i) {
    auto c = Context{};
    const auto o = problem.solve(c, t, dt);
    if (!o.status) {
      r.failed_iterations = o.iterations;
      return r;
    }
    r.iterations += o.iterations;
    problem.update(ctx) | or_die;
    t += dt;
  }
  const auto& u = problem.getUnknowns(ets);
  r.u.assign(u.HostRead(), u.HostRead() + u.Size());
  r.status = true;
  return r;
}  // end of compute

/*!
 * \return the maximum difference between two arrays, relative to the
 * maximum absolute value of the first one
 */
static double getRelativeDifference(const std::vector<mfem_mgis::real>& a,
                                    const std::vector<mfem_mgis::real>& b) {
  if ((a.size() != b.size()) || a.empty()) {
    return 1;
  }
  auto d = mfem_mgis::real{0};
  auto m = mfem_mgis::real{0};
  for (std::size_t i = 0; i != a.size(); ++i) {
    d = std::max(d, std::abs(a[i] - b[i]));
    m = std::max(m, std::abs(a[i]));
  }
  return m <= 0 ? d : d / m;
}  // end of getRelativeDifference

/*!
 * \brief an invalid option of KINSOL is refused
 * \param[in] mesh: mesh
 * \param[in] options: options of KINSOL
 * \param[in] message: part of the expected error message
 */
static void checkRefusal(const char* const mesh,
                         const mfem_mgis::Parameters& options,
                         const std::string& message) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  const auto o = construct<NonLinearEvolutionProblem>(
      ctx, Parameters{{"MeshFileName", mesh},
                      {"FiniteElementFamily", "H1"},
                      {"FiniteElementOrder", 1},
                      {"UnknownsSize", 3},
                      {"Hypothesis", "Tridimensional"},
                      {"NonLinearSolver", Parameters{{"KINSOL", options}}},
                      {"Parallel", parallel}});
  check(isInvalid(o), message + ": not refused");
  // reading the error message clears it
  const auto e = ctx.getErrorMessage();
  check(e.find(message) != std::string::npos,
        message + ": unexpected error message (" + e + ")");
}  // end of checkRefusal

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* mesh = "cube.mesh";
  const char* library = "src/libBehaviourTest.so";
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&mesh, "-m", "--mesh", "mesh of a cube");
  args.AddOption(&library, "-l", "--library", "behaviour library");
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "parallel computations");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  const auto reference = compute(mesh, library, Parameters{});
  check(reference.status, "Newton: resolution failed");
  std::cout << "Newton: " << reference.iterations << " iterations\n";
  // a computation with KINSOL, compared to the Newton solver
  auto run = [mesh, library, &reference](const std::string& name,
                                         const Parameters& options) {
    const auto r = compute(mesh, library, Parameters{{"KINSOL", options}});
    const auto d = getRelativeDifference(reference.u, r.u);
    std::cout << name << ": " << r.iterations
              << " iterations, relative difference with Newton " << d << '\n';
    check(r.status, name + ": resolution failed");
    check(d < solution_tolerance, name + ": solution differs from Newton's");
    return r;
  };
  const auto line_search = run("KINSOL, line search", Parameters{});
  check(line_search.iterations <= reference.iterations + 2,
        "KINSOL, line search: more iterations than the Newton solver");
  const auto newton = run("KINSOL, Newton", Parameters{{"Strategy", "Newton"}});
  check(newton.iterations <= reference.iterations + 2,
        "KINSOL, Newton: more iterations than the Newton solver");
  const auto modified = run("KINSOL, tangent updated every 3 iterations",
                            Parameters{{"MaximumSetupCalls", 3}});
  check(modified.iterations > newton.iterations,
        "KINSOL, modified Newton: not more iterations than Newton");
  const auto picard = run("KINSOL, Picard", Parameters{{"Strategy", "Picard"}});
  check(picard.iterations > newton.iterations,
        "KINSOL, Picard: not more iterations than Newton");
  const auto anderson =
      run("KINSOL, Picard, Anderson",
          Parameters{{"Strategy", "Picard"},
                     {"AndersonAccelerationSubspaceSize", 5}});
  check(anderson.iterations < picard.iterations,
        "KINSOL, Picard: Anderson's acceleration does not reduce the "
        "iterations");
  const auto damped = run("KINSOL, Picard, damped",
                          Parameters{{"Strategy", "Picard"}, {"Damping", 0.9}});
  check(damped.iterations > picard.iterations,
        "KINSOL, Picard: the damping does not slow down the iterations");
  const auto updated =
      run("KINSOL, Picard, tangent updated every 5 iterations",
          Parameters{{"Strategy", "Picard"}, {"MaximumSetupCalls", 5}});
  check(updated.iterations < picard.iterations,
        "KINSOL, Picard: the updates of the tangent operator do not reduce "
        "the iterations");
  // a relative tolerance below the floor of the rounding errors: the
  // resolution is stopped by the stagnation of the residual, and accepted
  {
    const auto r = compute(mesh, library, Parameters{{"KINSOL", Parameters{}}},
                           nullptr, {}, 1e-20);
    const auto d = getRelativeDifference(reference.u, r.u);
    std::cout << "KINSOL, stagnation: " << r.iterations
              << " iterations, relative difference with Newton " << d << '\n';
    check(r.status && (r.iterations < 50) && (d < solution_tolerance),
          "KINSOL: resolution not accepted at the floor of the rounding "
          "errors");
  }
  // a failure of the linear solver (the conjugate gradient stopped after
  // two iterations) stops the resolution at the first iteration
  {
    const auto r =
        compute(mesh, library, Parameters{{"KINSOL", Parameters{}}}, "CGSolver",
                Parameters{{"VerbosityLevel", 0},
                           {"RelativeTolerance", 1e-12},
                           {"MaximumNumberOfIterations", 2}});
    std::cout << "KINSOL, failure of the linear solver: "
              << (r.status ? "converged" : "failed") << " after "
              << r.failed_iterations << " iteration(s)\n";
    check((!r.status) && (r.failed_iterations <= 1),
          "KINSOL: failure of the linear solver not reported");
  }
  // refusals
  checkRefusal(mesh, Parameters{{"Strategy", "FixedPoint"}},
               "invalid strategy");
  checkRefusal(mesh, Parameters{{"Damping", 0.5}},
               "requires the strategy 'Picard'");
  checkRefusal(mesh, Parameters{{"MaximumSetupCalls", 0}},
               "maximum number of setup calls");
  checkRefusal(mesh, Parameters{{"ScaledStepTolerance", -1.}},
               "scaled step tolerance");
  checkRefusal(
      mesh,
      Parameters{{"Strategy", "Picard"}, {"AndersonAccelerationDamping", 2.}},
      "damping factors");
  checkRefusal(mesh, Parameters{{"StagnationRelativeTolerance", 1.}},
               "StagnationRelativeTolerance");
  checkRefusal(mesh, Parameters{{"MaximumNumberOfIterations", 0}},
               "maximum number of iterations");
  checkRefusal(mesh, Parameters{{"MaximumNewtonStep", -1.}},
               "maximum Newton step");
  checkRefusal(mesh,
               Parameters{{"Strategy", "Picard"}, {"ScaledStepTolerance", 1.}},
               "is not used by the strategy 'Picard'");
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
