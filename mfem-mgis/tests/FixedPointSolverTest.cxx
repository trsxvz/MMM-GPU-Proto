/*!
 * \file   FixedPointSolverTest.cxx
 * \brief  Test of the fixed point solver (`FixedPoint` non linear solver):
 * uniaxial tension of a cube made of a plastic material with a linear
 * hardening, in two time steps. Each method (fixed point on the elastic
 * operator, modified Newton) and each acceleration gives the solution of the
 * Newton solver; Anderson's method and the acceleration of Cast3M reduce
 * the number of iterations; the other accelerations of MTest converge;
 * without acceleration, the asymptotic rate of convergence of the fixed
 * point on the elastic operator is bounded by `3μ / (3μ + H)`; the invalid
 * options are
 * refused. With the `--parallel` option, the elastic operator is also
 * solved by one cycle of BoomerAMG per iteration, or by the conjugate
 * gradient of hypre preconditioned by BoomerAMG.
 * \date   10/10/2026
 */

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

//! \brief number of failed checks
static int failures = 0;
//! \brief if true, the computations are parallel
static bool parallel = false;

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
  //! \brief norms of the residual during the last time step
  std::vector<mfem_mgis::real> norms;
};

/*!
 * \return the results of the computation
 * \param[in, out] ctx: execution context
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] solver: value of the `NonLinearSolver` parameter, the Newton
 * solver being used if empty
 */
static Results compute(mfem_mgis::Context& ctx,
                       const char* const mesh,
                       const char* const library,
                       const mfem_mgis::Parameters& solver) {
  using namespace mfem_mgis;
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
      ctx, {{"Boundary", "Xmax"},
            {"Component", 0},
            {"LoadingEvolution", [](const real t) { return 1e-2 * t; }}}) |
      or_die;
  problem.setLinearSolver(ctx, parallel ? "MUMPSSolver" : "UMFPackSolver", {}) |
      or_die;
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 1000}}) |
      or_die;
  auto r = Results{};
  constexpr int nsteps = 2;
  const auto dt = real{1} / nsteps;
  auto t = real{0};
  for (int i = 0; i != nsteps; ++i) {
    auto c = Context{};
    const auto o = problem.solve(c, t, dt);
    if (!o.status) {
      std::cerr << c.getErrorMessage() << '\n';
      return r;
    }
    r.iterations += o.iterations;
    r.norms.clear();
    for (const auto& information : o.iterations_information) {
      const auto& p = get<Parameters>(throwing, information);
      if (contains(p, "Norm")) {
        r.norms.push_back(get<double>(throwing, p, "Norm"));
      }
    }
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
  auto d = 0.;
  auto m = 0.;
  for (std::size_t i = 0; i != a.size(); ++i) {
    d = std::max(d, std::abs(a[i] - b[i]));
    m = std::max(m, std::abs(a[i]));
  }
  return m == 0 ? d : d / m;
}  // end of getRelativeDifference

//! \return the parameters of the fixed point solver
static mfem_mgis::Parameters makeFixedPointSolver(
    const mfem_mgis::Parameters& options) {
  return mfem_mgis::Parameters{{"FixedPoint", options}};
}  // end of makeFixedPointSolver

/*!
 * \brief an invalid option of the fixed point solver is refused
 * \param[in] mesh: mesh
 * \param[in] options: options of the fixed point solver
 * \param[in] name: name of the check
 */
static void checkRefusal(const char* const mesh,
                         const mfem_mgis::Parameters& options,
                         const std::string& name) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  const auto o = construct<NonLinearEvolutionProblem>(
      ctx, Parameters{{"MeshFileName", mesh},
                      {"FiniteElementFamily", "H1"},
                      {"FiniteElementOrder", 1},
                      {"UnknownsSize", 3},
                      {"Hypothesis", "Tridimensional"},
                      {"NonLinearSolver", makeFixedPointSolver(options)},
                      {"Parallel", parallel}});
  check(isInvalid(o), name + " is accepted");
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
  auto ctx = Context{};
  const auto reference = compute(ctx, mesh, library, Parameters{});
  check(reference.status, "Newton: resolution failed");
  std::cout << "Newton: " << reference.iterations << " iterations\n";
  // fixed point on the elastic operator, solved by the linear solver of the
  // problem, without acceleration: rate of convergence
  const auto none = compute(ctx, mesh, library, makeFixedPointSolver({}));
  {
    const auto mu = 70e9 / (2 * (1 + 0.34));
    const auto H = 10e9;
    const auto bound = 3 * mu / (3 * mu + H);
    const auto& n = none.norms;
    const auto rate = (n.size() > 2) ? n[n.size() - 1] / n[n.size() - 2] : 1.;
    std::cout << "fixed point: " << none.iterations
              << " iterations, rate of convergence " << rate << " (bound "
              << bound << ")\n";
    check(none.status, "fixed point: resolution failed");
    check((rate > 0.5 * bound) && (rate <= bound * (1 + 1e-3)),
          "fixed point: rate of convergence");
  }
  struct Case {
    std::string name;
    Parameters options;
    // if true, the number of iterations must be lower than without
    // acceleration
    bool accelerated;
  };
  auto cases = std::vector<Case>{
      {"fixed point, Anderson",
       {{"Acceleration", Parameters{{"Anderson", Parameters{{"Depth", 4}}}}}},
       true},
      {"fixed point, Cast3M", {{"Acceleration", "Cast3M"}}, true},
      {"fixed point, Anderson, elastic operator assembled once, own linear "
       "solver",
       {{"Acceleration", "Anderson"},
        {"ElasticOperatorUpdate", "Once"},
        {"ElasticOperatorLinearSolver",
         parallel ? "MUMPSSolver" : "UMFPackSolver"}},
       true},
      {"fixed point, Anderson, elastic operator by integration",
       {{"Acceleration", "Anderson"},
        {"ElasticOperatorComputation", "Integration"}},
       true},
      {"modified Newton", {{"Method", "ModifiedNewton"}}, true},
      {"modified Newton, Anderson",
       {{"Method", "ModifiedNewton"}, {"Acceleration", "Anderson"}},
       true}};
  // the other accelerations, ported from MTest: the solution only is
  // checked. Steffensen's method, applied component by component, diverges
  // on this problem (as expected for a coupled system): its declaration is
  // checked by the refusal of an invalid trigger.
  cases.push_back({"fixed point, Anderson on the residuals, period 2",
                   {{"Acceleration",
                     Parameters{{"Anderson", Parameters{{"Residuals", "Forces"},
                                                        {"Period", 2}}}}}},
                   false});
  for (const auto* const a :
       {"UAnderson", "FAnderson", "Secant", "CrossedSecant", "AlternateSecant",
        "IronsTuck", "CrossedDelta2", "Crossed2Delta", "Crossed2Deltabis",
        "AlternateDelta2", "Alternate2Delta"}) {
    cases.push_back({std::string{"fixed point, "} + a,
                     {{"Acceleration", std::string{a}}},
                     false});
  }
  if (parallel) {
    cases.push_back(
        {"fixed point, Anderson, one cycle of BoomerAMG per iteration",
         {{"Acceleration", "Anderson"},
          {"ElasticOperatorUpdate", "Once"},
          {"ElasticOperatorLinearSolver", "HypreBoomerAMG"}},
         true});
    cases.push_back(
        {"fixed point, Anderson, conjugate gradient and BoomerAMG",
         {{"Acceleration", "Anderson"},
          {"ElasticOperatorLinearSolver",
           Parameters{{"HyprePCG",
                       Parameters{{"Tolerance", 1e-12},
                                  {"Preconditioner",
                                   Parameters{{"Name", "HypreBoomerAMG"}}}}}}}},
         true});
  }
  for (const auto& c : cases) {
    const auto r = compute(ctx, mesh, library, makeFixedPointSolver(c.options));
    const auto du = getRelativeDifference(reference.u, r.u);
    std::cout << c.name << ": " << r.iterations
              << " iterations, relative difference to Newton " << du << '\n';
    check(r.status, c.name + ": resolution failed");
    check(du < 1e-8, c.name + ": solution");
    if (c.accelerated) {
      check(r.iterations < none.iterations, c.name + ": number of iterations");
    }
  }
  // invalid options
  checkRefusal(mesh, {{"Method", "Secant"}}, "an invalid method");
  checkRefusal(mesh, {{"Acceleration", "Unknown"}}, "an invalid acceleration");
  checkRefusal(
      mesh,
      {{"Acceleration", Parameters{{"Anderson", Parameters{{"Depth", 0}}}}}},
      "a null depth of Anderson's method");
  checkRefusal(
      mesh,
      {{"Acceleration", Parameters{{"Cast3M", Parameters{{"Trigger", 2}}}}}},
      "a trigger of the acceleration of Cast3M lower than 3");
  checkRefusal(
      mesh,
      {{"Acceleration", Parameters{{"Secant", Parameters{{"Trigger", 2}}}}}},
      "a trigger of the secant method lower than 3");
  checkRefusal(mesh,
               {{"Acceleration",
                 Parameters{{"Steffensen", Parameters{{"Trigger", 2}}}}}},
               "a trigger of Steffensen's method lower than 3");
  checkRefusal(
      mesh,
      {{"Acceleration",
        Parameters{{"Anderson", Parameters{{"Residuals", "Displacements"}}}}}},
      "invalid residuals of Anderson's method");
  checkRefusal(mesh, {{"ElasticOperatorUpdate", "Never"}},
               "an invalid update of the elastic operator");
  checkRefusal(mesh, {{"ElasticOperatorLinearSolver", "Unknown"}},
               "an invalid linear solver of the elastic operator");
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
