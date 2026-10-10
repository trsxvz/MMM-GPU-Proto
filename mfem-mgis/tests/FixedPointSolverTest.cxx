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
#include <map>
#include <cstdlib>
#include <memory>
#include <utility>
#include <string>
#include <string_view>
#include <vector>
#include <iostream>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
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
 * \brief linear solver counting its resolutions, used to check which linear
 * solver makes the predictions
 */
struct CountingLinearSolver final : mfem::Solver {
  //! \brief number of resolutions
  static inline int count = 0;
  //! \param[in] s: linear solver
  explicit CountingLinearSolver(std::unique_ptr<mfem::Solver> s)
      : solver(std::move(s)) {}
  void SetOperator(const mfem::Operator& op) override {
    this->height = op.Height();
    this->width = op.Width();
    this->solver->SetOperator(op);
  }
  void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
    ++count;
    this->solver->Mult(x, y);
  }

 private:
  //! \brief linear solver
  std::unique_ptr<mfem::Solver> solver;
};

/*!
 * \brief declares the `CountingLinearSolver`, wrapping MUMPS in parallel
 * and UMFPack in sequential
 */
static void declareCountingLinearSolver() {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  if (parallel) {
#ifdef MFEM_USE_MPI
    LinearSolverFactory<true>::getFactory().add(
        ctx, "CountingLinearSolver",
        [](Context& c, FiniteElementSpace<true>& fes,
           const Parameters& opts) -> LinearSolverHandler {
          auto h = LinearSolverFactory<true>::getFactory().generate(
              c, "MUMPSSolver", fes, opts);
          return {std::make_unique<CountingLinearSolver>(
                      std::move(h.linear_solver)),
                  {}};
        }) |
        or_die;
#endif /* MFEM_USE_MPI */
  } else {
    LinearSolverFactory<false>::getFactory().add(
        ctx, "CountingLinearSolver",
        [](Context& c, FiniteElementSpace<false>& fes,
           const Parameters& opts) -> LinearSolverHandler {
          auto h = LinearSolverFactory<false>::getFactory().generate(
              c, "UMFPackSolver", fes, opts);
          return {std::make_unique<CountingLinearSolver>(
                      std::move(h.linear_solver)),
                  {}};
        }) |
        or_die;
  }
}  // end of declareCountingLinearSolver

/*!
 * \return the results of the computation
 * \param[in, out] ctx: execution context
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] solver: value of the `NonLinearSolver` parameter, the Newton
 * solver being used if empty
 * \param[in] prediction: if true, an elastic prediction is made at the
 * beginning of each time step
 * \param[in] linear_solver: linear solver of the problem, MUMPS in parallel
 * and UMFPack in sequential if null
 */
static Results compute(mfem_mgis::Context& ctx,
                       const char* const mesh,
                       const char* const library,
                       const mfem_mgis::Parameters& solver,
                       const bool prediction = false,
                       const char* const linear_solver = nullptr) {
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
  const auto* const ls = (linear_solver != nullptr)
                             ? linear_solver
                             : (parallel ? "MUMPSSolver" : "UMFPackSolver");
  problem.setLinearSolver(ctx, ls, {}) | or_die;
  if (prediction) {
    problem.setPredictionPolicy(
        {.strategy = PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION,
         .prediction_operator = PredictionOperator::ELASTIC});
  }
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
    // if true, an elastic prediction is made at the beginning of each time
    // step
    bool prediction = false;
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
  // Newton's method, the tangent operator being inverted by the conjugate
  // gradient preconditioned by the linear solver of the elastic operator
  cases.push_back(
      {"Newton, conjugate gradient preconditioned by the elastic operator",
       {{"Method", "NewtonK0"},
        {"ElasticOperatorUpdate", "Once"},
        {"ElasticOperatorLinearSolver",
         parallel ? "MUMPSSolver" : "UMFPackSolver"},
        {"TangentLinearSolver", Parameters{{"RelativeTolerance", 1e-12}}}},
       true});
  if (parallel) {
    cases.push_back(
        {"Newton, conjugate gradient preconditioned by BoomerAMG "
         "built on the elastic operator",
         {{"Method", "NewtonK0"},
          {"ElasticOperatorUpdate", "Once"},
          {"ElasticOperatorLinearSolver", "HypreBoomerAMG"},
          {"TangentLinearSolver", Parameters{{"RelativeTolerance", 1e-12}}}},
         true});
  }
  // elastic prediction, by the linear solver of the problem or by the
  // conjugate gradient preconditioned by the linear solver of the elastic
  // operator
  for (const auto* const pls : {"Problem", "ElasticOperator"}) {
    cases.push_back(
        {std::string{"fixed point, Anderson, prediction by the linear solver "
                     "of the "} +
             (std::string_view{pls} == "Problem" ? "problem"
                                                 : "elastic operator"),
         {{"Acceleration", "Anderson"},
          {"ElasticOperatorUpdate", "Once"},
          {"ElasticOperatorLinearSolver",
           parallel ? "MUMPSSolver" : "UMFPackSolver"},
          {"PredictionLinearSolver", std::string{pls}}},
         true,
         true});
  }
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
        {"fixed point, Anderson, two cycles of BoomerAMG per iteration",
         {{"Acceleration", "Anderson"},
          {"ElasticOperatorUpdate", "Once"},
          {"ElasticOperatorLinearSolver",
           Parameters{
               {"HypreBoomerAMG", Parameters{{"VerbosityLevel", 0},
                                             {"NumberOfCycles", 2},
                                             {"RelaxationSweeps", 1},
                                             {"StrengthThreshold", 0.25}}}}}},
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
  // maximum numbers of iterations of the accelerations, about 15 % above
  // the observed ones, to detect a wrong acceleration which still converges
  const auto maximum_numbers_of_iterations = std::map<std::string, size_type>{
      {"fixed point, Anderson", 43},
      {"Newton, conjugate gradient preconditioned by the elastic operator", 10},
      {"Newton, conjugate gradient preconditioned by BoomerAMG built on the "
       "elastic operator",
       10},
      {"fixed point, Cast3M", 65},
      {"modified Newton", 70},
      {"modified Newton, Anderson", 38},
      {"fixed point, Anderson on the residuals, period 2", 44},
      {"fixed point, UAnderson", 44},
      {"fixed point, FAnderson", 44},
      {"fixed point, Secant", 75},
      {"fixed point, CrossedSecant", 48},
      {"fixed point, AlternateSecant", 71},
      {"fixed point, IronsTuck", 56},
      {"fixed point, CrossedDelta2", 161},
      {"fixed point, Crossed2Delta", 81},
      {"fixed point, Crossed2Deltabis", 431},
      {"fixed point, AlternateDelta2", 69},
      {"fixed point, Alternate2Delta", 53}};
  // number of iterations of the fixed point with one cycle of BoomerAMG per
  // iteration, compared to two cycles per iteration (a setting of BoomerAMG
  // which must survive the setup of the elastic operator)
  auto amg_iterations = size_type{};
  for (const auto& c : cases) {
    const auto r = compute(ctx, mesh, library, makeFixedPointSolver(c.options),
                           c.prediction);
    if (c.name ==
        "fixed point, Anderson, one cycle of BoomerAMG per iteration") {
      amg_iterations = r.iterations;
    }
    if (c.name ==
        "fixed point, Anderson, two cycles of BoomerAMG per iteration") {
      check(r.iterations < amg_iterations,
            c.name + ": number of iterations (" + std::to_string(r.iterations) +
                " with two cycles, " + std::to_string(amg_iterations) +
                " with one cycle)");
    }
    const auto du = getRelativeDifference(reference.u, r.u);
    std::cout << c.name << ": " << r.iterations
              << " iterations, relative difference to Newton " << du << '\n';
    check(r.status, c.name + ": resolution failed");
    check(du < 1e-8, c.name + ": solution");
    if (const auto p = maximum_numbers_of_iterations.find(c.name);
        p != maximum_numbers_of_iterations.end()) {
      check(r.iterations <= p->second,
            c.name + ": number of iterations (" + std::to_string(r.iterations) +
                ", at most " + std::to_string(p->second) + " expected)");
    }
    if (c.accelerated) {
      check(r.iterations < none.iterations, c.name + ": number of iterations");
    }
  }
  // linear solver of the predictions: with the conjugate gradient
  // preconditioned by the linear solver of the elastic operator, the linear
  // solver of the problem only makes the first prediction (the elastic
  // operator does not exist yet)
  declareCountingLinearSolver();
  for (const auto* const pls : {"Problem", "ElasticOperator"}) {
    CountingLinearSolver::count = 0;
    const auto r = compute(
        ctx, mesh, library,
        makeFixedPointSolver({{"Acceleration", "Anderson"},
                              {"ElasticOperatorUpdate", "Once"},
                              {"ElasticOperatorLinearSolver",
                               parallel ? "MUMPSSolver" : "UMFPackSolver"},
                              {"PredictionLinearSolver", std::string{pls}}}),
        true, "CountingLinearSolver");
    const auto expected = std::string_view{pls} == "Problem" ? 2 : 1;
    check(r.status && (getRelativeDifference(reference.u, r.u) < 1e-8),
          std::string{"predictions by the linear solver '"} + pls +
              "': solution");
    check(CountingLinearSolver::count == expected,
          std::string{"predictions by the linear solver '"} + pls +
              "': " + std::to_string(CountingLinearSolver::count) +
              " resolutions by the linear solver of the problem, " +
              std::to_string(expected) + " expected");
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
  checkRefusal(mesh,
               {{"Method", "ModifiedNewton"},
                {"ElasticOperatorLinearSolver",
                 parallel ? "MUMPSSolver" : "UMFPackSolver"}},
               "a linear solver of the elastic operator with the modified "
               "Newton method");
  checkRefusal(mesh, {{"Method", "NewtonK0"}},
               "the NewtonK0 method without linear solver of the elastic "
               "operator");
  checkRefusal(mesh,
               {{"Method", "NewtonK0"},
                {"ElasticOperatorLinearSolver",
                 parallel ? "MUMPSSolver" : "UMFPackSolver"},
                {"Acceleration", "Anderson"}},
               "an accelerated NewtonK0 method");
  checkRefusal(mesh, {{"PredictionLinearSolver", "Tangent"}},
               "an invalid linear solver of the predictions");
  checkRefusal(mesh, {{"ElasticOperatorUpdate", "Never"}},
               "an invalid update of the elastic operator");
  checkRefusal(mesh, {{"ElasticOperatorLinearSolver", "Unknown"}},
               "an invalid linear solver of the elastic operator");
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
