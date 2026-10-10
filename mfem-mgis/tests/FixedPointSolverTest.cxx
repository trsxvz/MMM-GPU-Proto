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
 * refused. Newton's method with the tangent operator inverted by an
 * iterative solver preconditioned by the linear solver of the elastic
 * operator (`NewtonK0`) makes the iterations of the Newton solver, with
 * each iterative solver of MFEM, with an assembled tangent operator or
 * without assembly. The options may be changed between two time steps.
 * With the `--parallel` option, the elastic operator is also solved by one
 * cycle of BoomerAMG per iteration, or by the conjugate gradient of hypre
 * preconditioned by BoomerAMG.
 * \date   10/10/2026
 */

#include <cmath>
#include <map>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <cstddef>
#include <optional>
#include <tuple>
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
  //! \brief numbers of iterations of each time step
  std::vector<mfem_mgis::size_type> iterations_per_step;
  //! \brief total number of iterations of the iterative linear solvers
  mfem_mgis::size_type linear_iterations = 0;
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

//! \brief options of a computation
struct ComputationOptions {
  /*!
   * \brief prediction strategy at the beginning of each time step, if any,
   * with the elastic operator
   */
  std::optional<mfem_mgis::PredictionStrategy> prediction = {};
  /*!
   * \brief linear solver of the problem, MUMPS in parallel and UMFPack in
   * sequential if null
   */
  const char* linear_solver = nullptr;
  //! \brief if true, the tangent operator is applied without assembly
  bool matrix_free = false;
  //! \brief options of the non linear solver changed after the first step
  mfem_mgis::Parameters changes = {};
};

/*!
 * \return the results of the computation
 * \param[in, out] ctx: execution context
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] solver: value of the `NonLinearSolver` parameter, the Newton
 * solver being used if empty
 * \param[in] opts: options of the computation
 */
static Results compute(mfem_mgis::Context& ctx,
                       const char* const mesh,
                       const char* const library,
                       const mfem_mgis::Parameters& solver,
                       const ComputationOptions& opts = {}) {
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
  if (opts.matrix_free) {
    params.insert(throwing, {{"BehaviourIntegration", "Batch"},
                             {"TangentOperatorComputation", "MatrixFree"}});
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
  const auto* const ls = (opts.linear_solver != nullptr)
                             ? opts.linear_solver
                             : (parallel ? "MUMPSSolver" : "UMFPackSolver");
  problem.setLinearSolver(ctx, ls, {}) | or_die;
  if (opts.prediction.has_value()) {
    problem.setPredictionPolicy(
        {.strategy = *(opts.prediction),
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
    if ((i == 1) && (!opts.changes.empty())) {
      problem.setSolverParameters(ctx, opts.changes) | or_die;
    }
    auto c = Context{};
    const auto o = problem.solve(c, t, dt);
    if (!o.status) {
      std::cerr << c.getErrorMessage() << '\n';
      return r;
    }
    r.iterations += o.iterations;
    r.iterations_per_step.push_back(o.iterations);
    r.norms.clear();
    for (const auto& information : o.iterations_information) {
      const auto& p = get<Parameters>(throwing, information);
      if (contains(p, "Norm")) {
        r.norms.push_back(get<double>(throwing, p, "Norm"));
      }
      if (contains(p, "LinearSolver")) {
        const auto& ls_information =
            get<Parameters>(throwing, p, "LinearSolver");
        r.linear_iterations += static_cast<size_type>(
            get<int>(throwing, ls_information, "NumberOfIterations"));
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
  constexpr auto bots = PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION;
  const auto reference_with_prediction =
      compute(ctx, mesh, library, Parameters{}, {.prediction = bots});
  check(reference_with_prediction.status,
        "Newton with prediction: resolution failed");
  // direct solver, used as linear solver of the elastic operator
  const auto* const direct_solver = parallel ? "MUMPSSolver" : "UMFPackSolver";
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
    // options of the computation
    ComputationOptions computation = {};
    /*!
     * if set, the number of iterations from this time step must be the one
     * of the Newton solver, up to one (the last residuals are close to the
     * tolerance, 1e-12, and depend on the rounding errors of the linear
     * solvers)
     */
    std::optional<std::size_t> newton_iterations_from = {};
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
        {"ElasticOperatorLinearSolver", direct_solver}},
       true},
      {"fixed point, Anderson, elastic operator by integration",
       {{"Acceleration", "Anderson"},
        {"ElasticOperatorComputation", "Integration"}},
       true},
      {"modified Newton", {{"Method", "ModifiedNewton"}}, true},
      {"modified Newton, Anderson",
       {{"Method", "ModifiedNewton"}, {"Acceleration", "Anderson"}},
       true}};
  // Newton's method, the tangent operator being inverted by an iterative
  // solver preconditioned by the linear solver of the elastic operator
  auto newton_k0 = [direct_solver](const char* const s) {
    return Parameters{
        {"Method", "NewtonK0"},
        {"ElasticOperatorUpdate", "Once"},
        {"ElasticOperatorLinearSolver", direct_solver},
        {"TangentLinearSolver",
         Parameters{{s, Parameters{{"RelativeTolerance", 1e-12},
                                   {"MaximumNumberOfIterations", 2000}}}}}};
  };
  cases.push_back(
      {"Newton, conjugate gradient preconditioned by the elastic operator",
       newton_k0("CGSolver"),
       true,
       {},
       0});
  // the stationary iterations (`SLISolver`) converge as the fixed point on
  // the elastic operator
  for (const auto* const s : {"GMRESSolver", "FGMRESSolver", "BiCGSTABSolver",
                              "MINRESSolver", "SLISolver"}) {
    cases.push_back({std::string{"Newton, "} + s +
                         " preconditioned by the elastic operator",
                     newton_k0(s),
                     true,
                     {},
                     0});
  }
  cases.push_back(
      {"Newton, conjugate gradient preconditioned by the elastic "
       "operator, tangent operator applied without assembly",
       newton_k0("CGSolver"),
       true,
       {.matrix_free = true},
       0});
  {
    auto options = newton_k0("CGSolver");
    options.insert(throwing, {{"PredictionLinearSolver", "ElasticOperator"}});
    cases.push_back(
        {"Newton, conjugate gradient preconditioned by the elastic "
         "operator, prediction by the elastic operator",
         options,
         true,
         {.prediction = bots},
         0});
  }
  if (parallel) {
    auto options = newton_k0("CGSolver");
    options.replaceOrInsert("ElasticOperatorLinearSolver", "HypreBoomerAMG");
    cases.push_back(
        {"Newton, conjugate gradient preconditioned by BoomerAMG "
         "built on the elastic operator",
         options,
         true,
         {},
         0});
  }
  // options changed between the two time steps: from the accelerated fixed
  // point to NewtonK0 (the acceleration being removed) and back, and from
  // the linear solver of the elastic operator to the one of the problem
  {
    auto changes = newton_k0("CGSolver");
    changes.insert(throwing, {{"Acceleration", "None"}});
    cases.push_back({"fixed point, Anderson, then NewtonK0",
                     {{"Acceleration", "Anderson"},
                      {"ElasticOperatorUpdate", "Once"},
                      {"ElasticOperatorLinearSolver", direct_solver}},
                     true,
                     {.changes = changes},
                     1});
  }
  cases.push_back({"NewtonK0, then fixed point, Anderson",
                   newton_k0("CGSolver"),
                   true,
                   {.changes = {{"Method", "ElasticOperator"},
                                {"Acceleration", "Anderson"}}}});
  cases.push_back(
      {"fixed point, Anderson, own linear solver, then the "
       "linear solver of the problem",
       {{"Acceleration", "Anderson"},
        {"ElasticOperatorLinearSolver", direct_solver}},
       true,
       {.changes = {{"ElasticOperatorLinearSolver", "None"}}}});
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
          {"ElasticOperatorLinearSolver", direct_solver},
          {"PredictionLinearSolver", std::string{pls}}},
         true,
         {.prediction = bots}});
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
  // numbers of iterations of the iterative solvers of the tangent operator
  // of the `NewtonK0` method: the stationary iterations, converging as the
  // fixed point on the elastic operator, make many more iterations than
  // the conjugate gradient
  auto tangent_linear_iterations = std::map<std::string, size_type>{};
  for (const auto& c : cases) {
    const auto r = compute(ctx, mesh, library, makeFixedPointSolver(c.options),
                           c.computation);
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
    if (c.newton_iterations_from == 0) {
      tangent_linear_iterations[c.name] = r.linear_iterations;
    }
    const auto du = getRelativeDifference(reference.u, r.u);
    std::cout << c.name << ": " << r.iterations << " iterations ("
              << r.linear_iterations
              << " of the linear solver), relative difference to Newton " << du
              << '\n';
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
    if (c.newton_iterations_from.has_value()) {
      const auto& newton = c.computation.prediction.has_value()
                               ? reference_with_prediction
                               : reference;
      // numbers of iterations from the given time step
      auto count = [&c](const Results& results) {
        const auto& n = results.iterations_per_step;
        if (n.size() <= *(c.newton_iterations_from)) {
          return size_type{};
        }
        return std::accumulate(n.begin() + static_cast<std::ptrdiff_t>(
                                               *(c.newton_iterations_from)),
                               n.end(), size_type{});
      };
      const auto ni = count(r);
      const auto nr = count(newton);
      check((ni > 0) && (ni <= nr + 1),
            c.name + ": number of iterations (" + std::to_string(ni) + ", " +
                std::to_string(nr) + " with the Newton solver)");
    }
  }
  {
    const auto cg = tangent_linear_iterations.at(
        "Newton, conjugate gradient preconditioned by the elastic operator");
    const auto sli = tangent_linear_iterations.at(
        "Newton, SLISolver preconditioned by the elastic operator");
    check((cg > 0) && (sli > 2 * cg),
          "iterations of the iterative solvers of the tangent operator (" +
              std::to_string(cg) + " for the conjugate gradient, " +
              std::to_string(sli) + " for the stationary iterations)");
  }
  // linear solver of the predictions: with the conjugate gradient
  // preconditioned by the linear solver of the elastic operator, the linear
  // solver of the problem only makes the first prediction (the elastic
  // operator does not exist yet). The predictions at constant gradients
  // are always made by the linear solver of the problem.
  declareCountingLinearSolver();
  constexpr auto cgip =
      PredictionStrategy::CONSTANT_GRADIENTS_INTEGRATION_PREDICTION;
  for (const auto& [pls, strategy, expected] :
       {std::tuple{"Problem", bots, 2}, std::tuple{"ElasticOperator", bots, 1},
        std::tuple{"ElasticOperator", cgip, 2}}) {
    CountingLinearSolver::count = 0;
    const auto r = compute(
        ctx, mesh, library,
        makeFixedPointSolver({{"Acceleration", "Anderson"},
                              {"ElasticOperatorUpdate", "Once"},
                              {"ElasticOperatorLinearSolver", direct_solver},
                              {"PredictionLinearSolver", std::string{pls}}}),
        {.prediction = strategy, .linear_solver = "CountingLinearSolver"});
    const auto name = std::string{"predictions "} +
                      (strategy == cgip ? "at constant gradients " : "") +
                      "by the linear solver '" + pls + "'";
    check(r.status && (getRelativeDifference(reference.u, r.u) < 1e-8),
          name + ": solution");
    check(CountingLinearSolver::count == expected,
          name + ": " + std::to_string(CountingLinearSolver::count) +
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
                {"ElasticOperatorLinearSolver", direct_solver}},
               "a linear solver of the elastic operator with the modified "
               "Newton method");
  checkRefusal(mesh, {{"Method", "NewtonK0"}},
               "the NewtonK0 method without linear solver of the elastic "
               "operator");
  {
    auto options = newton_k0("CGSolver");
    options.insert(throwing, {{"Acceleration", "Anderson"}});
    checkRefusal(mesh, options, "an accelerated NewtonK0 method");
  }
  checkRefusal(mesh, {{"PredictionLinearSolver", "Tangent"}},
               "an invalid linear solver of the predictions");
  checkRefusal(mesh, {{"PredictionLinearSolver", direct_solver}},
               "a direct solver as iterative solver of the predictions");
  checkRefusal(
      mesh,
      {{"PredictionLinearSolver",
        Parameters{{"Problem", Parameters{{"RelativeTolerance", 1e-8}}}}}},
      "options of the linear solver of the problem for the predictions");
  {
    auto options = newton_k0("CGSolver");
    options.replaceOrInsert("TangentLinearSolver", direct_solver);
    checkRefusal(mesh, options,
                 "a direct solver as iterative solver of the tangent "
                 "operator");
    options.replaceOrInsert(
        "TangentLinearSolver",
        Parameters{{"CGSolver",
                    Parameters{{"Preconditioner",
                                Parameters{{"Name", "HypreBoomerAMG"}}}}}});
    checkRefusal(mesh, options,
                 "a preconditioner of the iterative solver of the tangent "
                 "operator");
  }
  checkRefusal(
      mesh, {{"Acceleration", Parameters{{"None", Parameters{{"Depth", 4}}}}}},
      "options of the removal of the acceleration");
  checkRefusal(
      mesh,
      {{"ElasticOperatorLinearSolver",
        Parameters{{"Name", "None"},
                   {"Options", Parameters{{"Symmetric", true}}}}}},
      "options of the removal of the linear solver of the elastic operator");
  checkRefusal(
      mesh,
      {{"Acceleration",
        Parameters{{"Anderson", Parameters{{"Regularization", 1e-20}}}}}},
      "a regularization of Anderson's method below the rounding "
      "errors");
  checkRefusal(mesh, {{"ElasticOperatorUpdate", "Never"}},
               "an invalid update of the elastic operator");
  checkRefusal(mesh, {{"ElasticOperatorLinearSolver", "Unknown"}},
               "an invalid linear solver of the elastic operator");
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
