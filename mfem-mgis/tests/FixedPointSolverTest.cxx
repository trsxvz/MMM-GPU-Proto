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
 * preconditioned by BoomerAMG. Newton's method with a relative tolerance of
 * the linear solver adapted to the decrease of the residual (Eisenstat and
 * Walker) gives the same solution with fewer iterations of the linear
 * solver.
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
#include "mfem/general/device.hpp"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearSolvers/AdaptiveLinearSolverTolerance.hxx"

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
  /*!
   * \brief norms of the residuals of the linear solves, if computed (first
   * choice of Eisenstat and Walker)
   */
  std::vector<mfem_mgis::real> linear_residual_norms;
  //! \brief norms of the residual during the last time step
  std::vector<mfem_mgis::real> norms;
  //! \brief information on an iteration
  struct Iteration {
    //! \brief norm of the residual at the end of the iteration
    mfem_mgis::real norm = 0;
    //! \brief relative tolerance of the linear solver, if adapted
    std::optional<mfem_mgis::real> tolerance;
    //! \brief norm of the residual of the linear solve, if computed
    std::optional<mfem_mgis::real> linear_residual_norm;
  };
  //! \brief information on the iterations of each time step
  std::vector<std::vector<Iteration>> iterations_information;
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
  //! \brief options of the linear solver of the problem
  mfem_mgis::Parameters linear_solver_options = {};
  //! \brief if true, the tangent operator is applied without assembly
  bool matrix_free = false;
  //! \brief options of the non linear solver changed before a time step
  mfem_mgis::Parameters changes = {};
  //! \brief index of the time step before which the options are changed
  int changes_step = 1;
  //! \brief number of time steps
  int number_of_steps = 2;
  //! \brief parameters of the non linear solver replacing the default ones
  mfem_mgis::Parameters solver_parameters = {};
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
  problem.setLinearSolver(ctx, ls, opts.linear_solver_options) | or_die;
  if (opts.prediction.has_value()) {
    problem.setPredictionPolicy(
        {.strategy = *(opts.prediction),
         .prediction_operator = PredictionOperator::ELASTIC});
  }
  auto solver_parameters = Parameters{{"VerbosityLevel", 0},
                                      {"RelativeTolerance", 1e-12},
                                      {"AbsoluteTolerance", 0.},
                                      {"MaximumNumberOfIterations", 1000}};
  for (const auto& [k, v] : opts.solver_parameters) {
    solver_parameters.replaceOrInsert(k, v);
  }
  problem.setSolverParameters(ctx, solver_parameters) | or_die;
  auto r = Results{};
  const auto nsteps = opts.number_of_steps;
  const auto dt = real{1} / nsteps;
  auto t = real{0};
  for (int i = 0; i != nsteps; ++i) {
    if ((i == opts.changes_step) && (!opts.changes.empty())) {
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
    auto& step = r.iterations_information.emplace_back();
    for (const auto& information : o.iterations_information) {
      const auto& p = get<Parameters>(throwing, information);
      auto& iteration = step.emplace_back();
      if (contains(p, "Norm")) {
        r.norms.push_back(get<double>(throwing, p, "Norm"));
        iteration.norm = r.norms.back();
      }
      if (contains(p, "LinearSolver")) {
        const auto& ls_information =
            get<Parameters>(throwing, p, "LinearSolver");
        r.linear_iterations += static_cast<size_type>(
            get<int>(throwing, ls_information, "NumberOfIterations"));
        if (contains(ls_information, "ResidualNorm")) {
          r.linear_residual_norms.push_back(
              get<double>(throwing, ls_information, "ResidualNorm"));
          iteration.linear_residual_norm = r.linear_residual_norms.back();
        }
        if (contains(ls_information, "RelativeTolerance")) {
          iteration.tolerance =
              get<double>(throwing, ls_information, "RelativeTolerance");
        }
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

/*!
 * \brief checks the norms of the residuals of the linear solves: given,
 * positive and finite with the first choice of Eisenstat and Walker, not
 * given otherwise
 * \param[in] r: results
 * \param[in] first_choice: if true, the first choice of Eisenstat and Walker
 * is used
 * \param[in] name: name of the computation
 */
static void checkLinearResidualNorms(const Results& r,
                                     const bool first_choice,
                                     const std::string& name) {
  if (!first_choice) {
    check(r.linear_residual_norms.empty(),
          name + ": norms of the residuals of the linear solves given");
    return;
  }
  check((!r.linear_residual_norms.empty()) &&
            std::ranges::all_of(r.linear_residual_norms,
                                [](const mfem_mgis::real v) noexcept {
                                  return std::isfinite(v) && (v > 0);
                                }),
        name + ": norms of the residuals of the linear solves");
}  // end of checkLinearResidualNorms

/*!
 * \brief checks the relative tolerances of the linear solves adapted to the
 * decrease of the residual: the tolerance of the iteration `k >= 2` of each
 * time step is computed again from the norms of the residual of the
 * iterations `k - 1` and `k - 2`, the norm of the residual of the linear
 * solve and the tolerance of the iteration `k - 1` (formulas of Eisenstat
 * and Walker, safeguards); the norm of the residual of a linear solve is
 * lower than the norm of the residual of the non linear problem
 * \param[in] r: results
 * \param[in] t: adaptive tolerance
 * \param[in] goal: absolute convergence criterion of the non linear solver
 * \param[in] name: name of the computation
 */
static void checkAdaptiveTolerances(
    const Results& r,
    const mfem_mgis::AdaptiveLinearSolverTolerance& t,
    const mfem_mgis::real goal,
    const std::string& name) {
  using Method = mfem_mgis::AdaptiveLinearSolverTolerance::Method;
  const auto golden_ratio = (1 + std::sqrt(5.)) / 2;
  const auto first_choice = t.method == Method::EISENSTAT_WALKER_1;
  auto n = 0;
  for (const auto& step : r.iterations_information) {
    for (std::size_t k = 1; k < step.size(); ++k) {
      // norm of the residual before the linear solve of the iteration k
      const auto fk = step[k - 1].norm;
      // the conjugate gradients of the test (without preconditioner, or
      // preconditioned by the exact inverse of the elastic operator) reduce
      // the norm of the residual: a norm above the one of the non linear
      // problem reveals a wrong computation of the residual of the linear
      // solve (its sign, for instance)
      if (step[k].linear_residual_norm.has_value()) {
        check(*(step[k].linear_residual_norm) < fk,
              name +
                  ": norm of the residual of a linear solve above the one "
                  "of the non linear problem");
      }
      if ((k < 2) || (!step[k].tolerance.has_value()) ||
          (!step[k - 1].tolerance.has_value())) {
        continue;
      }
      const auto fkm1 = step[k - 2].norm;
      const auto previous = *(step[k - 1].tolerance);
      auto eta = 0.;
      auto safeguard = 0.;
      if (first_choice) {
        eta =
            std::abs(fk - step[k - 1].linear_residual_norm.value_or(-1)) / fkm1;
        safeguard = std::pow(previous, golden_ratio);
      } else {
        eta = t.gamma * std::pow(fk / fkm1, t.alpha);
        safeguard = t.gamma * std::pow(previous, t.alpha);
      }
      if (safeguard > t.threshold) {
        eta = std::max(eta, safeguard);
      }
      eta = std::min(std::max(eta, goal / (2 * fk)), t.maximum_tolerance);
      const auto v = *(step[k].tolerance);
      check(std::abs(v - eta) <= 1e-12 * eta,
            name + ": relative tolerance of a linear solve (" +
                std::to_string(v) + " instead of " + std::to_string(eta) + ")");
      ++n;
    }
  }
  check(n > 0, name + ": no relative tolerance checked");
}  // end of checkAdaptiveTolerances

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
 * \param[in] message: if not null, part of the expected error message
 */
static void checkRefusal(const char* const mesh,
                         const mfem_mgis::Parameters& options,
                         const std::string& name,
                         const char* const message = nullptr) {
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
  if (message != nullptr) {
    // reading the error message clears it
    const auto e = ctx.getErrorMessage();
    check(e.find(message) != std::string::npos,
          name + ": unexpected error message (" + e + ")");
  }
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
  // results of the Newton solver with the prediction and the number of time
  // steps of a computation
  auto newton_references =
      std::map<std::pair<std::optional<PredictionStrategy>, int>, Results>{};
  auto getNewtonReference = [&ctx, mesh, library, &newton_references](
                                const ComputationOptions& o) -> const Results& {
    const auto key = std::pair{o.prediction, o.number_of_steps};
    auto p = newton_references.find(key);
    if (p == newton_references.end()) {
      const auto r = compute(
          ctx, mesh, library, Parameters{},
          {.prediction = o.prediction, .number_of_steps = o.number_of_steps});
      check(r.status, "Newton: resolution failed");
      p = newton_references.emplace(key, r).first;
    }
    return p->second;
  };
  declareCountingLinearSolver();
  // absolute convergence criterion of the computations with an adaptive
  // tolerance of the linear solver, about 1e-12 times the initial residual,
  // so that the tolerances can be computed again (see
  // `checkAdaptiveTolerances`)
  constexpr auto adaptive_goal = 1e-3;
  const auto adaptive_solver_parameters = Parameters{
      {"RelativeTolerance", 0.}, {"AbsoluteTolerance", adaptive_goal}};
  // adaptive tolerance of a method, with its default options
  auto getAdaptiveTolerance = [](const std::string& m) {
    auto c = Context{};
    const auto ot = AdaptiveLinearSolverTolerance::read(
        c, Parameters{{"AdaptiveLinearSolverTolerance", m}});
    if (!ot.has_value()) {
      std::cerr << c.getErrorMessage() << '\n';
      std::exit(EXIT_FAILURE);
    }
    return *ot;
  };
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
  // NewtonK0 with a relative tolerance of the conjugate gradient adapted to
  // the decrease of the residual (Eisenstat and Walker): fewer iterations of
  // the conjugate gradient (checked below)
  for (const auto* const m : {"EisenstatWalker1", "EisenstatWalker2"}) {
    auto options = newton_k0("CGSolver");
    options.insert(throwing, {{"AdaptiveLinearSolverTolerance", m}});
    cases.push_back({std::string{"NewtonK0, conjugate gradient, "} + m,
                     options,
                     true,
                     {.solver_parameters = adaptive_solver_parameters}});
  }
  // NewtonK0 with the prediction by the elastic operator, the linear solver
  // of the elastic operator being replaced before the third time step: the
  // preconditioners of the tangent operator and of the predictions apply
  // the new one. The conjugate gradient makes the predictions of the second
  // time step (the first one precedes the elastic operator) and of the
  // fourth one (the elastic operator is computed again at the third one).
  // A call of the replaced, destroyed, linear solver would be an undefined
  // behaviour: in practice, the test crashes (checked by mutation).
  {
    auto options = newton_k0("CGSolver");
    options.insert(throwing, {{"PredictionLinearSolver", "ElasticOperator"}});
    cases.push_back(
        {"NewtonK0, then another linear solver of the elastic operator",
         options,
         true,
         {.prediction = bots,
          .changes = {{"ElasticOperatorLinearSolver",
                       parallel ? "HypreBoomerAMG" : direct_solver}},
          .changes_step = 2,
          .number_of_steps = 4},
         0});
  }
  // NewtonK0 with the default iterative solver of the tangent operator (the
  // conjugate gradient, relative tolerance 1e-4)
  cases.push_back({"NewtonK0, default iterative solver of the tangent operator",
                   {{"Method", "NewtonK0"},
                    {"ElasticOperatorUpdate", "Once"},
                    {"ElasticOperatorLinearSolver", direct_solver}},
                   true});
  // prediction by GMRES preconditioned by the linear solver of the elastic
  // operator: one iteration of GMRES, the elastic operator being the one of
  // the prediction
  cases.push_back({"fixed point, Anderson, prediction by GMRESSolver",
                   {{"Acceleration", "Anderson"},
                    {"ElasticOperatorUpdate", "Once"},
                    {"ElasticOperatorLinearSolver", direct_solver},
                    {"PredictionLinearSolver", "GMRESSolver"}},
                   true,
                   {.prediction = bots}});
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
      {"NewtonK0, default iterative solver of the tangent operator", 12},
      {"fixed point, Anderson, prediction by GMRESSolver", 4},
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
  // numbers of iterations of the linear solvers of the cases. For the
  // iterative solvers of the tangent operator of the `NewtonK0` method, the
  // stationary iterations, converging as the fixed point on the elastic
  // operator, make many more iterations than the conjugate gradient, and the
  // adaptive tolerance reduces the number of iterations of the conjugate
  // gradient
  auto linear_iterations = std::map<std::string, size_type>{};
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
    linear_iterations[c.name] = r.linear_iterations;
    checkLinearResidualNorms(r, c.name.ends_with("EisenstatWalker1"), c.name);
    for (const auto* const m : {"EisenstatWalker1", "EisenstatWalker2"}) {
      if (c.name.ends_with(m)) {
        checkAdaptiveTolerances(r, getAdaptiveTolerance(m), adaptive_goal,
                                c.name);
      }
    }
    const auto& newton = getNewtonReference(c.computation);
    const auto du = getRelativeDifference(newton.u, r.u);
    std::cout << c.name << ": " << r.iterations << " iterations ("
              << r.linear_iterations
              << " of the linear solver), relative difference to Newton " << du
              << '\n';
    check(r.status, c.name + ": resolution failed");
    check(du < 1e-8, c.name + ": solution");
    if (const auto p = maximum_numbers_of_iterations.find(c.name);
        p != maximum_numbers_of_iterations.end()) {
      // on a GPU, the sums computed by atomic operations (assembly of the
      // tangent operator, setup of BoomerAMG by hypre) are not deterministic,
      // and the accelerations amplify these differences (up to 45% more
      // iterations observed with Crossed2Deltabis): the bounds are increased
      // by 50%
      const auto gpu = mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                                            mfem::Backend::HIP_MASK);
      const auto bound = gpu ? (3 * p->second) / 2 : p->second;
      check(r.iterations <= bound,
            c.name + ": number of iterations (" + std::to_string(r.iterations) +
                ", at most " + std::to_string(bound) + " expected)");
    }
    if (c.accelerated) {
      check(r.iterations < none.iterations, c.name + ": number of iterations");
    }
    if (c.newton_iterations_from.has_value()) {
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
    const auto cg = linear_iterations.at(
        "Newton, conjugate gradient preconditioned by the elastic operator");
    const auto sli = linear_iterations.at(
        "Newton, SLISolver preconditioned by the elastic operator");
    check((cg > 0) && (sli > 2 * cg),
          "iterations of the iterative solvers of the tangent operator (" +
              std::to_string(cg) + " for the conjugate gradient, " +
              std::to_string(sli) + " for the stationary iterations)");
    for (const auto* const m : {"EisenstatWalker1", "EisenstatWalker2"}) {
      const auto n = linear_iterations.at(
          std::string{"NewtonK0, conjugate gradient, "} + m);
      check(n < cg, std::string{"NewtonK0, conjugate gradient, "} + m +
                        ": iterations of the conjugate gradient (" +
                        std::to_string(n) + ", " + std::to_string(cg) +
                        " with a constant tolerance)");
    }
  }
  // Newton's method with the conjugate gradient of MFEM as linear solver,
  // whose relative tolerance is adapted to the decrease of the residual
  // (Eisenstat and Walker): same solution, fewer iterations of the
  // conjugate gradient than with a constant tolerance
  {
    const auto cg = ComputationOptions{
        .linear_solver = "CGSolver",
        .linear_solver_options = {{"RelativeTolerance", 1e-12},
                                  {"AbsoluteTolerance", 0.},
                                  {"MaximumNumberOfIterations", 10000}}};
    const auto constant = compute(ctx, mesh, library, Parameters{}, cg);
    check(constant.status, "Newton, conjugate gradient: resolution failed");
    for (const auto* const m : {"EisenstatWalker1", "EisenstatWalker2"}) {
      const auto name = std::string{"Newton, conjugate gradient, "} + m;
      auto cg_adaptive = cg;
      cg_adaptive.solver_parameters = adaptive_solver_parameters;
      const auto r = compute(
          ctx, mesh, library,
          Parameters{
              {"Newton", Parameters{{"AdaptiveLinearSolverTolerance", m}}}},
          cg_adaptive);
      const auto du = getRelativeDifference(reference.u, r.u);
      std::cout << name << ": " << r.iterations << " iterations ("
                << r.linear_iterations << " of the linear solver, "
                << constant.linear_iterations
                << " with a constant tolerance), relative difference to "
                   "Newton "
                << du << '\n';
      check(r.status, name + ": resolution failed");
      check(du < 1e-8, name + ": solution");
      checkLinearResidualNorms(r, std::string_view{m} == "EisenstatWalker1",
                               name);
      checkAdaptiveTolerances(r, getAdaptiveTolerance(m), adaptive_goal, name);
      check(r.linear_iterations < constant.linear_iterations,
            name + ": number of iterations of the linear solver");
    }
    // the adaptive tolerance requires an iterative linear solver
    const auto direct = compute(
        ctx, mesh, library,
        Parameters{{"Newton", Parameters{{"AdaptiveLinearSolverTolerance",
                                          "EisenstatWalker2"}}}});
    check(!direct.status,
          "Newton, adaptive tolerance with a direct solver: accepted");
  }
  // linear solver of the predictions: with the conjugate gradient
  // preconditioned by the linear solver of the elastic operator, the linear
  // solver of the problem only makes the first prediction (the elastic
  // operator does not exist yet). The predictions at constant gradients
  // are always made by the linear solver of the problem.
  constexpr auto cgip =
      PredictionStrategy::CONSTANT_GRADIENTS_INTEGRATION_PREDICTION;
  for (const auto& [method, pls, strategy, expected] :
       {std::tuple{"ElasticOperator", "Problem", bots, 2},
        std::tuple{"ElasticOperator", "ElasticOperator", bots, 1},
        std::tuple{"ElasticOperator", "ElasticOperator", cgip, 2},
        std::tuple{"NewtonK0", "ElasticOperator", bots, 1}}) {
    CountingLinearSolver::count = 0;
    auto options = Parameters{{"Method", method},
                              {"ElasticOperatorUpdate", "Once"},
                              {"ElasticOperatorLinearSolver", direct_solver},
                              {"PredictionLinearSolver", std::string{pls}}};
    if (std::string_view{method} == "ElasticOperator") {
      options.insert(throwing, {{"Acceleration", "Anderson"}});
    }
    const auto r = compute(
        ctx, mesh, library, makeFixedPointSolver(options),
        {.prediction = strategy, .linear_solver = "CountingLinearSolver"});
    const auto name = std::string{method} + ", predictions " +
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
        Parameters{
            {"CGSolver",
             Parameters{{"Preconditioner",
                         Parameters{{"Name", "OperatorJacobiSmoother"}}}}}});
    checkRefusal(mesh, options,
                 "a preconditioner of the iterative solver of the tangent "
                 "operator",
                 "option 'Preconditioner' not allowed");
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
  checkRefusal(mesh, {{"AdaptiveLinearSolverTolerance", "EisenstatWalker2"}},
               "an adaptive tolerance of the linear solver with the fixed "
               "point on the elastic operator",
               "only used by the 'NewtonK0' method");
  checkRefusal(mesh, {{"ElasticOperatorUpdate", "Never"}},
               "an invalid update of the elastic operator");
  checkRefusal(mesh, {{"ElasticOperatorLinearSolver", "Unknown"}},
               "an invalid linear solver of the elastic operator");
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
