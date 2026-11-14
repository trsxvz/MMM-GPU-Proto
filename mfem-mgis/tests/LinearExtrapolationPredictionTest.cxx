/*!
 * \file   LinearExtrapolationPredictionTest.cxx
 * \brief  Test of the prediction by linear extrapolation of the increment
 * of the unknowns of the previous time step
 * (`PredictionStrategy::LINEAR_EXTRAPOLATION_PREDICTION`).
 *
 * The face `Xmax` of a cube made of the
 * `IsotropicLinearHardeningPlasticity` behaviour is moved along its normal,
 * the faces `Xmin`, `Ymin` and `Zmin` being fixed along their normals: the
 * strain is uniform.
 *
 * In the elastic domain, the problem is linear. A resolution limited to
 * zero iteration leaves the prediction in the unknowns at the end of the
 * time step, and returns the norm of its residual and the reference of the
 * relative criterion. With a displacement which is not linear in time:
 * - at the first time step, the prediction is the default one;
 * - then, the prediction is the extrapolation of the previous increment
 *   weighted by the ratio of the time increments, the values imposed by the
 *   Dirichlet boundary conditions being the ones of the end of the time
 *   step;
 * - the reference of the relative criterion is the one of the default
 *   prediction, and the ratio between the norm of the residual of the
 *   prediction and this reference is the analytical one;
 * - after `revert`, the extrapolation is made again with the new time
 *   increment;
 * - a second resolution of the same time step starts from the solution of
 *   the first one;
 * - the default strategy releases the increment, and the linear
 *   extrapolation selected again does not use a stale one.
 *
 * With a displacement proportional to the time, the extrapolation is the
 * solution: no iteration is made from the second time step, also when the
 * time increment changes, whereas the default prediction requires at least
 * one.
 *
 * In the plastic domain, the solution is the one of the default prediction,
 * with at most as many iterations.
 *
 * With the `FailingIsotropicLinearElasticity` behaviour, whose integration
 * fails where the absolute value of the trace of the strain exceeds a
 * threshold, and a displacement which increases then decreases, the
 * behaviour can't be integrated at the extrapolated estimate of the time
 * step of the reversal, whereas it can be at the default prediction: the
 * resolution starts again from the default prediction, with a warning, and
 * gives the resolution of the default strategy. The extrapolation is used
 * again as soon as the behaviour can be integrated at its estimate. The
 * failures of the integration at the following estimates of a resolution
 * are left to the non linear solver.
 *
 * The checks are made with the Newton solver and with the `FixedPoint`
 * solver (fixed point on the elastic operator, `NewtonK0` method).
 * \date   14/11/2026
 */

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <iostream>
#include <algorithm>
#include <functional>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx"
#include "SequentialLinearSolver.hxx"
#include "Tolerance.hxx"

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

//! \return if two values are equal up to a relative tolerance
static bool areClose(const double a, const double b, const double tolerance) {
  return std::abs(a - b) <= tolerance * std::max(std::abs(a), std::abs(b));
}  // end of areClose

//! \return if two values are the same
static bool areEqual(const mfem_mgis::real a, const mfem_mgis::real b) {
  return (a <= b) && (b <= a);
}  // end of areEqual

/*!
 * \brief displacement, or uniaxial stress, imposed on the face `Xmax` as a
 * function of time
 */
using Loading = std::function<mfem_mgis::real(const mfem_mgis::real)>;

//! \brief non linear solver
struct Solver {
  //! \brief name of the solver in the messages
  std::string name;
  //! \brief value of the `NonLinearSolver` parameter, the default if empty
  mfem_mgis::Parameters parameters;
  //! \brief if true, the solver is a Newton method
  bool newton = true;
};

//! \brief relative tolerance of the non linear resolutions
static constexpr auto relative_tolerance =
    mfem_mgis::unit_tests::tolerance(1e-10, 1e-4);
//! \brief maximum number of iterations of the non linear resolutions
static constexpr auto maximum_number_of_iterations = 200;

//! \brief result of a resolution
struct Resolution {
  //! \brief true if the resolution converged
  bool status = false;
  //! \brief number of iterations
  mfem_mgis::size_type iterations = 0;
  //! \brief reference of the relative criterion
  double reference = 0;
  //! \brief norm of the residual at the end of the resolution
  double norm = 0;
  //! \brief unknowns at the end of the resolution
  std::vector<mfem_mgis::real> u;
  /*!
   * \brief number of warnings stating that the resolution started again
   * from the default prediction (printed by the main process only)
   */
  std::size_t fallbacks = 0;
};

//! \brief warning printed when the extrapolated estimate is discarded
static constexpr auto fallback_warning =
    "the resolution starts again from the default prediction";

/*!
 * \return a non linear evolution problem on a cube
 * \param[in, out] ctx: execution context
 * \param[in] mesh: mesh of a cube
 * \param[in] solver: non linear solver
 */
static mfem_mgis::NonLinearEvolutionProblem makeProblem(mfem_mgis::Context& ctx,
                                                        const char* const mesh,
                                                        const Solver& solver) {
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
  if (!solver.parameters.empty()) {
    params.insert(throwing, {{"NonLinearSolver", solver.parameters}});
  }
  return construct<NonLinearEvolutionProblem>(ctx, params) | or_die;
}  // end of makeProblem

/*!
 * \brief cube whose face `Xmax` is moved along its normal, or loaded by a
 * uniform stress along its normal
 */
struct Cube {
  /*!
   * \brief constructor
   * \param[in] mesh: mesh of a cube
   * \param[in] library: behaviour library
   * \param[in] solver: non linear solver
   * \param[in] loading: displacement imposed on the face `Xmax`, or
   * uniaxial stress
   * \param[in] strategy: prediction strategy
   * \param[in] behaviour: behaviour
   * \param[in] imposed_stress: if true, the loading is the uniaxial stress,
   * imposed by a pressure on the face `Xmax`
   */
  Cube(const char* const mesh,
       const char* const library,
       const Solver& solver,
       const Loading& loading,
       const mfem_mgis::PredictionStrategy strategy,
       const char* const behaviour = "IsotropicLinearHardeningPlasticity",
       const bool imposed_stress = false)
      : problem(makeProblem(ctx, mesh, solver)) {
    using namespace mfem_mgis;
    auto or_die = this->ctx.getFatalFailureHandler();
    this->problem.addBehaviourIntegrator(this->ctx, "Mechanics", "Cube",
                                         library, behaviour) |
        or_die;
    auto& m = this->problem.getMaterial(this->ctx, "Cube", 0) | or_die;
    for (auto* s : {&m.s0, &m.s1}) {
      mgis::behaviour::setExternalStateVariable(this->ctx, *s, "Temperature",
                                                293.15) |
          or_die;
    }
    for (const auto& [boundary, component] :
         {std::pair{"Xmin", 0}, std::pair{"Ymin", 1}, std::pair{"Zmin", 2}}) {
      this->problem.addUniformDirichletBoundaryCondition(
          this->ctx, {{"Boundary", boundary}, {"Component", component}}) |
          or_die;
    }
    if (imposed_stress) {
      // the pressure is the opposite of the stress
      this->problem.addBoundaryCondition(
          this->ctx, std::make_unique<UniformImposedPressureBoundaryCondition>(
                         this->problem.getFiniteElementDiscretization(), "Xmax",
                         [loading](const real t) { return -loading(t); })) |
          or_die;
    } else {
      this->problem.addUniformDirichletBoundaryCondition(
          this->ctx, {{"Boundary", "Xmax"},
                      {"Component", 0},
                      {"LoadingEvolution", loading}}) |
          or_die;
    }
    // direct solver: MUMPS in parallel if MFEM provides it, otherwise the
    // conjugate gradient of MFEM with the tightest tolerance (the tangent
    // operators of the behaviour are symmetric positive definite)
    if (parallel) {
#ifdef MFEM_MGIS_HAVE_MUMPS
      this->problem.setLinearSolver(this->ctx, "MUMPSSolver", {}) | or_die;
#else  /* MFEM_MGIS_HAVE_MUMPS */
      this->problem.setLinearSolver(this->ctx, "CGSolver",
                                    {{"VerbosityLevel", 0},
                                     {"RelativeTolerance", 1e-14},
                                     {"AbsoluteTolerance", 0.},
                                     {"MaximumNumberOfIterations", 10000}}) |
          or_die;
#endif /* MFEM_MGIS_HAVE_MUMPS */
    } else {
      this->problem.setLinearSolver(
          this->ctx, unit_tests::sequential_direct_solver,
          unit_tests::getSequentialDirectSolverParameters()) |
          or_die;
    }
    this->problem.setSolverParameters(
        this->ctx,
        {{"VerbosityLevel", 0},
         {"RelativeTolerance", static_cast<double>(relative_tolerance)},
         {"AbsoluteTolerance", 0.},
         {"MaximumNumberOfIterations", maximum_number_of_iterations}}) |
        or_die;
    this->setPredictionStrategy(strategy);
  }  // end of Cube
  //! \brief select the prediction strategy
  void setPredictionStrategy(const mfem_mgis::PredictionStrategy strategy) {
    this->problem.setPredictionPolicy({.strategy = strategy});
  }  // end of setPredictionStrategy
  //! \brief set the maximum number of iterations of the resolutions
  void setMaximumNumberOfIterations(const int n) {
    auto or_die = this->ctx.getFatalFailureHandler();
    this->problem.setSolverParameters(this->ctx,
                                      {{"MaximumNumberOfIterations", n}}) |
        or_die;
  }  // end of setMaximumNumberOfIterations
  //! \brief set the absolute tolerance of the non linear resolutions
  void setAbsoluteTolerance(const double tolerance) {
    auto or_die = this->ctx.getFatalFailureHandler();
    this->problem.setSolverParameters(this->ctx,
                                      {{"AbsoluteTolerance", tolerance}}) |
        or_die;
  }  // end of setAbsoluteTolerance
  /*!
   * \return the result of the resolution of a time step
   * \param[in] t: time at the beginning of the time step
   * \param[in] dt: time increment
   */
  Resolution solve(const mfem_mgis::real t, const mfem_mgis::real dt) {
    // the messages of the resolution are read, then printed
    auto log = std::ostringstream{};
    auto local = mfem_mgis::Context{};
    local.setLogStream(log);
    const auto o = this->problem.solve(local, t, dt);
    const auto messages = log.str();
    std::cout << messages;
    auto r = Resolution{.status = o.status,
                        .iterations = o.iterations,
                        .reference = o.initial_residual_norm,
                        .norm = o.final_residual_norm,
                        .u = this->getUnknowns()};
    for (auto pos = messages.find(fallback_warning); pos != std::string::npos;
         pos = messages.find(fallback_warning, pos + 1)) {
      ++(r.fallbacks);
    }
    return r;
  }  // end of solve
  /*!
   * \return the result of a resolution limited to zero iteration: the
   * unknowns are the prediction, and the norm of the residual is the one of
   * the prediction
   * \param[in] t: time at the beginning of the time step
   * \param[in] dt: time increment
   */
  Resolution predict(const mfem_mgis::real t, const mfem_mgis::real dt) {
    auto or_die = this->ctx.getFatalFailureHandler();
    this->problem.setSolverParameters(this->ctx,
                                      {{"MaximumNumberOfIterations", 0}}) |
        or_die;
    const auto r = this->solve(t, dt);
    this->problem.setSolverParameters(
        this->ctx,
        {{"MaximumNumberOfIterations", maximum_number_of_iterations}}) |
        or_die;
    return r;
  }  // end of predict
  //! \return the unknowns at the end of the time step
  std::vector<mfem_mgis::real> getUnknowns() {
    const auto& u = this->problem.getUnknowns(mfem_mgis::ets);
    return {u.HostRead(), u.HostRead() + u.Size()};
  }  // end of getUnknowns
  //! \brief accept the time step
  void update() {
    auto or_die = this->ctx.getFatalFailureHandler();
    this->problem.update(this->ctx) | or_die;
  }  // end of update
  //! \brief cancel the time step
  void revert() {
    auto or_die = this->ctx.getFatalFailureHandler();
    this->problem.revert(this->ctx) | or_die;
  }  // end of revert
  //! \brief execution context
  mfem_mgis::Context ctx;
  //! \brief non linear evolution problem
  mfem_mgis::NonLinearEvolutionProblem problem;
};

/*!
 * \return the maximum difference between two arrays, relative to the
 * maximum absolute value of the first one
 */
static double getRelativeDifference(const std::vector<mfem_mgis::real>& a,
                                    const std::vector<mfem_mgis::real>& b) {
  if (a.size() != b.size()) {
    return 1;
  }
  auto d = 0.;
  auto m = 0.;
  for (std::size_t i = 0; i != a.size(); ++i) {
    d = std::max(d, static_cast<double>(std::abs(a[i] - b[i])));
    m = std::max(m, static_cast<double>(std::abs(a[i])));
  }
  return m > 0 ? d / m : d;
}  // end of getRelativeDifference

/*!
 * \brief displacement of the face `Xmax` in the elastic domain, per unit
 * of time, or of square of the time
 */
static constexpr auto elastic_displacement = mfem_mgis::real{1e-4};

/*!
 * \brief checks of the prediction on a linear problem, with a displacement
 * which is not linear in time
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] solver: non linear solver
 */
static void checkPredictions(const char* const mesh,
                             const char* const library,
                             const Solver& solver) {
  using namespace mfem_mgis;
  using Vector = std::vector<real>;
  constexpr auto extrapolation =
      PredictionStrategy::LINEAR_EXTRAPOLATION_PREDICTION;
  // differences of the predictions, of the solutions, of the norms and of
  // the ratios of norms (in float, the linear solver of the sequential
  // tests is iterative)
  constexpr auto values_tolerance = unit_tests::tolerance(1e-13, 1e-5);
  constexpr auto solutions_tolerance = unit_tests::tolerance(1e-8, 1e-3);
  constexpr auto norms_tolerance = unit_tests::tolerance(1e-9, 1e-3);
  constexpr auto ratios_tolerance = unit_tests::tolerance(1e-7, 5e-3);
  const auto loading = Loading{
      [](const real t) noexcept { return elastic_displacement * t * t; }};
  // the same computation with the linear extrapolation and with the default
  // prediction
  auto e = Cube(mesh, library, solver, loading, extrapolation);
  auto d = Cube(mesh, library, solver, loading,
                PredictionStrategy::DEFAULT_PREDICTION);
  const auto dofs = e.problem.getEssentialDegreesOfFreedom();
  // unknowns at t = 1: the imposed displacements are not null
  auto moved = Vector{};
  /*!
   * checks of a prediction of the time step from t to t + dt, the previous
   * time step starting at tp: the unknowns are u0 + ratio * (u0 - up) (u0
   * and up: unknowns at t and tp), except the imposed ones, and the ratio
   * of the norm of its residual to the reference is the one of the
   * differences between the imposed displacement and the displacement
   * extrapolated or not
   */
  auto check_prediction = [&](const Resolution& p, const Vector& u0,
                              const Vector& up, const real tp, const real t,
                              const real dt, const real ratio,
                              const std::string& step) {
    const auto name = solver.name + ", " + step;
    check(!p.status, name + ": convergence without iteration");
    check(p.iterations == 0, name + ": iterations made");
    if ((p.u.size() != u0.size()) || (up.size() != u0.size())) {
      check(false, name + ": invalid number of unknowns");
      return;
    }
    auto expected = Vector(u0.size());
    for (std::size_t i = 0; i != u0.size(); ++i) {
      expected[i] = u0[i] + ratio * (u0[i] - up[i]);
    }
    auto exact = true;
    for (const auto dof : dofs) {
      const auto i = static_cast<std::size_t>(dof);
      const auto v = areEqual(moved[i], real{0}) ? real{0} : loading(t + dt);
      exact = exact && areEqual(p.u[i], v);
      expected[i] = v;
    }
    check(exact, name + ": imposed values not exact");
    const auto difference = getRelativeDifference(expected, p.u);
    check(difference <= values_tolerance,
          name + ": unexpected prediction (relative difference " +
              std::to_string(difference) + ")");
    const auto next = static_cast<double>(loading(t + dt));
    const auto current = static_cast<double>(loading(t));
    const auto previous = static_cast<double>(loading(tp));
    const auto expected_ratio =
        std::abs(next - (current + ratio * (current - previous))) /
        std::abs(next - current);
    check((p.reference > 0) &&
              areClose(p.norm / p.reference, expected_ratio, ratios_tolerance),
          name +
              ": unexpected ratio of the residual of the prediction to "
              "the reference (" +
              std::to_string(p.norm / p.reference) + ", expected " +
              std::to_string(expected_ratio) + ")");
  };
  /*!
   * the resolutions of the same time step with the linear extrapolation and
   * with the default prediction have the same reference
   */
  auto check_reference = [&](const Resolution& a, const Resolution& b,
                             const std::string& step) {
    check(areClose(a.reference, b.reference, norms_tolerance),
          solver.name + ", " + step + ": the reference is not the one of " +
              "the default prediction (" + std::to_string(a.reference) + ", " +
              std::to_string(b.reference) + ")");
  };
  const auto zeros = Vector(e.getUnknowns().size(), real{0});
  // first time step: default prediction
  const auto e1 = e.solve(0, 1);
  const auto d1 = d.solve(0, 1);
  if (!(e1.status && d1.status)) {
    check(false, solver.name + ", first time step: resolution failed");
    return;
  }
  check_reference(e1, d1, "first time step");
  check(e1.iterations == d1.iterations,
        solver.name +
            ", first time step: iterations of the default "
            "prediction expected");
  check(getRelativeDifference(d1.u, e1.u) <= solutions_tolerance,
        solver.name + ", first time step: unexpected solution");
  moved = e1.u;
  // a resolution which starts from a solution converges without iteration
  for (auto* const c : {&e, &d}) {
    c->setAbsoluteTolerance(100 * static_cast<double>(relative_tolerance) *
                            e1.reference);
    c->update();
  }
  // second time step, same time increment
  {
    const auto p = e.predict(1, 1);
    check_prediction(p, e1.u, zeros, 0, 1, 1, 1, "second time step");
    check_reference(p, d.predict(1, 1), "second time step");
    e.revert();
    d.revert();
  }
  // after `revert`, half of the time increment
  {
    const auto p = e.predict(1, real{1} / 2);
    check_prediction(p, e1.u, zeros, 0, 1, real{1} / 2, real{1} / 2,
                     "second time step, halved");
    check_reference(p, d.predict(1, real{1} / 2), "second time step, halved");
    e.revert();
    d.revert();
  }
  const auto e2 = e.solve(1, real{1} / 2);
  const auto d2 = d.solve(1, real{1} / 2);
  if (!(e2.status && d2.status)) {
    check(false, solver.name + ", second time step: resolution failed");
    return;
  }
  check_reference(e2, d2, "second time step, resolution");
  check(e2.iterations >= 1, solver.name + ", second time step: no iteration");
  check(getRelativeDifference(d2.u, e2.u) <= solutions_tolerance,
        solver.name + ", second time step: unexpected solution");
  // second resolution of the same time step: no extrapolation, the solution
  // is kept
  {
    const auto s = e.solve(1, real{1} / 2);
    check(s.status && (s.iterations == 0),
          solver.name + ", second resolution of a time step: iterations made");
    check(getRelativeDifference(e2.u, s.u) <= 0,
          solver.name + ", second resolution of a time step: solution changed");
  }
  e.update();
  // third time step: the time increment is doubled
  {
    const auto p = e.predict(real{3} / 2, 1);
    check_prediction(p, e2.u, e1.u, 1, real{3} / 2, 1, 2, "third time step");
    e.revert();
  }
  // the default strategy releases the increment
  e.setPredictionStrategy(PredictionStrategy::DEFAULT_PREDICTION);
  {
    const auto p = e.predict(real{3} / 2, 1);
    check_prediction(p, e2.u, e1.u, 1, real{3} / 2, 1, 0,
                     "third time step, default prediction");
    e.revert();
  }
  // the linear extrapolation selected again has no increment to extrapolate
  e.setPredictionStrategy(extrapolation);
  {
    const auto p = e.predict(real{3} / 2, 1);
    check_prediction(p, e2.u, e1.u, 1, real{3} / 2, 1, 0,
                     "third time step, linear extrapolation selected again");
    e.revert();
  }
  const auto e3 = e.solve(real{3} / 2, 1);
  if (!e3.status) {
    check(false, solver.name + ", third time step: resolution failed");
    return;
  }
  e.update();
  // fourth time step: the increment of the third one is extrapolated
  {
    const auto p = e.predict(real{5} / 2, 1);
    check_prediction(p, e3.u, e2.u, real{3} / 2, real{5} / 2, 1, 1,
                     "fourth time step");
    e.revert();
  }
}  // end of checkPredictions

/*!
 * \brief with a displacement proportional to the time, on a linear problem,
 * the extrapolation is the solution
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] solver: non linear solver
 */
static void checkExactPrediction(const char* const mesh,
                                 const char* const library,
                                 const Solver& solver) {
  using namespace mfem_mgis;
  constexpr auto solutions_tolerance = unit_tests::tolerance(1e-8, 1e-3);
  constexpr auto norms_tolerance = unit_tests::tolerance(1e-9, 1e-3);
  const auto loading =
      Loading{[](const real t) noexcept { return elastic_displacement * t; }};
  auto e = Cube(mesh, library, solver, loading,
                PredictionStrategy::LINEAR_EXTRAPOLATION_PREDICTION);
  auto d = Cube(mesh, library, solver, loading,
                PredictionStrategy::DEFAULT_PREDICTION);
  // time increments: the same one, then halved, then doubled
  auto t = real{0};
  auto step = 0;
  for (const auto dt : {real{1}, real{1}, real{1} / 2, real{1}}) {
    const auto name = solver.name + ", proportional displacement, step " +
                      std::to_string(step);
    const auto re = e.solve(t, dt);
    const auto rd = d.solve(t, dt);
    if (!(re.status && rd.status)) {
      check(false, name + ": resolution failed");
      return;
    }
    check(areClose(re.reference, rd.reference, norms_tolerance),
          name + ": the reference is not the one of the default prediction");
    check(rd.iterations >= 1, name + ": no iteration without extrapolation");
    if (step == 0) {
      check(re.iterations == rd.iterations,
            name + ": iterations of the default prediction expected");
      // a resolution which starts from the solution converges without
      // iteration
      for (auto* const c : {&e, &d}) {
        c->setAbsoluteTolerance(100 * static_cast<double>(relative_tolerance) *
                                re.reference);
      }
    } else {
      check(re.iterations == 0, name + ": " + std::to_string(re.iterations) +
                                    " iterations with the extrapolation");
    }
    check(getRelativeDifference(rd.u, re.u) <= solutions_tolerance,
          name + ": unexpected solution");
    e.update();
    d.update();
    t += dt;
    ++step;
  }
}  // end of checkExactPrediction

/*!
 * \brief in the plastic domain, the solution is the one of the default
 * prediction, with at most as many iterations
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] solver: non linear solver
 */
static void checkPlasticCube(const char* const mesh,
                             const char* const library,
                             const Solver& solver) {
  using namespace mfem_mgis;
  constexpr auto solutions_tolerance = unit_tests::tolerance(1e-8, 1e-2);
  // strain of 1 % at the end: the cube yields during the third time step
  const auto loading =
      Loading{[](const real t) noexcept { return real{1e-2} * t; }};
  auto e = Cube(mesh, library, solver, loading,
                PredictionStrategy::LINEAR_EXTRAPOLATION_PREDICTION);
  auto d = Cube(mesh, library, solver, loading,
                PredictionStrategy::DEFAULT_PREDICTION);
  constexpr auto nsteps = 5;
  const auto dt = real{1} / nsteps;
  auto iterations = std::pair<size_type, size_type>{0, 0};
  for (int i = 0; i != nsteps; ++i) {
    const auto name = solver.name + ", plastic cube, step " + std::to_string(i);
    const auto t = static_cast<real>(i) * dt;
    const auto re = e.solve(t, dt);
    const auto rd = d.solve(t, dt);
    if (!(re.status && rd.status)) {
      check(false, name + ": resolution failed");
      return;
    }
    iterations.first += re.iterations;
    iterations.second += rd.iterations;
    const auto difference = getRelativeDifference(rd.u, re.u);
    check(difference <= solutions_tolerance,
          name + ": unexpected solution (relative difference " +
              std::to_string(difference) + ")");
    e.update();
    d.update();
  }
  std::cout << solver.name << ", plastic cube: " << iterations.first
            << " iterations with the linear extrapolation, "
            << iterations.second << " with the default prediction\n";
  check(iterations.first <= iterations.second,
        solver.name +
            ", plastic cube: more iterations with the linear "
            "extrapolation");
}  // end of checkPlasticCube

/*!
 * \brief the behaviour can't be integrated at the extrapolated estimate,
 * whereas it can be at the default prediction: the resolution starts again
 * from the default prediction
 *
 * The integration of `FailingIsotropicLinearElasticity` fails where the
 * absolute value of the trace of the strain exceeds 1e-3. The cube is in
 * uniaxial stress (Poisson ratio 0.3) and has two elements along the x
 * axis. The displacement of the face `Xmax` is a t until t = 1, then
 * a (1 - 3 (t - 1) / 8), with a = 4e-4. In the elements next to the face
 * `Xmax`, the trace of the strain at the first estimate of a time step is:
 * - from t = 0 to t = 1, 2 a = 8e-4 (default prediction);
 * - from t = 1 to t = 3 (reversal), -1.1 a = -4.4e-4 at the default
 *   prediction and -4.3 a = -1.72e-3 at the extrapolated estimate (three
 *   times the displacements at t = 1, except on the face `Xmax`): the
 *   integration fails at this one only;
 * - from t = 1 to t = 2, -0.35 a at the default prediction and -1.95 a =
 *   -7.8e-4 at the extrapolated estimate, which is thus kept, although the
 *   norm of its residual is 11 / 3 times the reference;
 * - from t = 3 to t = 4, the extrapolated estimate is the solution.
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] solver: non linear solver
 */
static void checkIntegrationFailure(const char* const mesh,
                                    const char* const library,
                                    const Solver& solver) {
  using namespace mfem_mgis;
  constexpr auto solutions_tolerance = unit_tests::tolerance(1e-8, 1e-3);
  constexpr auto norms_tolerance = unit_tests::tolerance(1e-9, 1e-3);
  constexpr auto ratios_tolerance = unit_tests::tolerance(1e-7, 5e-3);
  constexpr auto behaviour = "FailingIsotropicLinearElasticity";
  // the warnings are printed by the main process
  const auto main_process = getMPIRankForDefaultCommunicator() == 0;
  const auto loading = Loading{[](const real t) noexcept {
    constexpr auto a = real{4e-4};
    return t <= 1 ? a * t : a * (1 - 3 * (t - 1) / 8);
  }};
  auto e = Cube(mesh, library, solver, loading,
                PredictionStrategy::LINEAR_EXTRAPOLATION_PREDICTION, behaviour);
  auto d = Cube(mesh, library, solver, loading,
                PredictionStrategy::DEFAULT_PREDICTION, behaviour);
  /*!
   * the resolution, or the first estimate of the resolution, is the one of
   * the default strategy
   */
  auto check_default = [&](const Resolution& re, const Resolution& rd,
                           const std::size_t fallbacks,
                           const std::string& step) {
    const auto name = solver.name + ", integration failure, " + step;
    check(re.status == rd.status, name + ": unexpected status");
    check(re.iterations == rd.iterations,
          name + ": " + std::to_string(re.iterations) + " iterations, " +
              std::to_string(rd.iterations) + " with the default prediction");
    check(areClose(re.reference, rd.reference, norms_tolerance),
          name + ": the reference is not the one of the default prediction");
    check(
        areClose(re.norm, rd.norm, norms_tolerance) || (re.status && rd.status),
        name + ": unexpected norm of the residual");
    const auto difference = getRelativeDifference(rd.u, re.u);
    check(difference <= solutions_tolerance,
          name + ": unexpected unknowns (relative difference " +
              std::to_string(difference) + ")");
    check((!main_process) || (re.fallbacks == fallbacks),
          name + ": " + std::to_string(re.fallbacks) + " warnings, " +
              std::to_string(fallbacks) + " expected");
    check(rd.fallbacks == 0,
          name + ": warning printed with the default prediction");
  };
  // first time step: default prediction
  const auto e1 = e.solve(0, 1);
  const auto d1 = d.solve(0, 1);
  if (!(e1.status && d1.status)) {
    check(false, solver.name +
                     ", integration failure, first time step: "
                     "resolution failed");
    return;
  }
  check_default(e1, d1, 0, "first time step");
  // a resolution which starts from a solution converges without iteration
  for (auto* const c : {&e, &d}) {
    c->setAbsoluteTolerance(100 * static_cast<double>(relative_tolerance) *
                            e1.reference);
    c->update();
  }
  // from t = 1 to t = 3: the first estimate is the default prediction (the
  // unknowns at t = 1 with the values imposed at t = 3), bit for bit, and
  // the norm of its residual is the reference
  {
    const auto pe = e.predict(1, 2);
    const auto pd = d.predict(1, 2);
    check_default(pe, pd, 1, "reversal, first estimate");
    check((!pe.status) && (pe.iterations == 0),
          solver.name + ", integration failure, reversal: iterations made");
    auto expected = e1.u;
    for (const auto dof : e.problem.getEssentialDegreesOfFreedom()) {
      const auto i = static_cast<std::size_t>(dof);
      if (i < expected.size()) {
        expected[i] = areEqual(e1.u[i], real{0}) ? real{0} : loading(3);
      }
    }
    check(getRelativeDifference(expected, pe.u) <= 0,
          solver.name +
              ", integration failure, reversal: the first "
              "estimate is not the default prediction");
    check((pe.reference > 0) &&
              areClose(pe.norm / pe.reference, 1., ratios_tolerance),
          solver.name +
              ", integration failure, reversal: unexpected ratio "
              "of the residual of the first estimate to the reference (" +
              std::to_string(pe.norm / pe.reference) + ")");
    e.revert();
    d.revert();
  }
  // from t = 1 to t = 2: the behaviour is integrated at the extrapolated
  // estimate, which is kept
  {
    const auto pe = e.predict(1, 1);
    check((!pe.status) && (pe.iterations == 0),
          solver.name +
              ", integration failure, half of the reversal: "
              "iterations made");
    check((pe.reference > 0) &&
              areClose(pe.norm / pe.reference, 11. / 3, ratios_tolerance),
          solver.name +
              ", integration failure, half of the reversal: "
              "unexpected ratio of the residual of the prediction to the "
              "reference (" +
              std::to_string(pe.norm / pe.reference) + ", expected 11 / 3)");
    check(pe.fallbacks == 0, solver.name +
                                 ", integration failure, half of the reversal: "
                                 "the extrapolated estimate was discarded");
    e.revert();
  }
  // from t = 1 to t = 3 again, after `revert`: resolution of the default
  // strategy
  const auto e2 = e.solve(1, 2);
  const auto d2 = d.solve(1, 2);
  check(e2.status, solver.name +
                       ", integration failure, reversal: resolution "
                       "failed");
  check_default(e2, d2, 1, "reversal");
  if (!(e2.status && d2.status)) {
    return;
  }
  // second resolution of the same time step: the solution is kept, without
  // warning
  {
    const auto s = e.solve(1, 2);
    check(s.status && (s.iterations == 0) && (s.fallbacks == 0) &&
              (getRelativeDifference(e2.u, s.u) <= 0),
          solver.name +
              ", integration failure, second resolution of the "
              "reversal: solution changed");
  }
  e.update();
  d.update();
  // from t = 3 to t = 4: the increment of the time step solved from the
  // default prediction is extrapolated, and is the solution
  const auto e3 = e.solve(3, 1);
  const auto d3 = d.solve(3, 1);
  if (!(e3.status && d3.status)) {
    check(false, solver.name +
                     ", integration failure, third time step: "
                     "resolution failed");
    return;
  }
  check(areClose(e3.reference, d3.reference, norms_tolerance),
        solver.name +
            ", integration failure, third time step: the reference "
            "is not the one of the default prediction");
  check((e3.iterations == 0) && (d3.iterations >= 1),
        solver.name + ", integration failure, third time step: " +
            std::to_string(e3.iterations) +
            " iterations with the "
            "extrapolation, " +
            std::to_string(d3.iterations) + " without");
  check(e3.fallbacks == 0, solver.name +
                               ", integration failure, third time step: the "
                               "extrapolated estimate was discarded");
  check(getRelativeDifference(d3.u, e3.u) <= solutions_tolerance,
        solver.name +
            ", integration failure, third time step: unexpected "
            "solution");
}  // end of checkIntegrationFailure

/*!
 * \brief the behaviour is integrated at the extrapolated estimate, but not
 * at the following estimates of the resolution: these failures are handled
 * by the non linear solver, and the resolution does not start again from
 * the default prediction
 *
 * The cube of `FailingIsotropicLinearElasticity` (trace of the strain
 * bounded to 1e-3, see `checkIntegrationFailure`) is loaded by a uniaxial
 * stress, the trace of the strain of the solution being the stress divided
 * by 5.2e11: 3.8e-4 at t = 1, and 1.5e-3 at t = 2, where the behaviour
 * can't be integrated. The trace of the strain of the extrapolated estimate
 * of the second time step is 7.7e-4.
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] solver: non linear solver
 */
static void checkFailureAfterTheFirstIntegration(const char* const mesh,
                                                 const char* const library,
                                                 const Solver& solver) {
  using namespace mfem_mgis;
  const auto name = solver.name + ", failure after the first integration";
  const auto loading = Loading{[](const real t) noexcept {
    return t <= 1 ? real{2e8} * t : real{2e8} + real{5.8e8} * (t - 1);
  }};
  for (const auto strategy :
       {PredictionStrategy::LINEAR_EXTRAPOLATION_PREDICTION,
        PredictionStrategy::DEFAULT_PREDICTION}) {
    auto c = Cube(mesh, library, solver, loading, strategy,
                  "FailingIsotropicLinearElasticity", true);
    const auto r1 = c.solve(0, 1);
    if (!r1.status) {
      check(false, name + ", first time step: resolution failed");
      return;
    }
    c.update();
    c.setMaximumNumberOfIterations(20);
    const auto r2 = c.solve(1, 1);
    check(!r2.status, name + ": convergence beyond the limit of the behaviour");
    check((r1.fallbacks == 0) && (r2.fallbacks == 0),
          name + ": the resolution started again from the default prediction");
  }
}  // end of checkFailureAfterTheFirstIntegration

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
  // direct solver, used as linear solver of the elastic operator: the one
  // of the problem
  const auto direct_solver = []() -> Parameter {
    if (parallel) {
#ifdef MFEM_MGIS_HAVE_MUMPS
      return "MUMPSSolver";
#else  /* MFEM_MGIS_HAVE_MUMPS */
      return Parameters{
          {"Name", "CGSolver"},
          {"Options", Parameters{{"VerbosityLevel", 0},
                                 {"RelativeTolerance", 1e-14},
                                 {"AbsoluteTolerance", 0.},
                                 {"MaximumNumberOfIterations", 10000}}}};
#endif /* MFEM_MGIS_HAVE_MUMPS */
    }
    if (unit_tests::has_sequential_direct_solver) {
      return unit_tests::sequential_direct_solver;
    }
    return Parameters{
        {"Name", unit_tests::sequential_direct_solver},
        {"Options", unit_tests::getSequentialDirectSolverParameters()}};
  }();
  const auto solvers = std::vector<Solver>{
      {"Newton solver", {}, true},
      {"fixed point on the elastic operator",
       Parameters{{"FixedPoint", Parameters{}}}, false},
      {"NewtonK0",
       Parameters{
           {"FixedPoint",
            Parameters{{"Method", "NewtonK0"},
                       {"ElasticOperatorUpdate", "Once"},
                       {"ElasticOperatorLinearSolver", direct_solver},
                       {"TangentLinearSolver",
                        Parameters{{"CGSolver",
                                    Parameters{{"RelativeTolerance", 1e-12},
                                               {"MaximumNumberOfIterations",
                                                2000}}}}}}}},
       true}};
  for (const auto& solver : solvers) {
    checkPredictions(mesh, library, solver);
    checkExactPrediction(mesh, library, solver);
    checkIntegrationFailure(mesh, library, solver);
    checkFailureAfterTheFirstIntegration(mesh, library, solver);
    if (solver.newton) {
      checkPlasticCube(mesh, library, solver);
    }
  }
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
