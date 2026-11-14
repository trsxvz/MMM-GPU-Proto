/*!
 * \file   NewtonLineSearchTest.cxx
 * \brief  Test of the line search of the Newton solver (option
 * `LineSearch`).
 *
 * A cube made of the `SaturatingElasticity` behaviour is loaded in uniaxial
 * stress close to the saturation, then partly unloaded, or loaded in the
 * other direction, in one time step: the strain is uniform, and the Newton
 * iterations are the ones of a scalar equation. Without line search, the
 * second resolution does not converge: the iterates oscillate between two
 * saturated states. With the line search, it converges to the analytical
 * solution, the norms of the residual decrease, and the step lengths are
 * the expected ones. A full step being accepted, the iterations are the
 * ones of the solver without line search. The option may be changed
 * between two time steps. With the first choice of Eisenstat and Walker,
 * the tolerance which follows a reduced step uses the norm of the linear
 * model at this step.
 *
 * A cube made of the `FailingIsotropicLinearElasticity` behaviour, whose
 * integration fails above a trace of the strain, is loaded beyond it: the
 * trials at which the integration fails are refused, the numbers of
 * reductions are the expected ones, and the resolution fails when no step
 * is acceptable.
 *
 * On the plastic cube of `FixedPointSolverTest`, the line search gives the
 * solution of the solver without line search. A resolution whose tolerance
 * can't be reached stops at the floor of the rounding errors and is
 * accepted.
 * \date   14/11/2026
 */

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <utility>
#include <iostream>
#include <optional>
#include <algorithm>
#include <functional>
#include "mfem/general/device.hpp"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MPI.hxx"
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
//! \brief if true, the numerical type of mfem-mgis is double
static constexpr auto double_precision =
    sizeof(mfem_mgis::real) == sizeof(double);

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

/*!
 * \brief material parameters of the `SaturatingElasticity` behaviour
 * (default values of the behaviour): bulk modulus, shear modulus,
 * saturation stress and residual shear modulus
 */
static constexpr auto bulk_modulus = 150e9;
static constexpr auto shear_modulus = 80e9;
static constexpr auto saturation_stress = 240e6;
static constexpr auto residual_shear_modulus = 4e9;

/*!
 * \return the axial and the lateral strains of the `SaturatingElasticity`
 * behaviour under an uniaxial stress
 * \param[in] p: stress
 */
static std::array<double, 2> getUniaxialStrains(const double p) {
  // equivalent stress as a function of the equivalent strain, inverted by
  // bisection
  auto seq = [](const double e) {
    return saturation_stress *
               std::tanh(3 * shear_modulus * e / saturation_stress) +
           3 * residual_shear_modulus * e;
  };
  auto lower = 0.;
  auto upper = 1.;
  for (int i = 0; i != 200; ++i) {
    const auto middle = (lower + upper) / 2;
    if (seq(middle) < std::abs(p)) {
      lower = middle;
    } else {
      upper = middle;
    }
  }
  const auto e = std::copysign((lower + upper) / 2, p);
  return {e + p / (9 * bulk_modulus), -e / 2 + p / (9 * bulk_modulus)};
}  // end of getUniaxialStrains

//! \brief information on an iteration
struct Iteration {
  //! \brief norm of the residual at the end of the iteration, if any
  std::optional<double> norm;
  //! \brief step length of the line search, if any
  std::optional<double> step_length;
  //! \brief number of reductions of the step by the line search
  mfem_mgis::size_type reductions = 0;
  //! \brief norm of the linear model at the step retained, if computed
  std::optional<double> linear_model_residual_norm;
  //! \brief relative tolerance of the linear solver, if adapted
  std::optional<double> tolerance;
  //! \brief norm of the residual of the linear solve, if computed
  std::optional<double> linear_residual_norm;
};

//! \brief information on a resolution
struct Resolution {
  //! \brief status
  bool status = false;
  //! \brief initial norm of the residual
  double initial_norm = 0;
  //! \brief number of iterations
  mfem_mgis::size_type number_of_iterations = 0;
  //! \brief iterations
  std::vector<Iteration> iterations;
};

//! \brief results of a computation
struct Results {
  //! \brief true if all the resolutions converged
  bool status = false;
  //! \brief resolutions of the time steps, the last one being the failed one
  std::vector<Resolution> resolutions;
  //! \brief smallest and largest unknowns at the end of the computation
  std::array<double, 2> extrema = {0, 0};
  //! \brief unknowns at the end of the computation
  std::vector<mfem_mgis::real> u;
};

//! \brief description of a computation
struct Computation {
  //! \brief behaviour
  const char* behaviour = "SaturatingElasticity";
  /*!
   * \brief uniaxial stresses imposed at the end of each time step, by a
   * pressure on the face `Xmax`. If empty, the displacement `1e-2 t` is
   * imposed on this face, during two time steps.
   */
  std::vector<double> stresses = {};
  //! \brief value of the `NonLinearSolver` parameter, the default if empty
  mfem_mgis::Parameters solver = {};
  //! \brief parameters of the non linear solver replacing the default ones
  mfem_mgis::Parameters solver_parameters = {};
  //! \brief options of the non linear solver changed before the second step
  mfem_mgis::Parameters changes = {};
  //! \brief linear solver, a direct solver if null
  const char* linear_solver = nullptr;
  //! \brief options of the linear solver
  mfem_mgis::Parameters linear_solver_options = {};
};

//! \return the line search and the other options of the Newton solver
static mfem_mgis::Parameters makeNewtonSolver(
    const mfem_mgis::Parameters& options) {
  return mfem_mgis::Parameters{{"Newton", options}};
}  // end of makeNewtonSolver

/*!
 * \return the results of a computation
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] c: computation
 */
static Results compute(const char* const mesh,
                       const char* const library,
                       const Computation& c) {
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
  if (!c.solver.empty()) {
    params.insert(throwing, {{"NonLinearSolver", c.solver}});
  }
  auto problem = construct<NonLinearEvolutionProblem>(ctx, params) | or_die;
  problem.addBehaviourIntegrator(ctx, "Mechanics", "Cube", library,
                                 c.behaviour) |
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
  if (c.stresses.empty()) {
    problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", "Xmax"},
              {"Component", 0},
              {"LoadingEvolution",
               [](const real t) noexcept { return 1e-2 * t; }}}) |
        or_die;
  } else {
    // the stress is linear over each time step, of unit duration; the
    // pressure is its opposite
    problem.addBoundaryCondition(
        ctx,
        std::make_unique<UniformImposedPressureBoundaryCondition>(
            problem.getFiniteElementDiscretization(), "Xmax",
            [stresses = c.stresses](const real t) {
              const auto n = static_cast<int>(stresses.size());
              const auto i =
                  std::clamp(static_cast<int>(std::ceil(t)) - 1, 0, n - 1);
              const auto s0 = (i == 0) ? 0. : stresses[i - 1];
              const auto s1 = stresses[i];
              return static_cast<real>(-(s0 + (s1 - s0) * (t - i)));
            })) |
        or_die;
  }
  // direct solver: MUMPS in parallel if MFEM provides it, otherwise the
  // conjugate gradient of MFEM with the tightest tolerance (the tangent
  // operators of the behaviours of the test are symmetric positive
  // definite)
#ifdef MFEM_MGIS_HAVE_MUMPS
  const auto* const parallel_solver = "MUMPSSolver";
  const auto parallel_solver_options = Parameters{};
#else  /* MFEM_MGIS_HAVE_MUMPS */
  const auto* const parallel_solver = "CGSolver";
  const auto parallel_solver_options =
      Parameters{{"VerbosityLevel", 0},
                 {"RelativeTolerance", 1e-14},
                 {"AbsoluteTolerance", 0.},
                 {"MaximumNumberOfIterations", 10000}};
#endif /* MFEM_MGIS_HAVE_MUMPS */
  const auto* const ls =
      (c.linear_solver != nullptr)
          ? c.linear_solver
          : (parallel ? parallel_solver : unit_tests::sequential_direct_solver);
  const auto ls_options =
      (c.linear_solver != nullptr)
          ? c.linear_solver_options
          : (parallel ? parallel_solver_options
                      : unit_tests::getSequentialDirectSolverParameters());
  problem.setLinearSolver(ctx, ls, ls_options) | or_die;
  auto solver_parameters = Parameters{{"VerbosityLevel", 0},
                                      {"RelativeTolerance", 1e-12},
                                      {"AbsoluteTolerance", 0.},
                                      {"MaximumNumberOfIterations", 50}};
  for (const auto& [k, v] : c.solver_parameters) {
    solver_parameters.replaceOrInsert(k, v);
  }
  problem.setSolverParameters(ctx, solver_parameters) | or_die;
  auto r = Results{};
  const auto nsteps =
      c.stresses.empty() ? 2 : static_cast<int>(c.stresses.size());
  for (int i = 0; i != nsteps; ++i) {
    if ((i == 1) && (!c.changes.empty())) {
      problem.setSolverParameters(ctx, c.changes) | or_die;
    }
    auto local = Context{};
    const auto o = problem.solve(local, real(i), 1);
    auto& resolution = r.resolutions.emplace_back();
    resolution.status = o.status;
    resolution.initial_norm = o.initial_residual_norm;
    resolution.number_of_iterations = o.iterations;
    for (const auto& information : o.iterations_information) {
      const auto& p = get<Parameters>(throwing, information);
      auto& iteration = resolution.iterations.emplace_back();
      if (contains(p, "Norm")) {
        iteration.norm = get<double>(throwing, p, "Norm");
      }
      if (contains(p, "LineSearch")) {
        const auto& l = get<Parameters>(throwing, p, "LineSearch");
        iteration.step_length = get<double>(throwing, l, "StepLength");
        iteration.reductions = get<int>(throwing, l, "NumberOfStepReductions");
        if (contains(l, "LinearModelResidualNorm")) {
          iteration.linear_model_residual_norm =
              get<double>(throwing, l, "LinearModelResidualNorm");
        }
      }
      if (contains(p, "LinearSolver")) {
        const auto& l = get<Parameters>(throwing, p, "LinearSolver");
        if (contains(l, "RelativeTolerance")) {
          iteration.tolerance = get<double>(throwing, l, "RelativeTolerance");
        }
        if (contains(l, "ResidualNorm")) {
          iteration.linear_residual_norm =
              get<double>(throwing, l, "ResidualNorm");
        }
      }
    }
    if (!o.status) {
      return r;
    }
    problem.update(ctx) | or_die;
  }
  const auto& u = problem.getUnknowns(ets);
  r.u.assign(u.HostRead(), u.HostRead() + u.Size());
  // the opposite of the smallest unknown and the largest one, on all the
  // processes
  auto extrema = std::array<double, 2>{0, 0};
  for (const auto v : r.u) {
    extrema[0] = std::max(extrema[0], -static_cast<double>(v));
    extrema[1] = std::max(extrema[1], static_cast<double>(v));
  }
#ifdef MFEM_USE_MPI
  if (parallel) {
    MPI_Allreduce(MPI_IN_PLACE, extrema.data(), 2, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
  }
#endif /* MFEM_USE_MPI */
  r.extrema = {-extrema[0], extrema[1]};
  r.status = true;
  return r;
}  // end of compute

/*!
 * \brief checks that the norms of the residual of the resolutions decrease
 * \param[in] r: results
 * \param[in] name: name of the computation
 */
static void checkDecreasingNorms(const Results& r, const std::string& name) {
  for (const auto& resolution : r.resolutions) {
    auto previous = resolution.initial_norm;
    for (const auto& iteration : resolution.iterations) {
      if (!iteration.norm.has_value()) {
        continue;
      }
      check(*(iteration.norm) < previous,
            name + ": the norm of the residual does not decrease");
      previous = *(iteration.norm);
    }
  }
}  // end of checkDecreasingNorms

/*!
 * \brief checks the solution of a cube of `SaturatingElasticity` under an
 * uniaxial stress: the displacements are linear, the largest and smallest
 * ones being the strains (the cube has a unit side)
 * \param[in] r: results
 * \param[in] p: uniaxial stress
 * \param[in] name: name of the computation
 */
static void checkUniaxialSolution(const Results& r,
                                  const double p,
                                  const std::string& name) {
  constexpr auto tolerance = mfem_mgis::unit_tests::tolerance(1e-8, 1e-3);
  const auto e = getUniaxialStrains(p);
  const auto expected =
      std::array<double, 2>{std::min(e[0], e[1]), std::max(e[0], e[1])};
  check(r.status && areClose(r.extrema[0], expected[0], tolerance) &&
            areClose(r.extrema[1], expected[1], tolerance),
        name + ": solution (strains " + std::to_string(r.extrema[0]) + " and " +
            std::to_string(r.extrema[1]) + ", " + std::to_string(expected[0]) +
            " and " + std::to_string(expected[1]) + " expected)");
}  // end of checkUniaxialSolution

/*!
 * \brief cube of `SaturatingElasticity` loaded close to the saturation,
 * then unloaded
 * \param[in] mesh: mesh
 * \param[in] library: behaviour library
 * \param[in] stresses: stresses at the end of the two time steps, divided
 * by the saturation stress
 * \param[in] step_lengths: step lengths expected for the first iterations
 * of the second time step, the following ones being full steps
 */
static void checkSaturation(const char* const mesh,
                            const char* const library,
                            const std::array<double, 2>& stresses,
                            const std::vector<double>& step_lengths) {
  using namespace mfem_mgis;
  const auto loading = std::vector<double>{stresses[0] * saturation_stress,
                                           stresses[1] * saturation_stress};
  const auto name = "saturation (" + std::to_string(stresses[0]) + ", " +
                    std::to_string(stresses[1]) + ")";
  const auto backtracking =
      makeNewtonSolver({{"LineSearch", "Backtracking"}});
  // without line search: the first resolution converges, the second one
  // does not (the case would not test anything otherwise)
  const auto none = compute(mesh, library, {.stresses = loading});
  check((!none.status) && (none.resolutions.size() == 2) &&
            (none.resolutions[0].status),
        name + ": convergence without line search");
  // with the line search
  const auto r =
      compute(mesh, library, {.stresses = loading, .solver = backtracking});
  check(r.status, name + ": no convergence with the line search");
  if ((!r.status) || (none.resolutions.size() != 2)) {
    return;
  }
  checkUniaxialSolution(r, loading[1], name);
  checkDecreasingNorms(r, name);
  // first resolution: the full steps are accepted, and the iterations are
  // the ones of the solver without line search (the last norms, at the
  // level of the rounding errors, are not compared). In float, the last
  // iterations, at the floor of the rounding errors, may be reduced.
  const auto& first = r.resolutions[0];
  const auto& reference = none.resolutions[0];
  // relative norm above which the norms of the two computations are
  // compared, and tolerance of the comparison
  constexpr auto norms_floor = mfem_mgis::unit_tests::tolerance(1e-6, 1e-3);
  if (double_precision) {
    check(first.number_of_iterations == reference.number_of_iterations,
          name + ": number of iterations of the first resolution (" +
              std::to_string(first.number_of_iterations) + ", " +
              std::to_string(reference.number_of_iterations) +
              " without line search)");
  }
  for (std::size_t i = 0;
       i != std::min(first.iterations.size(), reference.iterations.size());
       ++i) {
    const auto& iteration = first.iterations[i];
    check(iteration.step_length.has_value(),
          name + ": no information on the line search");
    if (double_precision) {
      check((iteration.step_length.value_or(0) == 1) &&
                (iteration.reductions == 0),
            name + ": step reduced during the first resolution");
    }
    const auto n = iteration.norm.value_or(0);
    if (n > norms_floor * first.initial_norm) {
      check(areClose(n, reference.iterations[i].norm.value_or(0),
                     norms_floor),
            name + ": norms of the first resolution");
    }
  }
  // second resolution: expected step lengths
  const auto& second = r.resolutions[1];
  check(second.iterations.size() >= step_lengths.size(),
        name + ": number of iterations of the second resolution");
  auto description = std::string{};
  for (std::size_t i = 0; i != second.iterations.size(); ++i) {
    const auto& iteration = second.iterations[i];
    const auto expected = (i < step_lengths.size()) ? step_lengths[i] : 1.;
    const auto theta = iteration.step_length.value_or(0);
    description += " " + std::to_string(theta);
    if ((i < step_lengths.size()) || double_precision) {
      check((theta == expected) &&
                (theta == std::ldexp(1., -iteration.reductions)),
            name + ": step length of the iteration " + std::to_string(i + 1) +
                " of the second resolution (" + std::to_string(theta) + ", " +
                std::to_string(expected) + " expected)");
    }
  }
  std::cout << name << ": " << second.number_of_iterations
            << " iterations of the second resolution, step lengths"
            << description << '\n';
  // line search enabled before the second time step only
  const auto changed = compute(
      mesh, library,
      {.stresses = loading, .changes = {{"LineSearch", "Backtracking"}}});
  check(changed.status &&
            (!changed.resolutions[0].iterations[0].step_length.has_value()) &&
            (changed.resolutions[1].iterations[0].step_length.has_value()),
        name + ": line search enabled before the second time step");
  // same state at the end of the first time step, same iterations (in
  // float, the last iterations of the first time step may differ)
  if (changed.status && double_precision) {
    check(changed.resolutions[1].number_of_iterations ==
              second.number_of_iterations,
          name + ": iterations after the line search was enabled");
  }
  if (changed.status) {
    checkUniaxialSolution(changed, loading[1], name + ", option changed");
  }
  // line search disabled before the second time step
  const auto disabled = compute(
      mesh, library,
      {.stresses = loading,
       .solver = backtracking,
       .changes = {{"LineSearch", "None"}}});
  check((!disabled.status) && (disabled.resolutions.size() == 2),
        name + ": line search disabled before the second time step");
}  // end of checkSaturation

/*!
 * \brief the tolerance of the linear solver which follows a reduced step
 * (first choice of Eisenstat and Walker) uses the norm of the linear model
 * at this step. The first linear solve of each resolution being precise,
 * the first step of the second resolution is halved as with a direct
 * solver, and the norm of the linear model is then half the norm of the
 * residual.
 * \param[in] mesh: mesh
 * \param[in] library: behaviour library
 */
static void checkAdaptiveTolerance(const char* const mesh,
                                   const char* const library) {
  using namespace mfem_mgis;
  const auto name = std::string{"saturation, EisenstatWalker1"};
  const auto loading =
      std::vector<double>{0.98 * saturation_stress, 0.4 * saturation_stress};
  constexpr auto initial_tolerance = 1e-10;
  // absolute convergence criterion, about 1e-11 times the initial residuals
  // (2e-4 in float), so that the tolerances can be computed again
  constexpr auto goal = mfem_mgis::unit_tests::tolerance(1e-3, 1e4);
  const auto r = compute(
      mesh, library,
      {.stresses = loading,
       .solver = makeNewtonSolver(
           {{"LineSearch", "Backtracking"},
            {"AdaptiveLinearSolverTolerance",
             Parameters{{"EisenstatWalker1",
                         Parameters{{"InitialTolerance",
                                     initial_tolerance}}}}}}),
       .solver_parameters = {{"RelativeTolerance", 0.},
                             {"AbsoluteTolerance", goal}},
       .linear_solver = "CGSolver",
       .linear_solver_options = {{"RelativeTolerance", 1e-12},
                                 {"AbsoluteTolerance", 0.},
                                 {"MaximumNumberOfIterations", 10000}}});
  check(r.status, name + ": no convergence");
  if (!r.status) {
    return;
  }
  constexpr auto tolerance = mfem_mgis::unit_tests::tolerance(1e-6, 1e-3);
  const auto e = getUniaxialStrains(loading[1]);
  check(areClose(r.extrema[0], std::min(e[0], e[1]), tolerance) &&
            areClose(r.extrema[1], std::max(e[0], e[1]), tolerance),
        name + ": solution");
  checkDecreasingNorms(r, name);
  const auto golden_ratio = (1 + std::sqrt(5.)) / 2;
  auto reduced = 0;
  auto checked = 0;
  // rounding errors of the norm of the linear model
  constexpr auto model_tolerance = mfem_mgis::unit_tests::tolerance(1e-6, 1e-4);
  for (const auto& resolution : r.resolutions) {
    const auto& iterations = resolution.iterations;
    for (std::size_t k = 0; k != iterations.size(); ++k) {
      const auto& iteration = iterations[k];
      // norm of the residual before this iteration
      const auto f = (k == 0) ? resolution.initial_norm
                              : iterations[k - 1].norm.value_or(0);
      const auto theta = iteration.step_length.value_or(0);
      const auto lr = iteration.linear_residual_norm.value_or(-1);
      check(lr >= 0, name + ": no norm of the residual of the linear solve");
      // norm of the linear model at the step retained
      auto model = lr;
      if (theta != 1) {
        ++reduced;
        check(iteration.linear_model_residual_norm.has_value(),
              name + ": no norm of the linear model after a reduced step");
        model = iteration.linear_model_residual_norm.value_or(-1);
        // norm of (1 - θ) r - θ lr
        check((model >= std::abs((1 - theta) * f - theta * lr) *
                               (1 - model_tolerance)) &&
                  (model <=
                   ((1 - theta) * f + theta * lr) * (1 + model_tolerance)),
              name + ": norm of the linear model (" + std::to_string(model) +
                  ")");
      } else {
        check(!iteration.linear_model_residual_norm.has_value(),
              name + ": norm of the linear model after a full step");
      }
      // tolerance of the next iteration
      if ((k + 1 == iterations.size()) ||
          (!iteration.tolerance.has_value()) ||
          (!iterations[k + 1].tolerance.has_value()) ||
          (!iteration.norm.has_value())) {
        continue;
      }
      const auto fn = *(iteration.norm);
      auto eta = std::abs(fn - model) / f;
      const auto safeguard = std::pow(*(iteration.tolerance), golden_ratio);
      if (safeguard > 0.1) {
        eta = std::max(eta, safeguard);
      }
      eta = std::min(std::max(eta, static_cast<double>(goal) / (2 * fn)), 0.9);
      const auto v = *(iterations[k + 1].tolerance);
      check(areClose(v, eta, mfem_mgis::unit_tests::tolerance(1e-10, 1e-4)),
            name + ": relative tolerance of a linear solve (" +
                std::to_string(v) + " instead of " + std::to_string(eta) +
                ")");
      ++checked;
    }
  }
  check(checked > 0, name + ": no relative tolerance checked");
  // first iteration of the second resolution
  const auto& second = r.resolutions[1];
  const auto& first_iteration = second.iterations[0];
  check((first_iteration.step_length.value_or(0) == 0.5) &&
            areClose(first_iteration.linear_model_residual_norm.value_or(0),
                     second.initial_norm / 2, tolerance),
        name + ": first iteration of the second resolution");
  std::cout << name << ": " << reduced << " reduced steps, " << checked
            << " tolerances checked\n";
}  // end of checkAdaptiveTolerance

/*!
 * \brief cube of `FailingIsotropicLinearElasticity` in uniaxial stress,
 * the trace of the strain of the solution being 1.5 times the one above
 * which the integration fails: the trials at which the integration fails
 * are refused, the estimates get closer to the limit, and the search
 * finally fails
 * \param[in] mesh: mesh
 * \param[in] library: behaviour library
 */
static void checkIntegrationFailures(const char* const mesh,
                                     const char* const library) {
  using namespace mfem_mgis;
  const auto name = std::string{"integration failures"};
  // maximum trace of the strain of the behaviour, and trace of the strain
  // of the solution, σ / (3 λ + 2 μ)
  constexpr auto maximum_trace = 1e-3;
  constexpr auto stress = 7.8e8;
  constexpr auto solution_trace = stress / (3 * 1.2e11 + 2 * 8e10);
  const auto r = compute(mesh, library,
                         {.behaviour = "FailingIsotropicLinearElasticity",
                          .stresses = {stress},
                          .solver = makeNewtonSolver(
                              {{"LineSearch", "Backtracking"}})});
  check((!r.status) && (r.resolutions.size() == 1),
        name + ": convergence beyond the limit of the behaviour");
  if (r.resolutions.empty()) {
    return;
  }
  checkDecreasingNorms(r, name);
  // direct evaluation: the problem is linear, the trace of the strain of
  // the trial of step length θ is t + θ (solution_trace - t)
  auto expected = std::vector<size_type>{};
  auto trace = 0.;
  auto failed = false;
  while (!failed) {
    auto theta = 1.;
    auto reductions = size_type{};
    while (trace + theta * (solution_trace - trace) > maximum_trace) {
      if (reductions == 10) {
        failed = true;
        break;
      }
      theta /= 2;
      ++reductions;
    }
    expected.push_back(reductions);
    trace += theta * (solution_trace - trace);
  }
  const auto& iterations = r.resolutions[0].iterations;
  auto description = std::string{};
  for (const auto& iteration : iterations) {
    description += " " + std::to_string(iteration.reductions);
  }
  std::cout << name << ": reductions" << description << '\n';
  // in float, the last reductions depend on the rounding errors: the first
  // three ones are compared
  const auto n = double_precision ? expected.size() : std::size_t{3};
  check(iterations.size() >= n, name + ": number of iterations");
  for (std::size_t i = 0; i != std::min(n, iterations.size()); ++i) {
    check(iterations[i].reductions == expected[i],
          name + ": " + std::to_string(iterations[i].reductions) +
              " reductions at the iteration " + std::to_string(i + 1) + ", " +
              std::to_string(expected[i]) + " expected");
  }
  // the search which fails made the maximum number of reductions
  if (!iterations.empty()) {
    const auto& last = iterations.back();
    check((!last.norm.has_value()) && (last.reductions == 10) &&
              (last.step_length.value_or(0) == std::ldexp(1., -10)),
          name + ": failure of the last search");
  }
  // maximum number of reductions given by the option: the first iteration
  // requires one reduction, the second one two reductions
  const auto bounded = compute(
      mesh, library,
      {.behaviour = "FailingIsotropicLinearElasticity",
       .stresses = {stress},
       .solver = makeNewtonSolver(
           {{"LineSearch",
             Parameters{
                 {"Backtracking",
                  Parameters{{"MaximumNumberOfStepReductions", 1}}}}}})});
  check((!bounded.status) && (bounded.resolutions.size() == 1) &&
            (bounded.resolutions[0].iterations.size() == 2) &&
            (bounded.resolutions[0].iterations[1].reductions == 1) &&
            (!bounded.resolutions[0].iterations[1].norm.has_value()),
        name + ": maximum number of reductions");
}  // end of checkIntegrationFailures

/*!
 * \brief plastic cube of `FixedPointSolverTest`: the line search gives the
 * solution of the solver without line search, with the same iterations if
 * no step is reduced
 * \param[in] mesh: mesh
 * \param[in] library: behaviour library
 */
static void checkPlasticCube(const char* const mesh,
                             const char* const library) {
  using namespace mfem_mgis;
  const auto name = std::string{"plastic cube"};
  const auto none = compute(
      mesh, library, {.behaviour = "IsotropicLinearHardeningPlasticity"});
  const auto r =
      compute(mesh, library,
              {.behaviour = "IsotropicLinearHardeningPlasticity",
               .solver = makeNewtonSolver({{"LineSearch", "Backtracking"}})});
  check(none.status && r.status, name + ": resolution failed");
  if (!(none.status && r.status)) {
    return;
  }
  auto reduced = 0;
  for (const auto& resolution : r.resolutions) {
    for (const auto& iteration : resolution.iterations) {
      if (iteration.reductions != 0) {
        ++reduced;
      }
    }
  }
  auto difference = 0.;
  auto maximum = 0.;
  if (none.u.size() == r.u.size()) {
    for (std::size_t i = 0; i != r.u.size(); ++i) {
      difference = std::max(
          difference, static_cast<double>(std::abs(r.u[i] - none.u[i])));
      maximum = std::max(maximum, static_cast<double>(std::abs(none.u[i])));
    }
  } else {
    difference = 1;
    maximum = 1;
  }
  std::cout << name << ": " << reduced
            << " reduced steps, relative difference to the solver without "
               "line search "
            << difference / maximum << '\n';
  check(difference <=
            mfem_mgis::unit_tests::tolerance(1e-8, 3e-4) * maximum,
        name + ": solution");
  if (reduced == 0) {
    for (std::size_t i = 0; i != r.resolutions.size(); ++i) {
      check(r.resolutions[i].number_of_iterations ==
                none.resolutions[i].number_of_iterations,
            name + ": number of iterations of the time step " +
                std::to_string(i + 1));
    }
  }
}  // end of checkPlasticCube

/*!
 * \brief a tolerance which can't be reached: the residual reaches the
 * floor of the rounding errors, where no step decreases it, and the
 * resolution is accepted (see the option `StagnationRelativeTolerance`)
 * \param[in] mesh: mesh
 * \param[in] library: behaviour library
 */
static void checkStagnation(const char* const mesh,
                            const char* const library) {
  using namespace mfem_mgis;
  const auto name = std::string{"stagnation"};
  const auto loading =
      std::vector<double>{0.98 * saturation_stress, 0.4 * saturation_stress};
  const auto r =
      compute(mesh, library,
              {.stresses = loading,
               .solver = makeNewtonSolver({{"LineSearch", "Backtracking"}}),
               .solver_parameters = {{"RelativeTolerance", 0.},
                                     {"AbsoluteTolerance", 0.}}});
  check(r.status, name + ": resolution refused at the floor of the residual");
  if (r.status) {
    checkUniaxialSolution(r, loading[1], name);
    checkDecreasingNorms(r, name);
  }
}  // end of checkStagnation

/*!
 * \brief an invalid option of the Newton solver is refused
 * \param[in] mesh: mesh
 * \param[in] options: options of the Newton solver
 * \param[in] name: name of the check
 * \param[in] message: part of the expected error message
 */
static void checkRefusal(const char* const mesh,
                         const mfem_mgis::Parameters& options,
                         const std::string& name,
                         const char* const message) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  const auto o = construct<NonLinearEvolutionProblem>(
      ctx, Parameters{{"MeshFileName", mesh},
                      {"FiniteElementFamily", "H1"},
                      {"FiniteElementOrder", 1},
                      {"UnknownsSize", 3},
                      {"Hypothesis", "Tridimensional"},
                      {"NonLinearSolver", makeNewtonSolver(options)},
                      {"Parallel", parallel}});
  check(isInvalid(o), name + " is accepted");
  // reading the error message clears it
  const auto e = ctx.getErrorMessage();
  check(e.find(message) != std::string::npos,
        name + ": unexpected error message (" + e + ")");
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
  // unloading from 0.98 to 0.4 times the saturation stress: the first step
  // is halved; loading in the other direction, from 0.95 to -0.5 times the
  // saturation stress: the second step is divided by four
  checkSaturation(mesh, library, {0.98, 0.4}, {0.5});
  checkSaturation(mesh, library, {0.95, -0.5}, {1, 0.25});
  checkAdaptiveTolerance(mesh, library);
  checkIntegrationFailures(mesh, library);
  checkPlasticCube(mesh, library);
  checkStagnation(mesh, library);
  checkRefusal(mesh, {{"LineSearch", "Unknown"}}, "an invalid line search",
               "invalid method");
  checkRefusal(
      mesh,
      {{"LineSearch",
        Parameters{{"Backtracking",
                    Parameters{{"SufficientDecrease", 2.}}}}}},
      "an invalid sufficient decrease", "sufficient decrease");
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
