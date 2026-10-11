/*!
 * \file   IntegrationRuleOrderTest.cxx
 * \brief  Test of the `IntegrationRuleOrder` parameter of the finite
 * element discretization: the behaviour integrators use the rule of the
 * given order (8 points on the hexahedra of order 1 for the order 3,
 * instead of 27 by default) for their quadrature spaces and their
 * integrations; on a mesh of affine hexahedra, whose elastic stiffness is
 * integrated exactly by both rules, the solutions are the same, with the
 * integration point by point and with the computations by kernels; a
 * negative order is refused.
 * \date   11/10/2026
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>
#include <vector>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

//! \brief results of a computation
struct Results {
  //! \brief number of integration points of each element of the material
  mfem_mgis::size_type points_per_element = 0;
  //! \brief number of points of the rule used by the behaviour integrator
  int rule_points = 0;
  //! \brief unknowns at the end of the computation
  std::vector<mfem_mgis::real> u;
};

/*!
 * \return the parameters of the problem
 * \param[in] mesh: mesh of the unit cube (hexahedra)
 * \param[in] order: order of the integration rules, if any
 * \param[in] kernels: if true, the behaviours are integrated by batches and
 * the gradients, the inner forces and the tangent operator are computed by
 * kernels
 */
static mfem_mgis::Parameters getParameters(const char* const mesh,
                                           const std::optional<int> order,
                                           const bool kernels) {
  using namespace mfem_mgis;
  const auto* const k = kernels ? "Kernels" : "ElementByElement";
  auto p = Parameters{
      {"MeshFileName", mesh},
      {"FiniteElementFamily", "H1"},
      {"FiniteElementOrder", 1},
      {"UnknownsSize", 3},
      {"NumberOfUniformRefinements", 1},
      {"Materials", Parameters{{"Cube", 1}}},
      {"Boundaries",
       Parameters{{"Ymin", 1}, {"Zmin", 2}, {"Xmax", 3}, {"Xmin", 5}}},
      {"Hypothesis", "Tridimensional"},
      {"BehaviourIntegration", kernels ? "Batch" : "PointByPoint"},
      {"GradientsComputation", k},
      {"InnerForcesComputation", k},
      {"TangentOperatorComputation", k},
      {"Parallel", false}};
  if (order.has_value()) {
    auto ctx = Context{};
    p.insert(ctx, "IntegrationRuleOrder", *order) |
        ctx.getFatalFailureHandler();
  }
  return p;
}  // end of getParameters

/*!
 * \brief elastic computation of a cube clamped on its face x = 0, the
 * face x = 1 being displaced along x and y (small displacements: the
 * plastic behaviour stays elastic)
 * \return the results
 * \param[in] mesh: mesh of the unit cube (hexahedra)
 * \param[in] library: behaviour library
 * \param[in] order: order of the integration rules, if any
 * \param[in] kernels: computations by kernels
 */
static Results compute(const char* const mesh,
                       const char* const library,
                       const std::optional<int> order,
                       const bool kernels) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  auto problem = construct<NonLinearEvolutionProblem>(
                     ctx, getParameters(mesh, order, kernels)) |
                 or_die;
  // the displacements stay in the elastic domain of the behaviour
  problem.addBehaviourIntegrator(ctx, "Mechanics", "Cube", library,
                                 "IsotropicLinearHardeningPlasticity") |
      or_die;
  auto& m = problem.getMaterial(ctx, "Cube", 0) | or_die;
  for (auto* s : {&m.s0, &m.s1}) {
    mgis::behaviour::setExternalStateVariable(ctx, *s, "Temperature", 293.15) |
        or_die;
  }
  for (const auto c : {0, 1, 2}) {
    problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", "Xmin"}, {"Component", c}}) |
        or_die;
  }
  problem.addUniformDirichletBoundaryCondition(
      ctx,
      {{"Boundary", "Xmax"},
       {"Component", 0},
       {"LoadingEvolution", [](const real t) noexcept { return 1e-5 * t; }}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(
      ctx,
      {{"Boundary", "Xmax"},
       {"Component", 1},
       {"LoadingEvolution", [](const real t) noexcept { return 2e-5 * t; }}}) |
      or_die;
  problem.setLinearSolver(ctx, "UMFPackSolver", {}) | or_die;
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  auto r = Results{};
  const auto& qs = m.getPartialQuadratureSpace();
  r.points_per_element =
      qs.getNumberOfIntegrationPoints() / qs.getNumberOfElements();
  auto& fes =
      problem.getFiniteElementDiscretization().getFiniteElementSpace<false>();
  auto& bi = problem.getBehaviourIntegrator(ctx, 1, 0) | or_die;
  r.rule_points =
      bi.getIntegrationRule(*(fes.GetFE(0)), *(fes.GetElementTransformation(0)))
          .GetNPoints();
  if (static_cast<size_type>(r.rule_points) != r.points_per_element) {
    // the integration would read the arrays of the quadrature space out of
    // their bounds
    check(false,
          "the behaviour integrator does not use the rule of its quadrature "
          "space");
    return r;
  }
  const auto o = problem.solve(ctx, 0, 1);
  check(o.status, "resolution failed (" + ctx.getErrorMessage() + ")");
  const auto& u = problem.getUnknowns(ets);
  const auto* const pu = u.HostRead();
  r.u.assign(pu, pu + u.Size());
  return r;
}  // end of compute

//! \return the relative difference of two vectors
static mfem_mgis::real getRelativeDifference(
    const std::vector<mfem_mgis::real>& a,
    const std::vector<mfem_mgis::real>& b) {
  using namespace mfem_mgis;
  if ((a.size() != b.size()) || (a.empty())) {
    return 1;
  }
  auto d = real{};
  auto n = real{};
  for (std::size_t i = 0; i != a.size(); ++i) {
    d = std::max(d, std::abs(a[i] - b[i]));
    n = std::max(n, std::abs(a[i]));
  }
  return d / n;
}  // end of getRelativeDifference

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* mesh = "cube.mesh";
  const char* library = "src/libBehaviourTest.so";
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&mesh, "-m", "--mesh", "mesh of the unit cube (hexahedra)");
  args.AddOption(&library, "-l", "--library", "behaviour library");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  // default rule: 27 points on the hexahedra of order 1
  const auto r = compute(mesh, library, std::nullopt, false);
  check((r.points_per_element == 27) && (r.rule_points == 27),
        "default rule: " + std::to_string(r.points_per_element) +
            " points per element, " + std::to_string(r.rule_points) +
            " points in the rule of the behaviour integrator");
  // rule of order 3: 8 points, same solution on affine hexahedra, point by
  // point and by kernels
  for (const auto kernels : {false, true}) {
    const auto name = std::string{kernels ? "kernels" : "point by point"};
    const auto r3 = compute(mesh, library, 3, kernels);
    check((r3.points_per_element == 8) && (r3.rule_points == 8),
          "rule of order 3, " + name + ": " +
              std::to_string(r3.points_per_element) + " points per element, " +
              std::to_string(r3.rule_points) +
              " points in the rule of the behaviour integrator");
    const auto d = getRelativeDifference(r.u, r3.u);
    std::cout << "rule of order 3, " << name
              << ": relative difference of the unknowns " << d << '\n';
    check(d < 1e-10, "rule of order 3, " + name + ": different solution");
  }
  // a negative order is refused
  {
    auto ctx = Context{};
    const auto p = construct<NonLinearEvolutionProblem>(
        ctx, getParameters(mesh, -1, false));
    const auto msg = std::string{ctx.getErrorMessage()};
    check(isInvalid(p) && (msg.find("invalid order of the integration rules") !=
                           std::string::npos),
          "negative order not refused (" + msg + ")");
  }
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
