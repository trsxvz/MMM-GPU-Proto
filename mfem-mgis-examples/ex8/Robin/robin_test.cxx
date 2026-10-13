/*!
 * \file   robin_test.cxx
 * \brief
 * Steady state heat conduction in a bar of length L with a constant
 * conductivity k. The temperature T0 is imposed at x = 0 and a Robin boundary
 * condition, j.n = h (T - Tinf), is applied at x = L. The exact solution is
 * linear and exactly represented by the finite element solution:
 *
 *   T(x) = T0 + a x, with a = -h (T0 - Tinf) / (k + h L)
 */

#include <type_traits>
#include <span>
#include <vector>
#include <memory>
#include <string>
#include <cstdlib>
#include <string_view>

#include "mfem/general/optparser.hpp"
#include "mfem/mesh/mesh.hpp"
#include "mfem/mesh/pmesh.hpp"

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/ExecutionSettings.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"

#include "../headers/RobinBC.hxx"

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  // boundaries at x = 0 and x = L
  constexpr int left = 1;
  constexpr int right = 2;
  // length, conductivity (see StationaryLinearHeatTransfer.mfront), imposed
  // temperature, ambient temperature and heat transfer coefficient
  constexpr double L = 1.0;
  constexpr double k = 10.0;
  constexpr double T0 = 500.0;
  constexpr double Tinf = 300.0;
  constexpr double h = 100.0;
  constexpr double a = -h * (T0 - Tinf) / (k + h * L);
  // options
  int nx = 10;
  int order = 2;
  const char* library = "";
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&nx, "-nx", "--nx", "Number of elements along the bar.");
  args.AddOption(&order, "-o", "--order", "Finite element order.");
  args.AddOption(&library, "-l", "--library",
                 "Behaviour library. By default, the library built with the "
                 "test for the CPU or for the GPU, as the option --execution "
                 "places the behaviours.");
  // options of mfem-mgis, such as --execution, handled by
  // mfem_mgis::initialize
  mfem_mgis::declareDefaultOptions(args);
  args.Parse();
  if (args.Help()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::finalize();
    std::exit(EXIT_SUCCESS);
  }
  if (!args.Good()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  const auto behaviour_library =
      std::string_view{library}.empty()
          ? mfem_mgis::getBehavioursLibrary(
                MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY)
          : std::string{library};
  // mesh of the bar
  auto mesh = mfem::Mesh::MakeCartesian3D(nx, 2, 2, mfem::Element::HEXAHEDRON,
                                          L, 1.0, 1.0);
  for (int i = 0; i != mesh.GetNBE(); ++i) {
    // center of the boundary element
    auto ip = mfem::IntegrationPoint{};
    ip.Set2(0.5, 0.5);
    auto c = mfem::Vector(3);
    mesh.GetBdrElementTransformation(i)->Transform(ip, c);
    const auto b = (c[0] < 1e-10) ? left : ((c[0] > L - 1e-10) ? right : 3);
    mesh.GetBdrElement(i)->SetAttribute(b);
  }
  mesh.SetAttributes();
  auto fed = mfem_mgis::make_shared<mfem_mgis::FiniteElementDiscretization>(
                 ctx, std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh),
                 mfem_mgis::Parameters{{"FiniteElementFamily", "H1"},
                                       {"FiniteElementOrder", order},
                                       {"UnknownsSize", 1}}) |
             or_die;
  auto problem = mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
                     ctx, fed, mgis::behaviour::Hypothesis::TRIDIMENSIONAL) |
                 or_die;
  problem.getUnknowns(mfem_mgis::bts) = T0;
  problem.getUnknowns(mfem_mgis::ets) = T0;
  // material
  problem.addBehaviourIntegrator(ctx, "StationaryNonLinearHeatTransfer", 1,
                                 behaviour_library,
                                 "StationaryLinearHeatTransfer") |
      or_die;
  auto& m = problem.getMaterial(ctx, 1, 0) | or_die;
  // the temperature at the integration points must be stored externally
  auto temperatures = std::vector<mfem_mgis::real>(m.n, T0);
  for (auto& s : {&m.s0, &m.s1}) {
    mgis::behaviour::setExternalStateVariable(
        ctx, *s, "Temperature", std::span<mfem_mgis::real>(temperatures),
        mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
        mgis::behaviour::MaterialStateManager::NOUPDATE) |
        or_die;
  }
  // boundary conditions
  problem.addBoundaryCondition(
      ctx, mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, fed, left, 0, [](const auto) noexcept { return T0; }) |
               or_die) |
      or_die;
  problem.addBoundaryCondition(ctx, mfem_mgis::make_unique<mfem_mgis::RobinBC>(
                                        ctx, fed, right, h, Tinf, nullptr) |
                                        or_die) |
      or_die;
  // resolution
  problem.setLinearSolver(
      ctx, "HypreGMRES",
      {{"Tolerance", 1e-12},
       {"MaximumNumberOfIterations", 1000},
       {"Preconditioner", mfem_mgis::Parameters{{"Name", "HypreBoomerAMG"}}}}) |
      or_die;
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 1},
                                    {"RelativeTolerance", 1e-10},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  problem.solve(ctx, 0, 1) | or_die;
  // comparison to the exact solution. With float, the L2 error is 9.3e-4 K
  // (about 3e-6 of the temperatures)
  constexpr auto threshold =
      std::is_same_v<mfem_mgis::real, double> ? 1e-8 : 1e-2;
  const auto ok =
      mfem_mgis::compareToAnalyticalSolution(
          ctx, problem,
          [](mfem::Vector& T, const mfem::Vector& x) { T(0) = T0 + a * x(0); },
          {{"CriterionThreshold", threshold}}) |
      or_die;
  if (!ok) {
    mfem_mgis::getErrorStream()
        << "the temperature does not match the exact solution\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
