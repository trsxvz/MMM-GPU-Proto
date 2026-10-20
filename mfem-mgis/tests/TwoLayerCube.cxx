/*!
 * \file   tests/TwoLayerCube.cxx
 * \brief
 * This example models a periodic unit cube made of two layers, split at
 * x = 0.5, under an imposed macroscopic strain. The solution is compared to
 * the analytical solution of this case, for the loading case selected by the
 * --test-case option.
 *
 * The cube is meshed by cube_2mat_per.mesh (4x4x4 hexahedra) and, more
 * finely, by Box.med (8x8x8 hexahedra), whose periodicity is described by
 * Box.per. Reading Box.med requires MFEM built with MED support:
 *
 *   --mesh Box.med
 *
 * Mechanical strain:
 *                 eps = E + grad_s v
 *
 *           with  E the given macroscopic strain
 *                 v the periodic displacement fluctuation
 * Displacement:
 *                   u = U + v
 *
 *           with  U the given displacement associated to E
 *                   E = grad_s U
 * The local microscopic strain is equal, on average, to the macroscopic strain:
 *           <eps> = <E>
 * \author Thomas Helfer, Guillaume Latu
 * \date   14/10/2020
 */

#include <limits>
#include <string>
#include <memory>
#include <string_view>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/ExecutionSettings.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "Tolerance.hxx"

constexpr double xmax = 1.;

void (*getSolution(const std::size_t i))(mfem::Vector&, const mfem::Vector&) {
  constexpr const auto xthr = xmax / 2.;
  std::array<void (*)(mfem::Vector&, const mfem::Vector&), 6u> solutions = {
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(1) = gradx * x(0);
        } else {
          u(1) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(2) = gradx * x(0);
        } else {
          u(2) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector&) { u = mfem_mgis::real{}; }};
  return solutions[i];
}

[[nodiscard]] static bool setLinearSolver(
    mfem_mgis::Context& ctx,
    mfem_mgis::AbstractNonLinearEvolutionProblem& p,
    const std::string_view s) noexcept {
  const auto ilu = mfem_mgis::Parameters{
      {"Name", "HypreILU"},
      {"Options",
       mfem_mgis::Parameters{
           // ILU(0) when hypre computes on a GPU, the only level of fill
           // supported there
           {"HypreILULevelOfFill", mfem_mgis::isHypreUsingGPU() ? 0 : 1}}}};
  const auto diagscale = mfem_mgis::Parameters{
      {"Name", "HypreDiagScale"},
      {"Options", mfem_mgis::Parameters{{"VerbosityLevel", 0}}}};
  if (s == "GMRESSolver") {
    return p.setLinearSolver(ctx, "GMRESSolver",
                             {{"VerbosityLevel", 1},
                              {"AbsoluteTolerance", 1e-12},
                              {"RelativeTolerance", 1e-12},
                              {"MaximumNumberOfIterations", 5000}});
  } else if (s == "CGSolver") {
    return p.setLinearSolver(ctx, "CGSolver",
                             {{"VerbosityLevel", 1},
                              {"AbsoluteTolerance", 1e-12},
                              {"RelativeTolerance", 1e-12},
                              {"MaximumNumberOfIterations", 5000}});
#ifdef MFEM_MGIS_HAVE_SUITESPARSE
  } else if (s == "UMFPackSolver") {
    return p.setLinearSolver(ctx, "UMFPackSolver", {});
#endif
#ifdef MFEM_MGIS_HAVE_MUMPS
  } else if (s == "MUMPSSolver") {
    return p.setLinearSolver(ctx, "MUMPSSolver",
                             {{"Symmetric", true}, {"PositiveDefinite", true}});
#endif
  } else if (s == "HypreFGMRES") {
    return p.setLinearSolver(ctx, "HypreFGMRES",
                             {{"VerbosityLevel", 1},
                              {"Tolerance", 1e-12},
                              {"Preconditioner", ilu},
                              {"MaximumNumberOfIterations", 5000}});
  } else if (s == "HyprePCG") {
    return p.setLinearSolver(ctx, "HyprePCG",
                             {{"VerbosityLevel", 1},
                              {"Tolerance", 1e-12},
                              {"Preconditioner", diagscale},
                              {"MaximumNumberOfIterations", 5000}});
  } else if (s == "HypreGMRES") {
    return p.setLinearSolver(ctx, "HypreGMRES",
                             {{"VerbosityLevel", 1},
                              {"Tolerance", 1e-12},
                              {"MaximumNumberOfIterations", 5000}});
  }
  return ctx.registerErrorMessage("unsupported linear solver '" +
                                  std::string{s} + "'");
}

static bool setSolverParameters(
    mfem_mgis::Context& ctx,
    mfem_mgis::AbstractNonLinearEvolutionProblem& problem) {
  // the absolute tolerance is the criterion of the case 5, whose initial
  // residual is only made of roundings: in float, they are about 4e-6, the
  // residuals of the other cases being about 1e2
  return problem.setSolverParameters(
      ctx,
      {{"VerbosityLevel", 0},
       {"RelativeTolerance", 1e-12},
       {"AbsoluteTolerance", mfem_mgis::unit_tests::tolerance(1e-12, 1e-4)},
       {"MaximumNumberOfIterations", 10}});
}  // end of setSolverParmeters

/*!
 * \return the threshold of the L2 error to the analytical solution: 1e-7
 * with double; with float, about ten times the largest error measured with
 * the given linear solver over the cases, sequential and parallel (the L2
 * norms of the solutions are 0.038 to 0.096). The relative tolerances of the
 * iterative solvers are bounded to 1.2e-5 by mfem-mgis in float. A solver
 * which was not measured is given the largest threshold.
 * \param[in] s: linear solver
 */
static mfem_mgis::real getThreshold(const std::string_view s) noexcept {
  // largest errors measured in float: 1.2e-4 (HypreFGMRES), 3.7e-5
  // (GMRESSolver, HypreGMRES), 3.5e-6 (CGSolver, HyprePCG, MUMPSSolver)
  if ((s == "CGSolver") || (s == "HyprePCG") || (s == "MUMPSSolver") ||
      (s == "UMFPackSolver")) {
    return mfem_mgis::unit_tests::tolerance(1e-7, 4e-5);
  }
  if ((s == "GMRESSolver") || (s == "HypreGMRES")) {
    return mfem_mgis::unit_tests::tolerance(1e-7, 4e-4);
  }
  return mfem_mgis::unit_tests::tolerance(1e-7, 1.2e-3);
}  // end of getThreshold

/*!
 * \brief prints the L2 error to the analytical solution, without comparing
 * it to a threshold
 * \return true, or an empty value on error
 */
std::optional<bool> printError(mfem_mgis::Context& ctx,
                               mfem_mgis::NonLinearEvolutionProblem& problem,
                               const std::size_t i) {
  const auto ob = mfem_mgis::compareToAnalyticalSolution(
      ctx, problem, getSolution(i),
      {{"VerbosityLevel", 2},
       {"CriterionThreshold", std::numeric_limits<mfem_mgis::real>::max()}});
  if (mfem_mgis::isInvalid(ob)) {
    return {};
  }
  return true;
}

std::optional<bool> checkSolution(mfem_mgis::Context& ctx,
                                  mfem_mgis::NonLinearEvolutionProblem& problem,
                                  const std::size_t i,
                                  const std::string_view s) {
  // the L2 error is printed (verbosity level 2)
  const auto ob = mfem_mgis::compareToAnalyticalSolution(
      ctx, problem, getSolution(i),
      {{"VerbosityLevel", 2}, {"CriterionThreshold", getThreshold(s)}});
  if (mfem_mgis::isInvalid(ob)) {
    return {};
  }
  if (!(*ob)) {
    mfem_mgis::getErrorStream() << "Error is greater than threshold\n";
    return false;
  }
  mfem_mgis::getErrorStream() << "Error is lower than threshold\n";
  return true;
}

/*!
 * \return the default behaviour library: the library built with the example
 * by `mfem_mgis_add_behaviours_library` for the processor of the behaviours
 * (option --execution), or else `src/libBehaviour.so`
 */
static std::string getDefaultBehaviourLibrary() {
#ifdef MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY
  return mfem_mgis::getBehavioursLibrary(
      MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY);
#else  /* MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY */
  return "src/libBehaviour.so";
#endif /* MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY */
}  // end of getDefaultBehaviourLibrary

struct TestParameters {
  const char* mesh_file = "cube_2mat_per.mesh";
  // behaviour library given by the user, if any
  const char* library = "";
  /*!
   * behaviour library used: the one given by the user, or else the default
   * one (see `getDefaultBehaviourLibrary`)
   */
  std::string behaviour_library;
  int order = 1;
  int refinement = 0;
  int tcase = 1;
  const char* linearsolver = "CGSolver";
  mfem_mgis::real xmax = 1.;
  mfem_mgis::real ymax = 1.;
  mfem_mgis::real zmax = 1.;
  bool parallel = true;
  bool check = true;
  //! \brief export the results for ParaView
  bool post_processing = true;
  //! \brief print the time table of the profiler at the end
  bool profiling = false;
};

TestParameters parseCommandLineOptions(int& argc, char* argv[]) {
  TestParameters p;

  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library",
                 "Material library. By default, the library built with the "
                 "example for the CPU or for the GPU, as the option "
                 "--execution places the behaviours.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "Number of uniform refinements of the mesh.");
  args.AddOption(&p.xmax, "-xm", "--xmax",
                 "x coordinate of the upper corner of the cube, which must "
                 "match the mesh.");
  args.AddOption(&p.ymax, "-ym", "--ymax",
                 "y coordinate of the upper corner of the cube, which must "
                 "match the mesh.");
  args.AddOption(&p.zmax, "-zm", "--zmax",
                 "z coordinate of the upper corner of the cube, which must "
                 "match the mesh.");
  args.AddOption(&p.tcase, "-t", "--test-case",
                 "identifier of the case : Exx->0, Eyy->1, Ezz->2, Exy->3, "
                 "Exz->4, Eyz->5");
  args.AddOption(&p.linearsolver, "-ls", "--linearsolver",
                 "Linear solver: GMRESSolver, CGSolver, UMFPackSolver "
                 "(sequential only), MUMPSSolver, HypreFGMRES with the "
                 "HypreILU preconditioner, HyprePCG with the HypreDiagScale "
                 "preconditioner or HypreGMRES (parallel only).");
  args.AddOption(&p.parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "Perform parallel computations.");
  args.AddOption(&p.check, "-c", "--check", "-no-c", "--no-check",
                 "Compare the solution to the analytical solution of the "
                 "two-layer cube, only valid for the provided meshes. "
                 "Without the comparison, the L2 error to the analytical "
                 "solution is still printed.");
  args.AddOption(&p.post_processing, "-pp", "--post-processing", "-no-pp",
                 "--no-post-processing",
                 "Export the results for ParaView.");
  args.AddOption(&p.profiling, "-prof", "--profiling", "-no-prof",
                 "--no-profiling",
                 "Print the time table of the profiler at the end.");
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
  args.PrintOptions(mfem_mgis::getOutputStream());
  p.behaviour_library = std::string_view{p.library}.empty()
                            ? getDefaultBehaviourLibrary()
                            : std::string{p.library};
  if ((p.tcase < 0) || (p.tcase > 5)) {
    mfem_mgis::getErrorStream() << "Invalid test case\n";
    mfem_mgis::abort(EXIT_FAILURE);
  }
  return p;
}

int executeMFEMMGISTest(mfem_mgis::Context& ctx, const TestParameters& p) {
  auto or_die = ctx.getFatalFailureHandler();
  constexpr const auto dim = mfem_mgis::size_type{3};
  // creating the finite element workspace

  auto fed =
      mfem_mgis::construct<mfem_mgis::FiniteElementDiscretization>(
          ctx,
          mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
                                {"FiniteElementFamily", "H1"},
                                {"FiniteElementOrder", p.order},
                                {"UnknownsSize", dim},
                                {"NumberOfUniformRefinements", p.refinement},
                                {"Parallel", p.parallel}}) |
      or_die;
  // building the non linear problem
  std::vector<mfem_mgis::real> corner1({0., 0., 0.});
  std::vector<mfem_mgis::real> corner2({p.xmax, p.ymax, p.zmax});
  auto problem =
      mfem_mgis::construct<mfem_mgis::PeriodicNonLinearEvolutionProblem>(
          ctx, fed, corner1, corner2) |
      or_die;
  problem.addBehaviourIntegrator(ctx, "Mechanics", 1, p.behaviour_library,
                                 "IsotropicLinearElasticity") |
      or_die;
  problem.addBehaviourIntegrator(ctx, "Mechanics", 2, p.behaviour_library,
                                 "IsotropicLinearElasticity") |
      or_die;
  // materials
  auto& m1 = problem.getMaterial(ctx, 1, 0) | or_die;
  auto& m2 = problem.getMaterial(ctx, 2, 0) | or_die;
  // setting the material properties
  auto set_properties = [&ctx, &or_die](auto& m, const double l,
                                        const double mu) {
    mgis::behaviour::setMaterialProperty(ctx, m.s0, "FirstLameCoefficient", l) |
        or_die;
    mgis::behaviour::setMaterialProperty(ctx, m.s0, "ShearModulus", mu) |
        or_die;
    mgis::behaviour::setMaterialProperty(ctx, m.s1, "FirstLameCoefficient", l) |
        or_die;
    mgis::behaviour::setMaterialProperty(ctx, m.s1, "ShearModulus", mu) |
        or_die;
  };

  std::array<mfem_mgis::real, 2> lambda({100, 200});
  std::array<mfem_mgis::real, 2> mu({75, 150});
  set_properties(m1, lambda[0], mu[0]);
  set_properties(m2, lambda[1], mu[1]);
  //
  auto set_temperature = [&ctx, &or_die](auto& m) {
    mgis::behaviour::setExternalStateVariable(ctx, m.s0, "Temperature",
                                              293.15) |
        or_die;
    mgis::behaviour::setExternalStateVariable(ctx, m.s1, "Temperature",
                                              293.15) |
        or_die;
  };
  set_temperature(m1);
  set_temperature(m2);

  // macroscopic strain
  std::vector<mfem_mgis::real> e(6, mfem_mgis::real{});
  if (p.tcase < 3) {
    e[p.tcase] = 1;
  } else {
    e[p.tcase] = 1.41421356237309504880 / 2;
  }
  problem.setMacroscopicGradientsEvolution([e](const double) { return e; });
  //
  setLinearSolver(ctx, problem, p.linearsolver) | or_die;
  setSolverParameters(ctx, problem) | or_die;

  // Add postprocessing and outputs
  if (p.post_processing) {
    problem.addPostProcessing(ctx, "ParaviewExportResults",
                              {{"OutputFileName", "TwoLayerCubeOutput-" +
                                                      std::to_string(p.tcase)}}) |
        or_die;
    std::vector<mfem_mgis::Parameter> materials_out{1, 2};
    problem.addPostProcessing(
        ctx, "ParaviewExportIntegrationPointResultsAtNodes",
        {{"OutputFileName",
          "TwoLayerCubeOutput-Strain-" + std::to_string(p.tcase)},
         {"Materials", {materials_out}},
         {"Results", "Strain"}}) |
        or_die;
    problem.addPostProcessing(
        ctx, "ParaviewExportIntegrationPointResultsAtNodes",
        {{"OutputFileName",
          "TwoLayerCubeOutput-Stress-" + std::to_string(p.tcase)},
         {"Materials", {materials_out}},
         {"Results", "Stress"}}) |
        or_die;
  }
  // solving the problem
  problem.solve(ctx, 0, 1) | or_die;
  problem.executePostProcessings(ctx, 0, 1) | or_die;
  // the L2 error to the analytical solution is printed, and compared to the
  // threshold of the linear solver if requested
  const auto b =
      p.check ? (checkSolution(ctx, problem, p.tcase, p.linearsolver) | or_die)
              : (printError(ctx, problem, p.tcase) | or_die);
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return b ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(int argc, char* argv[]) {
  auto ctx = mfem_mgis::Context{};
  mfem_mgis::initialize(argc, argv);
  const auto p = parseCommandLineOptions(argc, argv);
  ctx.enableProfiling(p.profiling);
  return executeMFEMMGISTest(ctx, p);
}
