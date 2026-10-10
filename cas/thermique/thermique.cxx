/*!
 * \file   thermique.cxx
 * \brief  Stationary non linear heat transfer in a beam (conductivity
 * depending on the temperature, temperatures imposed at both ends), solved
 * by the Newton solver of mfem-mgis, point by point or by kernels, on the
 * host or on the GPU.
 */

#include <cmath>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>
#include <string_view>
#include <vector>
#include "mfem/general/optparser.hpp"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* mesh = "beam-tet.mesh";
  const char* library = "src/libBehaviour.so";
  const char* linear_solver = "gmres-amg";
  const char* behaviour_integration = "";
  const char* adaptive_linear_solver_tolerance = "None";
  auto parallel = true;
  auto refinements = 3;
  auto steps = 5;
  auto max_iterations = 50;
  auto profiling = false;
  auto heating = real{600};
  auto relative_tolerance = real{1e-8};
  auto linear_tolerance = real{1e-10};
  auto expected_norm = real{-1};
  auto norm_tolerance = real{1e-6};
  mfem::OptionsParser args(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&mesh, "-m", "--mesh", "Mesh file.");
  args.AddOption(&library, "-l", "--library", "Behaviour library.");
  args.AddOption(&linear_solver, "-ls", "--linear-solver",
                 "Linear solver: mumps (parallel), umfpack (sequential), "
                 "gmres-amg (parallel, GMRES of hypre and BoomerAMG) or "
                 "gmres-jacobi (GMRES of MFEM and Jacobi preconditioner, "
                 "usable with the tangent operator applied without "
                 "assembly). The tangent operator of the heat transfer "
                 "behaviours is negative definite and the diagonal of the "
                 "imposed temperatures is set to one: the conjugate gradient "
                 "does not apply.");
  args.AddOption(&behaviour_integration, "-bi", "--behaviour-integration",
                 "Integration of the behaviour: PointByPoint or Batch. By "
                 "default, the value given by the environment variable "
                 "MFEM_MGIS_BEHAVIOUR_INTEGRATION, or PointByPoint.");
  args.AddOption(&adaptive_linear_solver_tolerance, "-alst",
                 "--adaptive-linear-solver-tolerance",
                 "Relative tolerance of the iterative linear solver adapted "
                 "to the decrease of the residual: None, EisenstatWalker1 or "
                 "EisenstatWalker2.");
  args.AddOption(&linear_tolerance, "-lt", "--linear-tolerance",
                 "Relative tolerance of the iterative linear solver.");
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "Parallel computation.");
  args.AddOption(&refinements, "-r", "--refinements",
                 "Number of uniform refinements of the mesh.");
  args.AddOption(&steps, "-n", "--steps", "Number of time steps.");
  args.AddOption(&max_iterations, "-mi", "--max-iterations",
                 "Maximum number of Newton iterations per step.");
  args.AddOption(&heating, "-dT", "--heating",
                 "Increase of the temperature imposed at the end of the "
                 "beam.");
  args.AddOption(&relative_tolerance, "-rt", "--relative-tolerance",
                 "Relative tolerance of the Newton solver.");
  args.AddOption(&expected_norm, "-en", "--expected-norm",
                 "Expected norm of the final temperature (no check if "
                 "negative).");
  args.AddOption(&norm_tolerance, "-nt", "--norm-tolerance",
                 "Relative tolerance of the comparison to the expected norm.");
  args.AddOption(&profiling, "-prof", "--profiling", "-no-prof",
                 "--no-profiling",
                 "Print the time spent in each phase of the computation.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    abort(EXIT_FAILURE);
  }
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  ctx.enableProfiling(profiling);
  auto success = true;
  {
    // the other options of mfem-mgis (computations by kernels) are given by
    // their environment variables
    auto parameters = Parameters{{"MeshFileName", mesh},
                                 {"FiniteElementFamily", "H1"},
                                 {"FiniteElementOrder", 1},
                                 {"UnknownsSize", 1},
                                 {"NumberOfUniformRefinements", refinements},
                                 {"Hypothesis", "Tridimensional"},
                                 {"Parallel", parallel}};
    if (!std::string_view{behaviour_integration}.empty()) {
      parameters.insert(ctx, "BehaviourIntegration",
                        std::string{behaviour_integration}) |
          or_die;
    }
    auto problem =
        construct<NonLinearEvolutionProblem>(ctx, parameters) | or_die;
    problem.getUnknowns(bts) = 293.15;
    problem.getUnknowns(ets) = 293.15;
    // temperatures at the integration points of the two element attributes
    // of beam-tet.mesh, external state variables of the behaviour, computed
    // from the unknowns
    // (the buffers of the vectors are kept when the outer vector grows)
    auto temperatures = std::vector<std::vector<real>>{};
    for (const auto m : {1, 2}) {
      problem.addBehaviourIntegrator(ctx, "StationaryNonLinearHeatTransfer",
                                     m, library,
                                     "StationaryNonLinearHeatTransfer") |
          or_die;
      auto& mat = problem.getMaterial(ctx, m, 0) | or_die;
      for (auto* s : {&mat.s0, &mat.s1}) {
        auto& T = temperatures.emplace_back(mat.n, 293.15);
        mgis::behaviour::setExternalStateVariable(ctx, *s, "Temperature",
                                                  std::span<real>(T)) |
            or_die;
      }
    }
    // temperature imposed at both ends
    problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", 1},
              {"Component", 0},
              {"LoadingEvolution",
               [](const real) noexcept { return 293.15; }}}) |
        or_die;
    problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", 2},
              {"Component", 0},
              {"LoadingEvolution",
               [heating](const real t) noexcept {
                 return 293.15 + heating * t;
               }}}) |
        or_die;
    const auto ls = std::string{linear_solver};
    if (ls == "mumps") {
      problem.setLinearSolver(ctx, "MUMPSSolver", {}) | or_die;
    } else if (ls == "umfpack") {
      problem.setLinearSolver(ctx, "UMFPackSolver", {}) | or_die;
    } else if (ls == "gmres-amg") {
      problem.setLinearSolver(
          ctx, "HypreGMRES",
          {{"VerbosityLevel", 0},
           {"Tolerance", linear_tolerance},
           {"MaximumNumberOfIterations", 2000},
           {"Preconditioner",
            Parameters{{"Name", "HypreBoomerAMG"},
                       {"Options", Parameters{{"VerbosityLevel", 0}}}}}}) |
          or_die;
    } else if (ls == "gmres-jacobi") {
      problem.setLinearSolver(
          ctx, "GMRESSolver",
          {{"VerbosityLevel", 0},
           {"RelativeTolerance", linear_tolerance},
           {"AbsoluteTolerance", 0.},
           {"MaximumNumberOfIterations", 100000},
           {"KDim", 100},
           {"Preconditioner",
            Parameters{{"Name", "OperatorJacobiSmoother"}}}}) |
          or_die;
    } else {
      std::cerr << "invalid linear solver '" << ls << "'\n";
      abort(EXIT_FAILURE);
    }
    auto solver_parameters =
        Parameters{{"VerbosityLevel", 0},
                   {"RelativeTolerance", relative_tolerance},
                   {"AbsoluteTolerance", 0.},
                   {"MaximumNumberOfIterations", max_iterations}};
    if (std::string_view{adaptive_linear_solver_tolerance} != "None") {
      solver_parameters.insert(ctx, "AdaptiveLinearSolverTolerance",
                               adaptive_linear_solver_tolerance) |
          or_die;
    }
    problem.setSolverParameters(ctx, solver_parameters) | or_die;
    // loop over the time steps
    const auto dt = real{1} / steps;
    auto t = real{0};
    auto iterations = std::vector<int>{};
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i != steps; ++i) {
      const auto r = problem.solve(ctx, t, dt);
      if (!r.status) {
        std::cerr << "step " << i + 1 << ": not converged\n";
        abort(EXIT_FAILURE);
      }
      iterations.push_back(r.iterations);
      problem.update(ctx) | or_die;
      t += dt;
    }
    auto total =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    // norm of the final temperature, to compare the solvers
    const auto& u = problem.getUnknowns(bts);
    auto unorm = mfem::InnerProduct(u, u);
#ifdef MFEM_USE_MPI
    if (parallel) {
      MPI_Allreduce(MPI_IN_PLACE, &unorm, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);
    }
#endif /* MFEM_USE_MPI */
    const auto& fed = problem.getFiniteElementDiscretization();
    const auto ndofs =
        fed.describesAParallelComputation()
            ? static_cast<long long>(
                  fed.getFiniteElementSpace<true>().GlobalTrueVSize())
            : static_cast<long long>(
                  fed.getFiniteElementSpace<false>().GetTrueVSize());
    if (isMainProcess(fed)) {
      auto sum = 0;
      std::cout << "dofs " << ndofs << " | ls " << ls << " | ||T|| "
                << std::setprecision(10) << std::sqrt(unorm)
                << " | iterations per step:";
      for (const auto n : iterations) {
        std::cout << ' ' << n;
        sum += n;
      }
      std::cout << " | total " << sum << '\n'
                << std::fixed << std::setprecision(3) << "time " << total
                << " s\n";
    }
    if ((expected_norm >= 0) && (std::abs(std::sqrt(unorm) - expected_norm) >
                                 norm_tolerance * expected_norm)) {
      if (isMainProcess(fed)) {
        std::cerr << "the norm of the solution differs from the expected one ("
                  << std::setprecision(10) << expected_norm << ")\n";
      }
      success = false;
    }
    if (profiling) {
      Profiler::OutputManager::printTimeTable(ctx);
    }
  }  // the problem is destroyed before the finalization of MPI
  finalize();
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
