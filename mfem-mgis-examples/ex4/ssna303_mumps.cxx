/*!
 * \file   ssna303_mumps.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#include <string>
#include <memory>
#include <string_view>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/linalg/hypre.hpp"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#ifdef MFEM_USE_MUMPS
#include "mfem/linalg/mumps.hpp"
#endif /* MFEM_USE_MUMPS */

#include "mfem/fem/datacollection.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/ExecutionSettings.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "CheckResultantForce.hxx"

int main(int argc, char** argv) {
  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  mfem_mgis::initialize(argc, argv);
  constexpr const auto dim = mfem_mgis::size_type{3};
  const char* mesh_file = "ssna303_3d.msh";
  const char* behaviour = "IsotropicLinearHardeningPlasticity";
  const char* library = "";
  bool use_fbar = false;
  // not null, since mfem::OptionsParser::PrintUsage stops at the first null
  // string
  const char* reference_file = "";
  const char* standard_reference_file = "";
  auto parallel = true;
  auto order = 1;
  auto refinement = 0;
  // order of the integration rules of the behaviour integrator, the default
  // rule of mfem-mgis if negative
  auto integration_rule_order = -1;
  auto nbsteps = 50;
  auto end_time = mfem_mgis::real{1};
  auto profiling = false;

  // options treatment
  mfem::OptionsParser args(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "Perform parallel computations.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&refinement, "-r", "--refinement",
                 "Number of uniform refinements of the mesh (the reference "
                 "files are the ones of the original mesh).");
  args.AddOption(&integration_rule_order, "-iro", "--integration-rule-order",
                 "Order of the integration rules (for instance 3 for 8 points "
                 "by element instead of 27), default rule if negative.");
  args.AddOption(&nbsteps, "-ns", "--nbsteps", "Number of time steps.");
  args.AddOption(
      &end_time, "-et", "--end-time",
      "End time. The displacement of the upper boundary is 6e-3 * t.");
  args.AddOption(&reference_file, "-rf", "--reference-file",
                 "Reference values of the resultant force on the upper "
                 "boundary, no comparison if empty.");
  args.AddOption(&profiling, "-prof", "--profiling", "-no-prof",
                 "--no-profiling",
                 "Print the time spent in each phase of the computation.");
#ifdef MGIS_HAVE_TFEL
  args.AddOption(&use_fbar, "-fb", "--use-fbar", "-no-fb", "--no-use-fbar",
                 "Use the FBar formulation.");
  args.AddOption(&standard_reference_file, "-srf", "--standard-reference-file",
                 "Reference values of the resultant force on the upper "
                 "boundary computed without FBar, compared with a larger "
                 "tolerance, no comparison if empty.");
#endif /* MGIS_HAVE_TFEL */
  args.AddOption(&library, "-l", "--library",
                 "Behaviour library. By default, the library built with the "
                 "example for the CPU or for the GPU, as the option "
                 "--execution places the behaviours.");
  args.Parse();
  if (args.Help()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::finalize();
    return EXIT_SUCCESS;
  }
  if (!args.Good()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  args.PrintOptions(mfem_mgis::getOutputStream());
  // behaviour library: given by --library, or else the library built with
  // the example for the processor of the behaviours (option --execution)
  const auto behaviour_library =
      std::string_view{library}.empty()
          ? mfem_mgis::getBehavioursLibrary(
                MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY)
          : std::string{library};
  ctx.enableProfiling(profiling);
  const auto* const output_file = use_fbar ? "force-fbar.txt" : "force.txt";
  // the non linear problem
  auto parameters = mfem_mgis::Parameters{
      {"MeshFileName", mesh_file},
      {"Materials", mfem_mgis::Parameters{{"NotchedBeam", 1}}},
      {"FiniteElementFamily", "H1"},
      {"FiniteElementOrder", order},
      {"NumberOfUniformRefinements", refinement},
      {"UnknownsSize", dim},
      {"Hypothesis", "Tridimensional"},
      {"Parallel", parallel}};
  if (integration_rule_order >= 0) {
    parameters.insert(ctx, "IntegrationRuleOrder", integration_rule_order) |
        or_die;
  }
  auto problem = mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
                     ctx, parameters) |
                 or_die;

  // 2 1 "Volume"
#ifdef MGIS_HAVE_TFEL
  if (use_fbar) {
    problem.addBehaviourIntegrator(
        ctx, "Mechanics", "NotchedBeam", behaviour_library, behaviour,
        {{"Regularization",
          mfem_mgis::Parameters{{"FBar", mfem_mgis::Parameters{}}}}}) |
        or_die;
  } else {
    problem.addBehaviourIntegrator(ctx, "Mechanics", "NotchedBeam",
                                   behaviour_library, behaviour) |
        or_die;
  }
#else  /* MGIS_HAVE_TFEL */
  problem.addBehaviourIntegrator(ctx, "Mechanics", "NotchedBeam",
                                 behaviour_library, behaviour) |
      or_die;
#endif /* MGIS_HAVE_TFEL */
  // materials
  auto& m1 = problem.getMaterial(ctx, "NotchedBeam", 0) | or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s0, "Temperature", 293.15) |
      or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s1, "Temperature", 293.15) |
      or_die;

  // boundary conditions

  // 3 LowerBoundary
  problem.addBoundaryCondition(
      ctx, mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 3, 1) |
               or_die) |
      or_die;
  // 4 SymmetryPlane1
  problem.addBoundaryCondition(
      ctx, mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 4, 0) |
               or_die) |
      or_die;
  // 5 SymmetryPlane2
  problem.addBoundaryCondition(
      ctx, mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 5, 2) |
               or_die) |
      or_die;
  // 2 UpperBoundary
  problem.addBoundaryCondition(
      ctx, mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 2, 1,
               [](const auto t) {
                 const auto u = 6e-3 * t;
                 return u;
               }) |
               or_die) |
      or_die;

  // solving the problem without petsc
  if (!mfem_mgis::usePETSc()) {
    // the default prediction concentrates the increment of the imposed
    // displacement in the elements next to the upper boundary
    problem.setPredictionPolicy(
        {.strategy =
             mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});
    problem.setSolverParameters(ctx, {{"VerbosityLevel", 2},
                                      {"RelativeTolerance", 1e-6},
                                      {"AbsoluteTolerance", 0.},
                                      {"MaximumNumberOfIterations", 10}}) |
        or_die;
    if (parallel) {
      problem.setLinearSolver(ctx, "MUMPSSolver", {}) | or_die;
    } else {
      problem.setLinearSolver(ctx, "UMFPackSolver", {}) | or_die;
    }
  }

  // vtk export
  problem.addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName", std::string("ssna303-displacements")}}) |
      or_die;
  problem.addPostProcessing(
      ctx, "ComputeResultantForceOnBoundary",
      {{"Boundary", 2}, {"OutputFileName", output_file}}) |
      or_die;

  // loop over time step
  const auto nsteps = mfem_mgis::size_type(nbsteps);
  const auto dt = end_time / nsteps;
  auto t = mfem_mgis::real{0};
  auto iteration = mfem_mgis::size_type{};
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    std::cout << "iteration " << iteration << " from " << t << " to " << t + dt
              << '\n';
    // resolution
    auto ct = t;
    auto dt2 = dt;
    auto nsteps = mfem_mgis::size_type{1};
    auto niter = mfem_mgis::size_type{0};
    while (nsteps != 0) {
      bool converged = problem.solve(ctx, ct, dt2);
      if (converged) {
        --nsteps;
        ct += dt2;
        if (nsteps == 0) {
          problem.executePostProcessings(ctx, t, dt) | or_die;
        }
        problem.update(ctx) | or_die;
      } else {
        // reason of the failure, if reported (the messages are cleared)
        const auto msg = ctx.getErrorMessage();
        std::cout << "\nsubstep: " << niter;
        if (!msg.empty()) {
          std::cout << " (" << msg << ")";
        }
        std::cout << '\n';
        nsteps *= 2;
        dt2 /= 2;
        ++niter;
        problem.revert(ctx) | or_die;
        if (niter == 10) {
          mfem_mgis::abort("maximum number of substeps");
        }
      }
    }
    t += dt;
    ++iteration;
    std::cout << '\n';
  }
  // comparison to the reference values, only on the process writing the
  // resultant force. The tolerance is above the rounding of the forces, which
  // are written with 6 significant digits. It is larger for the reference
  // values computed without FBar, which are only close to the results with
  // FBar.
  if (mfem_mgis::isMainProcess(problem.getFiniteElementDiscretization())) {
    if ((!std::string_view{reference_file}.empty()) &&
        (!checkVerticalForce(output_file, reference_file, 1e-4))) {
      return EXIT_FAILURE;
    }
    if ((!std::string_view{standard_reference_file}.empty()) &&
        (!checkVerticalForce(output_file, standard_reference_file, 1e-3))) {
      return EXIT_FAILURE;
    }
  }
  if (profiling) {
    mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  }
  return EXIT_SUCCESS;
}
