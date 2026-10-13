/*!
 * \file   ExecutionSettingsTest.cxx
 * \brief  Test of the execution settings (`--execution cpu|gpu|mixed` and
 * the locations of the phases in the mixed mode): the settings read by
 * `initialize`, the device of MFEM and the memory of hypre, and a uniaxial
 * tension of a cube, whose stress is known, computed with the default
 * values of the parameters given by the settings: the finite elements, the
 * behaviours and the post-processings must be computed where the settings
 * place them, and a behaviour library or a linear solver computing
 * elsewhere must be refused. The contradictory options are refused by
 * `initialize` (tests registered with the expected message).
 * \date   13/10/2026
 */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <string>
#include <string_view>
#include "mfem/general/device.hpp"
#include "mfem/general/optparser.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/ExecutionSettings.hxx"
#include "MFEMMGIS/PostProcessingsComputation.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

//! \return a description of the execution settings
static std::string describe(const mfem_mgis::ExecutionSettings& s) {
  using namespace mfem_mgis;
  const auto mode = [&s]() -> std::string {
    switch (s.mode) {
      case ExecutionMode::CPU:
        return "cpu";
      case ExecutionMode::GPU:
        return "gpu";
      case ExecutionMode::MIXED:
        return "mixed";
      case ExecutionMode::UNSPECIFIED:
        break;
    }
    return "none";
  }();
  if (s.mode == ExecutionMode::UNSPECIFIED) {
    return mode;
  }
  return mode + "," + std::string{getExecutionLocationName(s.finite_elements)} +
         "," + std::string{getExecutionLocationName(s.behaviours)} + "," +
         std::string{getExecutionLocationName(s.linear_solver)} + "," +
         std::string{getExecutionLocationName(s.post_processings)};
}  // end of describe

//! \brief Young modulus of the behaviour IsotropicLinearHardeningPlasticity
static constexpr mfem_mgis::real young = 70e9;
//! \brief imposed strain, in the elastic domain of the behaviour
static constexpr mfem_mgis::real strain = 1e-5;

/*!
 * \return the problem of the uniaxial tension of the unit cube, with the
 * default values of the parameters choosing the computations
 * \param[in, out] ctx: execution context
 * \param[in] mesh: mesh of the unit cube
 */
static mfem_mgis::NonLinearEvolutionProblem makeProblem(
    mfem_mgis::Context& ctx, const char* const mesh) {
  using namespace mfem_mgis;
  auto or_die = ctx.getFatalFailureHandler();
  return construct<NonLinearEvolutionProblem>(
             ctx, Parameters{{"MeshFileName", mesh},
                             {"FiniteElementFamily", "H1"},
                             {"FiniteElementOrder", 1},
                             {"UnknownsSize", 3},
                             {"NumberOfUniformRefinements", 1},
                             {"Materials", Parameters{{"Cube", 1}}},
                             {"Boundaries", Parameters{{"Ymin", 1},
                                                       {"Zmin", 2},
                                                       {"Xmax", 3},
                                                       {"Xmin", 5}}},
                             {"Hypothesis", "Tridimensional"},
                             {"Parallel", true}}) |
         or_die;
}  // end of makeProblem

/*!
 * \brief uniaxial tension of the unit cube, computed with the default
 * values of the parameters given by the execution settings, and checks of
 * the locations of the computations
 * \param[in] mesh: mesh of the unit cube
 * \param[in] library: behaviour library integrated where the settings place
 * the behaviours
 */
static void computeUniaxialTension(const char* const mesh,
                                   const char* const library) {
  using namespace mfem_mgis;
  const auto& s = getExecutionSettings();
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  auto problem = makeProblem(ctx, mesh);
  problem.addBehaviourIntegrator(ctx, "Mechanics", "Cube", library,
                                 "IsotropicLinearHardeningPlasticity") |
      or_die;
  auto& m = problem.getMaterial(ctx, "Cube", 0) | or_die;
  for (auto* st : {&m.s0, &m.s1}) {
    mgis::behaviour::setExternalStateVariable(ctx, *st, "Temperature", 293.15) |
        or_die;
  }
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Xmin"}, {"Component", 0}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Ymin"}, {"Component", 1}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Zmin"}, {"Component", 2}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Xmax"},
            {"Component", 0},
            {"LoadingEvolution",
             [](const real t) noexcept { return strain * t; }}}) |
      or_die;
  // the solvers of hypre compute where its memory is: they are accepted in
  // all the execution modes
  problem.setLinearSolver(
      ctx, "HyprePCG",
      {{"VerbosityLevel", 0},
       {"Tolerance", 1e-14},
       {"MaximumNumberOfIterations", 1000},
       {"Preconditioner", Parameters{{"Name", "HypreBoomerAMG"}}}}) |
      or_die;
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  const auto o = problem.solve(ctx, 0, 1);
  check(o.status, "resolution failed (" + ctx.getErrorMessage() + ")");
  // uniaxial stress state: sxx = E e, the other components are null
  m.synchronizeOnHostForReading();
  const auto& sig = m.s1.thermodynamic_forces;
  const auto sxx = young * strain;
  auto dmax = real{};
  for (std::size_t i = 0; i != sig.size(); ++i) {
    const auto expected = (i % 6 == 0) ? sxx : real{0};
    dmax = std::max(dmax, std::abs(sig[i] - expected));
  }
  check(dmax < 1e-8 * sxx, "stress: largest difference " +
                               std::to_string(dmax / sxx) +
                               " (relative to the axial stress)");
  // location of the computations of the finite elements
  const auto& implementation = problem.getImplementation<true>();
  const auto fe_kernels =
      implementation.getNumberOfTangentOperatorAssembliesByKernels() != 0;
  const auto fe_gpu = (s.mode != ExecutionMode::UNSPECIFIED) &&
                      (s.finite_elements == ExecutionLocation::GPU);
  check(fe_kernels == fe_gpu,
        std::string{"finite elements computed "} +
            (fe_kernels ? "by kernels" : "element by element"));
  // location of the integration of the behaviours
  const auto bh_gpu = (s.mode != ExecutionMode::UNSPECIFIED) &&
                      (s.behaviours == ExecutionLocation::GPU);
  check(isIntegratedByBatchesOnDevice(m.b) == bh_gpu,
        "behaviour integrated on the wrong processor");
  // location of the post-processings: the measure of the material
  const auto n0 = getNumberOfPostProcessingsComputationsByKernels();
  auto& bi = problem.getBehaviourIntegrator(ctx, 1, 0) | or_die;
  const auto measure = computeMeasure(bi);
  check(std::abs(measure - 1) < 1e-12,
        "measure of the cube: " + std::to_string(measure));
  const auto pp_kernels =
      getNumberOfPostProcessingsComputationsByKernels() != n0;
  const auto pp_gpu = (s.mode != ExecutionMode::UNSPECIFIED) &&
                      (s.post_processings == ExecutionLocation::GPU);
  check(pp_kernels == pp_gpu,
        std::string{"post-processings computed "} +
            (pp_kernels ? "by kernels" : "element by element"));
}  // end of computeUniaxialTension

/*!
 * \brief checks that a behaviour library integrated elsewhere than where
 * the settings place the behaviours is refused
 * \param[in] mesh: mesh of the unit cube
 * \param[in] library: behaviour library
 */
static void checkBehaviourLibraryRefusal(const char* const mesh,
                                         const char* const library) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto problem = makeProblem(ctx, mesh);
  const auto r = problem.addBehaviourIntegrator(
      ctx, "Mechanics", "Cube", library, "IsotropicLinearHardeningPlasticity");
  const auto msg = std::string{ctx.getErrorMessage()};
  check(isInvalid(r) &&
            (msg.find("while the execution settings place the behaviours") !=
             std::string::npos),
        "behaviour library '" + std::string{library} + "' not refused (" + msg +
            ")");
}  // end of checkBehaviourLibraryRefusal

/*!
 * \brief checks that a linear solver is refused, or accepted
 * \param[in] mesh: mesh of the unit cube
 * \param[in] solver: name of the linear solver
 * \param[in] refused: expected result
 */
static void checkLinearSolver(const char* const mesh,
                              const char* const solver,
                              const bool refused) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto problem = makeProblem(ctx, mesh);
  const auto r = problem.setLinearSolver(ctx, solver, {});
  const auto msg = std::string{ctx.getErrorMessage()};
  if (refused) {
    check((!r) && (msg.find("while the execution settings place the linear "
                            "solvers") != std::string::npos),
          "linear solver '" + std::string{solver} + "' not refused (" + msg +
              ")");
  } else {
    check(r,
          "linear solver '" + std::string{solver} + "' refused (" + msg + ")");
  }
}  // end of checkLinearSolver

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  try {
    initialize(argc, argv);
  } catch (std::exception& e) {
    // contradictory options: the message is checked by the test, and the
    // program stops at once (MPI and hypre may be partially initialized)
    std::cerr << "refused: " << e.what() << std::endl;
    std::_Exit(EXIT_FAILURE);
  }
  const char* expected = "none";
  const char* mesh = "cube.mesh";
  const char* library = "src/libBehaviourTest.so";
  const char* gpu_library = "";
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&expected, "-e", "--expected-settings",
                 "expected settings: none, or the mode followed by the "
                 "locations of the finite elements, of the behaviours, of "
                 "the linear solvers and of the post-processings, separated "
                 "by commas (for instance mixed,cpu,gpu,gpu,cpu)");
  args.AddOption(&mesh, "-m", "--mesh", "mesh of the unit cube (hexahedra)");
  args.AddOption(&library, "-l", "--library",
                 "behaviour library computing on the host");
  args.AddOption(&gpu_library, "-lg", "--gpu-library",
                 "behaviour library generated by the generic-parallel "
                 "interface for the GPU (not tested if empty)");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  const auto& s = getExecutionSettings();
  check(describe(s) == expected, "settings '" + describe(s) + "' instead of '" +
                                     std::string{expected} + "'");
  // device of MFEM and memory of hypre
  const auto gpu_device =
      mfem::Device::Allows(mfem::Backend::CUDA_MASK | mfem::Backend::HIP_MASK);
  if (s.mode == ExecutionMode::CPU) {
    check(!gpu_device, "the device uses a GPU with --execution cpu");
  } else if (s.mode != ExecutionMode::UNSPECIFIED) {
    check(gpu_device, "the device does not use a GPU");
#if defined(MFEM_USE_MPI) && defined(HYPRE_USING_GPU)
    check(isHypreUsingGPU() == (s.linear_solver == ExecutionLocation::GPU),
          "hypre does not compute where the linear solvers are placed");
#endif /* defined(MFEM_USE_MPI) && defined(HYPRE_USING_GPU) */
  }
  const auto has_gpu_library = std::strlen(gpu_library) != 0;
  const auto bh_gpu = (s.mode != ExecutionMode::UNSPECIFIED) &&
                      (s.behaviours == ExecutionLocation::GPU);
  if ((!bh_gpu) || (has_gpu_library)) {
    computeUniaxialTension(mesh, bh_gpu ? gpu_library : library);
  } else {
    std::cout << "no library for the GPU, computation not tested\n";
  }
  // libraries and linear solvers computing elsewhere are refused
  if (s.mode != ExecutionMode::UNSPECIFIED) {
    if (bh_gpu) {
      checkBehaviourLibraryRefusal(mesh, library);
    } else if (has_gpu_library) {
      checkBehaviourLibraryRefusal(mesh, gpu_library);
    }
#ifdef MFEM_USE_MUMPS
    checkLinearSolver(mesh, "MUMPSSolver",
                      s.linear_solver == ExecutionLocation::GPU);
#endif /* MFEM_USE_MUMPS */
    // the solvers of MFEM compute where the device of MFEM computes
    checkLinearSolver(mesh, "CGSolver",
                      (s.mode == ExecutionMode::MIXED) &&
                          (s.linear_solver == ExecutionLocation::CPU));
  }
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
