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

#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
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
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"
#include "MFEMMGIS/BehaviourIntegratorFactory.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
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
  // all the execution modes, unless hypre is built without GPU support and
  // the linear solvers are placed on the GPU, where the solvers of MFEM are
  // used
  if ((s.mode != ExecutionMode::UNSPECIFIED) &&
      (s.linear_solver == ExecutionLocation::GPU) && (!isHypreUsingGPU())) {
    problem.setLinearSolver(
        ctx, "CGSolver",
        {{"VerbosityLevel", 0},
         {"RelativeTolerance", 1e-14},
         {"AbsoluteTolerance", 0.},
         {"MaximumNumberOfIterations", 10000},
         {"Preconditioner", Parameters{{"Name", "OperatorJacobiSmoother"}}}}) |
        or_die;
  } else {
    problem.setLinearSolver(
        ctx, "HyprePCG",
        {{"VerbosityLevel", 0},
         {"Tolerance", 1e-14},
         {"MaximumNumberOfIterations", 1000},
         {"Preconditioner", Parameters{{"Name", "HypreBoomerAMG"}}}}) |
        or_die;
  }
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  const auto t0 = getMaterialArraysTransfers();
  const auto o = problem.solve(ctx, 0, 1);
  check(o.status, "resolution failed (" + ctx.getErrorMessage() + ")");
  // copies of the arrays of the material during the resolution: none if
  // the finite elements and the behaviours are computed on the host, in
  // both directions if they are computed on different processors
  const auto t1 = getMaterialArraysTransfers();
  const auto h2d = t1.host_to_device_copies - t0.host_to_device_copies;
  const auto d2h = t1.device_to_host_copies - t0.device_to_host_copies;
  std::cout << "copies of the arrays of the material: " << h2d
            << " from the host to the device, " << d2h
            << " from the device to the host\n";
  // without execution mode, the computations may be placed by the
  // environment (test suites): their locations are not checked
  const auto specified = s.mode != ExecutionMode::UNSPECIFIED;
  const auto fe_on_gpu =
      specified && (s.finite_elements == ExecutionLocation::GPU);
  const auto bh_on_gpu = specified && (s.behaviours == ExecutionLocation::GPU);
  // the debug device of MFEM, accepted with --execution cpu, copies the
  // arrays to its own memory
  if ((specified) && (!fe_on_gpu) && (!bh_on_gpu) &&
      (!mfem::Device::Allows(mfem::Backend::DEBUG_DEVICE))) {
    check((h2d == 0) && (d2h == 0),
          "copies of the arrays with the computations on the host");
  }
  if ((specified) && (fe_on_gpu != bh_on_gpu)) {
    check((h2d != 0) && (d2h != 0),
          "no copies of the arrays between the finite elements and the "
          "behaviours computed on different processors");
  }
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
  // with --execution cpu, the kernels compute on the host
  if (specified) {
    check(fe_kernels == (fe_on_gpu || (s.mode == ExecutionMode::CPU)),
          std::string{"finite elements computed "} +
              (fe_kernels ? "by kernels" : "element by element"));
  }
  // location of the integration of the behaviours (on a GPU)
  const auto bh_gpu =
      isIntegratedByBatchesOnDevice(m.b) &&
      mfem::Device::Allows(mfem::Backend::CUDA_MASK | mfem::Backend::HIP_MASK);
  if (specified) {
    check(bh_gpu == bh_on_gpu, "behaviour integrated on the wrong processor");
  }
  // location of the post-processings: the measure of the material
  const auto n0 = getNumberOfPostProcessingsComputationsByKernels();
  auto& bi = problem.getBehaviourIntegrator(ctx, 1, 0) | or_die;
  const auto measure = computeMeasure(bi);
  check(std::abs(measure - 1) < 1e-12,
        "measure of the cube: " + std::to_string(measure));
  const auto pp_kernels =
      getNumberOfPostProcessingsComputationsByKernels() != n0;
  if (specified) {
    const auto pp_gpu = s.post_processings == ExecutionLocation::GPU;
    check(pp_kernels == pp_gpu,
          std::string{"post-processings computed "} +
              (pp_kernels ? "by kernels" : "element by element"));
  }
}  // end of computeUniaxialTension

/*!
 * \brief a behaviour integrator which does not support the integration by
 * batches, the default of `AbstractBehaviourIntegrator` (all the behaviour
 * integrators of mfem-mgis support it), whose computations throw
 */
struct NoBatchIntegrationBehaviourIntegrator final
    : mfem_mgis::BehaviourIntegratorBase {
  NoBatchIntegrationBehaviourIntegrator(
      const mfem_mgis::FiniteElementDiscretization& fed,
      const mfem_mgis::size_type m,
      std::unique_ptr<const mfem_mgis::Behaviour> behaviour)
      : mfem_mgis::BehaviourIntegratorBase(
            std::make_shared<mfem_mgis::PartialQuadratureSpace>(
                fed,
                mfem_mgis::MaterialIdentifier{.id = m},
                selectIntegrationRule),
            std::move(behaviour)) {}
  static const mfem::IntegrationRule& selectIntegrationRule(
      const mfem::FiniteElement& e, const mfem::ElementTransformation&) {
    return mfem::IntRules.Get(e.GetGeomType(), 2);
  }
  const mfem::IntegrationRule& getIntegrationRule(
      const mfem::FiniteElement& e,
      const mfem::ElementTransformation& tr) const override {
    return selectIntegrationRule(e, tr);
  }
  mfem_mgis::real getIntegrationPointWeight(
      mfem::ElementTransformation& tr,
      const mfem::IntegrationPoint& ip) const noexcept override {
    return ip.weight * tr.Weight();
  }
  bool integrate(const mfem::FiniteElement&,
                 mfem::ElementTransformation&,
                 const mfem::Vector&,
                 const mfem_mgis::IntegrationType) override {
    throw std::runtime_error("NoBatchIntegrationBehaviourIntegrator: unused");
  }
  void updateResidual(mfem::Vector&,
                      const mfem::FiniteElement&,
                      mfem::ElementTransformation&,
                      const mfem::Vector&) override {
    throw std::runtime_error("NoBatchIntegrationBehaviourIntegrator: unused");
  }
  void updateJacobian(mfem::DenseMatrix&,
                      const mfem::FiniteElement&,
                      mfem::ElementTransformation&,
                      const mfem::Vector&) override {
    throw std::runtime_error("NoBatchIntegrationBehaviourIntegrator: unused");
  }
  void computeInnerForces(mfem::Vector&,
                          const mfem::FiniteElement&,
                          mfem::ElementTransformation&) override {
    throw std::runtime_error("NoBatchIntegrationBehaviourIntegrator: unused");
  }
};  // end of struct NoBatchIntegrationBehaviourIntegrator

/*!
 * \brief with an execution mode, the behaviours are integrated by batches,
 * unless the behaviour integrators do not support it: a behaviour
 * integrator without batch integration is integrated point by point (its
 * computations throw: the resolution fails after the checks of the
 * integration, which refuse it without this replacement)
 * \param[in] mesh: mesh of the square
 * \param[in] library: behaviour library
 */
static void checkBatchIntegrationFallback(const char* const mesh,
                                          const char* const library) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  BehaviourIntegratorFactory::get(Hypothesis::PLANESTRAIN)
          .addGenerator(
              ctx, "NoBatchIntegration",
              [](Context&, const FiniteElementDiscretization& fed,
                 const size_type m, std::unique_ptr<const Behaviour> b,
                 const Parameters&)
                  -> std::unique_ptr<AbstractBehaviourIntegrator> {
                return std::make_unique<NoBatchIntegrationBehaviourIntegrator>(
                    fed, m, std::move(b));
              }) |
      or_die;
  auto problem = construct<NonLinearEvolutionProblem>(
                     ctx, Parameters{{"MeshFileName", mesh},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", 1},
                                     {"UnknownsSize", 2},
                                     {"Materials", Parameters{{"Square", 1}}},
                                     {"Hypothesis", "PlaneStrain"},
                                     {"Parallel", true}}) |
                 or_die;
  problem.addBehaviourIntegrator(ctx, "NoBatchIntegration", "Square", library,
                                 "IsotropicLinearHardeningPlasticity") |
      or_die;
  auto& m = problem.getMaterial(ctx, "Square", 0) | or_die;
  for (auto* st : {&m.s0, &m.s1}) {
    mgis::behaviour::setExternalStateVariable(ctx, *st, "Temperature", 293.15) |
        or_die;
  }
  auto c = Context{};
  const auto r = problem.solve(c, 0, 1);
  const auto msg = std::string{c.getErrorMessage()};
  check((!r.status) &&
            (msg.find("NoBatchIntegrationBehaviourIntegrator: unused") !=
             std::string::npos) &&
            (msg.find("does not support the integration of the behaviours by "
                      "batches") == std::string::npos),
        "behaviour integrator without batch integration: not integrated "
        "point by point (" +
            msg + ")");
}  // end of checkBatchIntegrationFallback

/*!
 * \brief with `--execution cpu`, the finite elements are computed by
 * kernels on the host, unless the behaviour integrators do not support
 * them: uniaxial tension (plane strain) of a square made of elements of
 * different geometries (different integration rules), computed element by
 * element
 * \param[in] mesh: mesh of the square
 * \param[in] library: behaviour library
 */
static void checkFiniteElementsKernelsFallback(const char* const mesh,
                                               const char* const library) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  auto problem = construct<NonLinearEvolutionProblem>(
                     ctx, Parameters{{"MeshFileName", mesh},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", 1},
                                     {"UnknownsSize", 2},
                                     {"Materials", Parameters{{"Square", 1}}},
                                     {"Boundaries", Parameters{{"Bottom", 1},
                                                               {"Right", 2},
                                                               {"Top", 3},
                                                               {"Left", 4}}},
                                     {"Hypothesis", "PlaneStrain"},
                                     {"Parallel", true}}) |
                 or_die;
  problem.addBehaviourIntegrator(ctx, "Mechanics", "Square", library,
                                 "IsotropicLinearHardeningPlasticity") |
      or_die;
  auto& m = problem.getMaterial(ctx, "Square", 0) | or_die;
  for (auto* st : {&m.s0, &m.s1}) {
    mgis::behaviour::setExternalStateVariable(ctx, *st, "Temperature", 293.15) |
        or_die;
  }
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Left"}, {"Component", 0}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Bottom"}, {"Component", 1}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Right"},
            {"Component", 0},
            {"LoadingEvolution",
             [](const real t) noexcept { return strain * t; }}}) |
      or_die;
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
  check(o.status, "elements of different geometries: resolution failed (" +
                      ctx.getErrorMessage() + ")");
  check(problem.getImplementation<true>()
                .getNumberOfTangentOperatorAssembliesByKernels() == 0,
        "elements of different geometries: computations by kernels");
  // plane strain: sxx = E e / (1 - nu * nu)
  m.synchronizeOnHostForReading();
  const auto& sig = m.s1.thermodynamic_forces;
  constexpr auto nu = mfem_mgis::real{0.34};
  const auto sxx = young * strain / (1 - nu * nu);
  auto dmax = real{};
  for (std::size_t i = 0; i < sig.size(); i += 4) {
    dmax = std::max(dmax, std::abs(sig[i] - sxx));
  }
  check(dmax < 1e-8 * sxx,
        "elements of different geometries: stress, "
        "largest difference " +
            std::to_string(dmax / sxx));
}  // end of checkFiniteElementsKernelsFallback

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
 * \param[in] preconditioner: name of the preconditioner, if any
 */
static void checkLinearSolver(const char* const mesh,
                              const char* const solver,
                              const bool refused,
                              const char* const preconditioner = nullptr) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto problem = makeProblem(ctx, mesh);
  auto parameters = Parameters{};
  if (preconditioner != nullptr) {
    parameters.insert(
        may_throw, {{"Preconditioner", Parameters{{"Name", preconditioner}}}});
  }
  const auto r = problem.setLinearSolver(ctx, solver, parameters);
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
  if (preconditioner != nullptr) {
    std::cout << "linear solver '" << solver << "' with preconditioner '"
              << preconditioner << "': " << (r ? "accepted" : "refused")
              << '\n';
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
  const char* behaviours = "ExecutionSettingsBehaviours";
  const char* host_library = "";
  const char* mixed_mesh = "";
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&expected, "-e", "--expected-settings",
                 "expected settings: none, or the mode followed by the "
                 "locations of the finite elements, of the behaviours, of "
                 "the linear solvers and of the post-processings, separated "
                 "by commas (for instance mixed,cpu,gpu,gpu,cpu)");
  args.AddOption(&mesh, "-m", "--mesh", "mesh of the unit cube (hexahedra)");
  args.AddOption(&behaviours, "-bd", "--behaviours-directory",
                 "directory of the behaviour libraries built by "
                 "mfem_mgis_add_behaviours_library (the library for the GPU "
                 "is not tested if it is not built)");
  args.AddOption(&mixed_mesh, "-mm", "--mixed-mesh",
                 "mesh of a square made of elements of different geometries "
                 "(not used if empty)");
  args.AddOption(&host_library, "-lh", "--host-library",
                 "library generated by the generic-parallel interface for "
                 "the host, used too if the behaviours are on the CPU (not "
                 "used if empty)");
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
  // library placed where the settings place the behaviours, and the
  // library of the other processor
  const auto bh_gpu = (s.mode != ExecutionMode::UNSPECIFIED) &&
                      (s.behaviours == ExecutionLocation::GPU);
  const auto library = getBehavioursLibrary(behaviours);
  const auto expected_library = std::string{behaviours} +
                                (bh_gpu ? "/gpu/" : "/cpu/") +
                                "libBehaviour.so";
  check(library == expected_library,
        "library '" + library + "' instead of '" + expected_library + "'");
  const auto cpu_library = std::string{behaviours} + "/cpu/libBehaviour.so";
  const auto gpu_library = std::string{behaviours} + "/gpu/libBehaviour.so";
  const auto has_gpu_library = std::filesystem::exists(gpu_library);
  if ((!bh_gpu) || (has_gpu_library)) {
    computeUniaxialTension(mesh, library.c_str());
  } else {
    check(false, "no library for the GPU");
  }
  // the libraries generated by the generic-parallel interface for the host,
  // integrated by batches on the CPU (the default with an execution mode)
  if ((s.mode != ExecutionMode::UNSPECIFIED) && (!bh_gpu) &&
      (std::strlen(host_library) != 0)) {
    computeUniaxialTension(mesh, host_library);
  }
  // computations element by element when the kernels are not supported
  if ((s.mode == ExecutionMode::CPU) && (std::strlen(mixed_mesh) != 0)) {
    checkFiniteElementsKernelsFallback(mixed_mesh, library.c_str());
    checkBatchIntegrationFallback(mixed_mesh, library.c_str());
  }
  // no preconditioner, for the solvers of hypre (where they are accepted)
  // and of MFEM
  const auto hypre_accepted = (s.mode == ExecutionMode::UNSPECIFIED) ||
                              (s.linear_solver == ExecutionLocation::CPU) ||
                              isHypreUsingGPU();
  if (hypre_accepted) {
    checkLinearSolver(mesh, "HyprePCG", false, "None");
  }
  if ((s.mode != ExecutionMode::MIXED) ||
      (s.linear_solver == ExecutionLocation::GPU)) {
    checkLinearSolver(mesh, "CGSolver", false, "None");
  }
  // libraries and linear solvers computing elsewhere are refused
  if (s.mode != ExecutionMode::UNSPECIFIED) {
    if (bh_gpu) {
      checkBehaviourLibraryRefusal(mesh, cpu_library.c_str());
    } else if (has_gpu_library) {
      checkBehaviourLibraryRefusal(mesh, gpu_library.c_str());
    }
#ifdef MFEM_USE_MUMPS
    checkLinearSolver(mesh, "MUMPSSolver",
                      s.linear_solver == ExecutionLocation::GPU);
#endif /* MFEM_USE_MUMPS */
    // the solvers of MFEM compute where the device of MFEM computes
    checkLinearSolver(mesh, "CGSolver",
                      (s.mode == ExecutionMode::MIXED) &&
                          (s.linear_solver == ExecutionLocation::CPU));
    // the preconditioner of hypre Euclid only computes on the host
    if (hypre_accepted) {
      checkLinearSolver(mesh, "HyprePCG",
                        s.linear_solver == ExecutionLocation::GPU,
                        "HypreEuclid");
    }
  }
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
