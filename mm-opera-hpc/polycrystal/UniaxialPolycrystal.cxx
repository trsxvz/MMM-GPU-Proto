/*!
 * \file   Mono_UO2_CosH_Jaco.cxx
 * \brief
 * This example is modelling ...
 * \author Raphael Prat, Maxence Wangermez
 * \date   06/2023
 */

#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <exception>
#include <csignal>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <functional>
#include <string_view>

#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Behaviour.hxx"
#include "MFEMMGIS/ExecutionSettings.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/RotationMatrix.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/MechanicalPostProcessings.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"
#include "MFEMMGIS/PartialQuadratureFunctionsSet.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/L2Projection.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MM_OPERA_HPC/Utilities.hxx"
#include "MM_OPERA_HPC/GrainOrientations.hxx"
#include "MM_OPERA_HPC/MacroscropicElasticMaterialProperties.hxx"
#include "MM_OPERA_HPC/UniaxialMacroscopicStressPeriodicSimulation.hxx"

/*!
 * \brief parameters related to the mesh
 */
struct MeshParameters {
  //! \brief default value for the mesh file
  const char *mesh_file = "mesh/5crystals.msh";
  //! \brief default value for the read mesh mode ("Restart")
  const char *mesh_mode = "FromScratch";
  //! \brief number of uniform refinement of the mesh
  int refinement = 0;
};

/*!
 * \brief parameters related to the numerical aspects of the simulation
 */
struct NumericalParameters {
  //! \brief polynomial order of the finite element used
  int order = 1;
  /*!
   * \brief secant acceleration of the fixed point on the lateral macroscopic
   * stresses (see UniaxialMacroscopicStressPeriodicSimulation)
   */
  bool macroscopic_secant_acceleration = false;
  //! \brief name of the linear solver
  const char *linear_solver = "HyprePCG";
  //! \brief name of the preconditioner for the linear solver
  const char *linear_solver_preconditioner = "HypreBoomerAMG";
  //! \brief default tolerance for the convergence of the linear solver
  mfem_mgis::real linear_solver_tolerance = mfem_mgis::real{1e-12};
  /*!
   * \brief non linear solver: `Newton`, or `NewtonK0` (Newton's method, the
   * conjugate gradient being preconditioned by BoomerAMG built once on the
   * elastic operator, see the `FixedPoint` solver of mfem-mgis)
   */
  const char *non_linear_solver = "Newton";
  /*!
   * \brief maximum number of iterations of the non linear solver for each
   * resolution of the microscopic equilibrium: a resolution which does not
   * converge within this number of iterations fails, and the time step is
   * divided by two
   */
  int non_linear_solver_maximum_iterations = 20;
  /*!
   * \brief relative tolerance of the linear solver of the non linear solver
   * adapted to the decrease of the residual: `None`, `EisenstatWalker1` or
   * `EisenstatWalker2` (see `AdaptiveLinearSolverTolerance` in mfem-mgis)
   */
  const char *adaptive_linear_solver_tolerance = "None";
  //! \brief initial adaptive tolerance, the default one if not given
  mfem_mgis::real adaptive_linear_solver_initial_tolerance = -1;
  /*!
   * \brief maximum adaptive tolerance, the default one (0.9) if not given:
   * the tolerance of Eisenstat and Walker rises to it when the residual
   * increases
   */
  mfem_mgis::real adaptive_linear_solver_maximum_tolerance = -1;
  /*!
   * \brief linear solver of the elastic operator of `NewtonK0`:
   * `HypreBoomerAMG` or `AmgX` (one cycle of the multigrid of AmgX, which
   * requires mfem-mgis compiled with AmgX)
   */
  const char* elastic_operator_linear_solver = "HypreBoomerAMG";
  /*!
   * \brief parameters of AmgX completing its default configuration in
   * mfem-mgis, as `name=value` separated by commas, with the names of AmgX
   */
  const char* amgx_parameters = "";
  //! \brief computations of AmgX in single precision
  bool amgx_single_precision = false;
  //! \brief reuse the analysis of the direct solvers (MUMPSSolver,
  //! CuDSSSolver) for the matrices of the same pattern
  bool reuse_analysis = false;
};

/*!
 * \brief parameters describing the elastic material properties
 */
struct ElasticMaterialPropertiesParameters {
  //! \brief Young's modulus of uranium dioxide
  mfem_mgis::real young_modulus = mfem_mgis::real{222.e9};
  //! \brief Poisson's ratio of uranium dioxide
  mfem_mgis::real poisson_ratio = mfem_mgis::real{0.27};
  //! \brief shear modulus of uranium dioxide
  mfem_mgis::real shear_modulus = mfem_mgis::real{54.e9};
};

struct MaterialParameters : ElasticMaterialPropertiesParameters {
  //! \brief default material library
  const char* library = "";
  /*!
   * behaviour library used: the one given by the user (option --library),
   * or else the one built with the case for the processor of the behaviours
   * (option --execution)
   */
  std::string behaviour_library;
  //! \brief default mechanical behaviour
  const char *behaviour = "Mono_UO2_Cosh_Jaco3";
  //! \brief default file name describing orientation vectors
  const char *vect_file = "mesh/vectors_5crystals.txt";
  //! \brief give the displacement gradient to the behaviour
  bool displacement_gradient = false;
  /*!
   * \brief one material for all the grains, whose orientations are given
   * point by point, instead of one material per grain: the behaviour is
   * integrated by one batch over all the grains
   */
  bool single_material = false;
};

/*!
 * \brief grains of the mesh when all the grains are one material (option
 * --single-material)
 */
struct Grains {
  //! \brief number of grains (attributes of the mesh)
  int number = 0;
  //! \brief grain of each element (former attribute, from 1)
  std::vector<int> of_elements;
  //! \brief grain of each integration point of the material (from 0)
  std::vector<int> of_points;
};

struct LoadingParameters {
  //! \brief default value of the temperature
  mfem_mgis::real temperature = mfem_mgis::real{1600};
  //! \brief default value of the imposed linear strain rate
  mfem_mgis::real linear_strain_rate = mfem_mgis::real{5e-4};
};

struct TimeDiscretizationParameters {
  //! \brief total duration of the simulated test
  mfem_mgis::real duration = mfem_mgis::real{200};
  //! \brief number of step describing the temporal sequences
  int nstep = 200;
};

struct PostProcessingParameters {
  //! \brief boolean state if post-processings shall be executed
  bool post_processings = true;
  //! \brief boolean stating if the displacement is exported to ParaView at
  //! each step (only used if the post-processings are executed; the averages
  //! of the stresses per grain are always computed)
  bool paraview_export = true;
  //! \brief boolean stating if the von Mises stress of the Cauchy stress is
  //! exported to VTK
  bool export_von_Mises_stress = false;
  //! \brief boolean stating if first eigen stress of the Cauchy stress is
  //! exported to VTK
  bool export_first_eigen_stress = false;
};

// main parameters
struct TestParameters : MeshParameters,
                        NumericalParameters,
                        MaterialParameters,
                        LoadingParameters,
                        TimeDiscretizationParameters,
                        PostProcessingParameters {
  //! \brief default output file for macroscopic results
  const char *output_file = "uniaxial-polycrystal.res";
  //! \brief default verbositiy level
  int verbosity_level = 0;  // default value : lower level
};

// a few utility functions defined after the main functions

static void parseCommandLineArguments(mfem::OptionsParser &, TestParameters &);
static mfem_mgis::Parameter getElasticOperatorLinearSolver(
    const NumericalParameters&);
static void addPostProcessings(mfem_mgis::attributes::MayAbort,
                               mfem_mgis::Context &,
                               mfem_mgis::PeriodicNonLinearEvolutionProblem &,
                               const std::string &,
                               const bool,
                               const Grains &);
static void addOptionalPostProcessings(
    mfem_mgis::attributes::MayAbort,
    mfem_mgis::Context &,
    mfem_mgis::PeriodicNonLinearEvolutionProblem &,
    const PostProcessingParameters &,
    mfem::ParaViewDataCollection *&);
static void setupMaterials(
    mfem_mgis::attributes::MayAbort,
    mfem_mgis::Context &,
    mfem_mgis::PeriodicNonLinearEvolutionProblem &,
    mm_opera_hpc::MacroscropicElasticMaterialProperties &,
    const TestParameters &,
    Grains &);
static Grains mergeGrains(mfem_mgis::FiniteElementDiscretization &);
static void addGrainsMeanStresses(mfem_mgis::attributes::MayAbort,
                                  mfem_mgis::Context &,
                                  mfem_mgis::PeriodicNonLinearEvolutionProblem &,
                                  const Grains &,
                                  const std::string &);
static void setLinearSolver(mfem_mgis::attributes::MayAbort,
                            mfem_mgis::Context &,
                            mfem_mgis::PeriodicNonLinearEvolutionProblem &,
                            const TestParameters &);

int main(int argc, char *argv[]) {
  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  ctx.enableProfiling(true);
  // mpi initialization here
  mfem_mgis::initialize(argc, argv);

  // get parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  parseCommandLineArguments(args, p);

  // creating the finite element workspace
  auto fed = mfem_mgis::construct<mfem_mgis::FiniteElementDiscretization>(
      ctx, mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
                                 {"MeshReadMode", p.mesh_mode},
                                 {"FiniteElementFamily", "H1"},
                                 {"FiniteElementOrder", p.order},
                                 {"UnknownsSize", mfem_mgis::size_type{3}},
                                 {"NumberOfUniformRefinements", p.refinement},
                                 {"Parallel", true}})|or_die;
  // one material for all the grains: the elements are given the attribute 1
  // before the creation of the problem, their grains are kept
  auto grains = Grains{};
  if (p.single_material) {
    grains = mergeGrains(fed);
  }
  auto problem_parameters = mfem_mgis::Parameters{};
  // options of the non linear solver: adaptive tolerance of its linear
  // solver, if required
  auto non_linear_solver_options = mfem_mgis::Parameters{};
  const auto initial_tolerance_given =
      p.adaptive_linear_solver_initial_tolerance != -1;
  const auto maximum_tolerance_given =
      p.adaptive_linear_solver_maximum_tolerance != -1;
  if ((std::string_view{p.adaptive_linear_solver_tolerance} == "None") &&
      (initial_tolerance_given || maximum_tolerance_given)) {
    mfem_mgis::Profiler::Utils::Message(
        "error: an initial or maximum tolerance is given without adaptive "
        "tolerance");
    std::exit(EXIT_FAILURE);
  }
  // options of the linear solver of the problem
  const auto direct_solver =
      (std::string_view{p.linear_solver} == "MUMPSSolver") ||
      (std::string_view{p.linear_solver} == "CuDSSSolver");
  const auto newton_k0 = std::string_view{p.non_linear_solver} == "NewtonK0";
  if (mfem_mgis::usePETSc()) {
    if ((direct_solver) || (p.reuse_analysis)) {
      mfem_mgis::Profiler::Utils::Message(
          "warning: the options --linear-solver and --reuse-analysis are "
          "ignored (the linear solver is given by the configuration file of "
          "PETSc)");
    }
  } else if (newton_k0) {
    if ((direct_solver) || (p.reuse_analysis)) {
      mfem_mgis::Profiler::Utils::Message(
          "warning: the options --linear-solver and --reuse-analysis are "
          "ignored by NewtonK0 (its tangent systems are solved by a "
          "conjugate gradient preconditioned by the elastic operator)");
    }
  } else if (direct_solver) {
    // the adaptive tolerance requires an iterative solver
    if (std::string_view{p.adaptive_linear_solver_tolerance} != "None") {
      mfem_mgis::Profiler::Utils::Message(
          "error: the adaptive tolerance of the linear solver "
          "(--adaptive-linear-solver-tolerance) requires an iterative "
          "linear solver");
      std::exit(EXIT_FAILURE);
    }
  } else if (p.reuse_analysis) {
    mfem_mgis::Profiler::Utils::Message(
        "warning: the option --reuse-analysis is ignored (the linear solver "
        "is not a direct solver)");
  }
  if (p.non_linear_solver_maximum_iterations < 1) {
    mfem_mgis::Profiler::Utils::Message(
        "error: the maximum number of iterations of the non linear solver "
        "(--non-linear-solver-maximum-iterations) must be at least 1");
    std::exit(EXIT_FAILURE);
  }
  if (std::string_view{p.adaptive_linear_solver_tolerance} != "None") {
    auto adaptive_options = mfem_mgis::Parameters{};
    if (initial_tolerance_given) {
      adaptive_options.insert(ctx, "InitialTolerance",
                              p.adaptive_linear_solver_initial_tolerance) |
          or_die;
    }
    if (maximum_tolerance_given) {
      adaptive_options.insert(ctx, "MaximumTolerance",
                              p.adaptive_linear_solver_maximum_tolerance) |
          or_die;
    }
    non_linear_solver_options.insert(
        ctx, "AdaptiveLinearSolverTolerance",
        mfem_mgis::Parameters{
            {p.adaptive_linear_solver_tolerance, adaptive_options}}) |
        or_die;
  }
  if (std::string_view{p.non_linear_solver} == "NewtonK0") {
    // the behaviour does not provide its prediction operator: the elastic
    // operator is given by an integration
    const auto tangent_linear_solver_options = mfem_mgis::Parameters{
        {"RelativeTolerance", p.linear_solver_tolerance},
        {"MaximumNumberOfIterations", 2000}};
    auto options = mfem_mgis::Parameters{
        {"Method", "NewtonK0"},
        {"ElasticOperatorComputation", "Integration"},
        {"ElasticOperatorUpdate", "Once"},
        {"ElasticOperatorLinearSolver", getElasticOperatorLinearSolver(p)},
        {"TangentLinearSolver",
         mfem_mgis::Parameters{{"CGSolver", tangent_linear_solver_options}}}};
    for (const auto &[k, v] : non_linear_solver_options) {
      options.insert(ctx, k, v) | or_die;
    }
    problem_parameters.insert(ctx, "NonLinearSolver",
                              mfem_mgis::Parameters{{"FixedPoint", options}}) |
        or_die;
  } else if (std::string_view{p.non_linear_solver} == "Newton") {
    if (!non_linear_solver_options.empty()) {
      problem_parameters.insert(
          ctx, "NonLinearSolver",
          mfem_mgis::Parameters{{"Newton", non_linear_solver_options}}) |
          or_die;
    }
  } else {
    mfem_mgis::Profiler::Utils::Message(
        "error: invalid non linear solver (expected Newton or NewtonK0)");
    std::exit(EXIT_FAILURE);
  }
  auto problem =
      mfem_mgis::construct<mfem_mgis::PeriodicNonLinearEvolutionProblem>(
          ctx, fed, mfem_mgis::FIX_XMIN, problem_parameters) |
      or_die;
  mm_opera_hpc::printMeshInformation(ctx, problem);

  // set problem
  mm_opera_hpc::MacroscropicElasticMaterialProperties mp;
  setupMaterials(mfem_mgis::may_abort, ctx, problem, mp, p, grains);
  setLinearSolver(mfem_mgis::may_abort, ctx, problem, p);
  // add post processings
  mfem::ParaViewDataCollection *paraview_exporter = nullptr;
  if (p.post_processings) {
    addPostProcessings(mfem_mgis::may_abort, ctx, problem,
                       "OutputFile-Uniaxial-polycrystal", p.paraview_export,
                       grains);
    addOptionalPostProcessings(mfem_mgis::may_abort, ctx, problem, p,
                               paraview_exporter);
  }
  // definition of the temporal sequences
  const auto te = p.duration;
  const auto nsteps = p.nstep;
  const auto temporal_sequences = [&te, &nsteps] {
    auto times = std::vector<mfem_mgis::real>{};
    const double dt = te / nsteps;
    times.reserve(nsteps + 1);
    for (std::size_t i = 0; i != nsteps + 1; ++i) {
      times.push_back(i * dt);
    }
    return times;
  }();
  // evolution of the axial component of the displacement gradient
  const auto Hzz = [&p](const mfem_mgis::real ets) {
    return p.linear_strain_rate * ets;
  };
  // main output file, written by the main process only
  std::ofstream out;
  if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) {
    out.open(p.output_file);
    if (!out) {
      mfem_mgis::abort(
          ("can't open output file '" + std::string{p.output_file} + "'")
              .c_str());
    }
    out.precision(14);
  }
  // setting the simulation and running it
  const auto np = mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation::
      NumericalParameters{
          .microscopic_equilibrium_maximum_number_of_iterations =
              p.non_linear_solver_maximum_iterations,
          .macroscopic_elastic_material_properties = mp,
          .macroscopic_stress_secant_acceleration =
              p.macroscopic_secant_acceleration};
  mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation s(
      problem, Hzz, np, p.post_processings, p.displacement_gradient);
  const auto success = s.run(ctx, out, temporal_sequences);

  if (paraview_exporter != nullptr) {
    delete paraview_exporter;
  }
  // print and write timetable
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

static void parseCommandLineArguments(mfem::OptionsParser &args,
                                      TestParameters &p) {
  // options of mfem-mgis (device of MFEM...)
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&p.output_file, "", "--macroscopic-stress-output-file",
                 "main output file containing the evolution of the diagonal "
                 "components of the deformation gradient and the  diagonal "
                 "components of the Cauchy stress.");
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.mesh_mode, "-mm", "--mesh_mode",
                 "Read Mesh Mode (FromScratch or Restart = splited file).");
  args.AddOption(&p.vect_file, "-f", "--vect", "Vector file to use.");
  args.AddOption(&p.library, "-l", "--library",
                 "Material library. By default, the library built with the "
                 "case for the CPU or for the GPU, as the option --execution "
                 "places the behaviours.");
  args.AddOption(&p.behaviour, "-b", "--behaviour", "Mechanical behaviour.");
  args.AddOption(&p.linear_solver, "", "--linear-solver",
                 "linear solver of the Newton method: a Krylov solver of "
                 "hypre (HyprePCG by default, HypreGMRES or HypreFGMRES, "
                 "preconditioned by --linear-solver-preconditioner), or a "
                 "direct solver of mfem-mgis, MUMPSSolver or CuDSSSolver "
                 "(on a CUDA GPU only), the matrices being given as general "
                 "(non symmetric) ones; ignored by NewtonK0");
  args.AddOption(&p.reuse_analysis, "", "--reuse-analysis", "",
                 "--no-reuse-analysis",
                 "reuse the analysis of the direct solver (ordering and "
                 "symbolic factorization) for the matrices of the same "
                 "pattern, as the tangent operators of the successive Newton "
                 "iterations (parameter ReuseAnalysis of mfem-mgis); "
                 "ignored, with a warning, without a direct solver");
  args.AddOption(&p.linear_solver_preconditioner, "",
                 "--linear-solver-preconditioner",
                 "preconditioner of the linear solver to be used");
  args.AddOption(&p.non_linear_solver, "", "--non-linear-solver",
                 "non linear solver: Newton (default) or NewtonK0");
  args.AddOption(&p.non_linear_solver_maximum_iterations, "",
                 "--non-linear-solver-maximum-iterations",
                 "maximum number of iterations of the non linear solver for "
                 "each resolution of the microscopic equilibrium (a "
                 "resolution which does not converge fails, and the time "
                 "step is divided by two)");
  args.AddOption(&p.adaptive_linear_solver_tolerance, "",
                 "--adaptive-linear-solver-tolerance",
                 "relative tolerance of the linear solver adapted to the "
                 "decrease of the residual: None (default), "
                 "EisenstatWalker1 or EisenstatWalker2");
  args.AddOption(&p.adaptive_linear_solver_initial_tolerance, "",
                 "--adaptive-linear-solver-initial-tolerance",
                 "initial adaptive tolerance, the default one if -1");
  args.AddOption(&p.adaptive_linear_solver_maximum_tolerance, "",
                 "--adaptive-linear-solver-maximum-tolerance",
                 "maximum adaptive tolerance, the default one (0.9) if -1");
  args.AddOption(&p.elastic_operator_linear_solver, "",
                 "--elastic-operator-linear-solver",
                 "linear solver of the elastic operator of NewtonK0: "
                 "HypreBoomerAMG (default) or AmgX (one cycle of the "
                 "multigrid of AmgX)");
  args.AddOption(&p.amgx_parameters, "", "--amgx-parameters",
                 "parameters of AmgX, as name=value separated by commas, "
                 "for instance strength_threshold=0.8,interpolator=MULTIPASS");
  args.AddOption(&p.amgx_single_precision, "", "--amgx-single-precision", "",
                 "--no-amgx-single-precision",
                 "computations of AmgX in single precision");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "refinement level of the mesh, default = 0");
  args.AddOption(&p.export_von_Mises_stress, "",
                 "--enable-export-von-Mises-stress", "",
                 "--disable-export-von-Mises-stress",
                 "export the von Mises stress", false);
  args.AddOption(&p.export_first_eigen_stress, "",
                 "--enable-export-first_eigen_stress", "",
                 "--disable-export-first_eigen_stress",
                 "export first eigen stress", false);
  args.AddOption(&p.post_processings, "", "--enable-post-processing", "",
                 "--disable-post-processing", "run post processing steps",
                 false);
  args.AddOption(&p.paraview_export, "", "--enable-paraview-export", "",
                 "--disable-paraview-export",
                 "export the displacement to ParaView at each step (the "
                 "averages of the stresses per grain are still computed)",
                 false);
  args.AddOption(&p.single_material, "", "--single-material", "",
                 "--material-per-grain",
                 "one material for all the grains, oriented point by point "
                 "(one batch of integration over all the grains), instead of "
                 "one material per grain");
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "choose the verbosity level");
  args.AddOption(&p.duration, "-d", "--duration",
                 "choose the duration");
  args.AddOption(&p.nstep, "-n", "--nstep",
                 "choose the number of steps");
  args.AddOption(&p.macroscopic_secant_acceleration, "",
                 "--macroscopic-secant-acceleration", "",
                 "--no-macroscopic-secant-acceleration",
                 "Secant acceleration (Broyden) of the fixed point on the "
                 "lateral macroscopic stresses, whose correction otherwise "
                 "uses the macroscopic elastic compliance.");
  args.AddOption(&p.displacement_gradient, "-dg", "--displacement-gradient",
                 "-no-dg", "--no-displacement-gradient",
                 "Give the displacement gradient instead of the deformation "
                 "gradient to the behaviour (parameter DisplacementGradient "
                 "of mfem-mgis). In single precision, the macroscopic "
                 "loading keeps its precision; the behaviour, written with "
                 "the deformation gradient, receives F = I + H.");

  args.Parse();

  if (!args.Good()) {
    if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) args.PrintUsage(std::cout);
    mfem_mgis::finalize();
    // the help (-h) is a success, an invalid option a failure, which a
    // test must see
    exit(args.Help() ? EXIT_SUCCESS : EXIT_FAILURE);
  }
  if (p.mesh_file == nullptr) {
    if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) {
      std::cout << "ERROR: Mesh file missing" << std::endl;
      args.PrintUsage(std::cout);
    }
    mfem_mgis::abort(EXIT_FAILURE);
  }
  if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) args.PrintOptions(std::cout);
  // behaviour library: given by --library, or else the library built with
  // the case for the processor of the behaviours (option --execution)
  p.behaviour_library = std::string_view{p.library}.empty()
                            ? mfem_mgis::getBehavioursLibrary(
                                  MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY)
                            : std::string{p.library};
}

static void setupMaterials(
    mfem_mgis::attributes::MayAbort,
    mfem_mgis::Context &ctx,
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
    mm_opera_hpc::MacroscropicElasticMaterialProperties &mp,
    const TestParameters &p,
    Grains &grains) {
  using namespace mgis::behaviour;
  auto or_die = ctx.getFatalFailureHandler();
  using real = mfem_mgis::real;

  CatchTimeSection(ctx, "set_mgis_stuff");

  // const int nMat = 8;
  const int nMat =
      getMaterialsAttributes(problem.getFiniteElementDiscretization())
          .Max();
  mfem_mgis::Profiler::Utils::Message("Nombre de matériaux : ", nMat);

  for (int i = 0; i < nMat; i++) {
    problem.addBehaviourIntegrator(
        ctx, "Mechanics", i + 1, p.behaviour_library, p.behaviour,
        {{mfem_mgis::displacement_gradient_parameter,
          p.displacement_gradient}}) |
        or_die;
  }

  // cubic symmetry elasticity
  const mfem_mgis::real young1 = p.young_modulus;
  const mfem_mgis::real young2 = young1;
  const mfem_mgis::real young3 = young1;
  const mfem_mgis::real poisson12 = p.poisson_ratio;
  const mfem_mgis::real poisson23 = poisson12;
  const mfem_mgis::real poisson13 = poisson12;
  const mfem_mgis::real shear12 = p.shear_modulus;
  const mfem_mgis::real shear23 = shear12;
  const mfem_mgis::real shear13 = shear12;

  //
  mp.update(young1, poisson12, shear12);
  // materials
  auto set_properties =
      [&ctx, &or_die](auto &m, const mfem_mgis::real yo1,
                      const mfem_mgis::real yo2, const mfem_mgis::real yo3,
                      const mfem_mgis::real po12, const mfem_mgis::real po23,
                      const mfem_mgis::real po13, const mfem_mgis::real sm12,
                      const mfem_mgis::real sm23, const mfem_mgis::real sm13) {
        for (auto s : {&m.s0, &m.s1}) {
          setMaterialProperty(ctx, *s, "YoungModulus1", yo1) | or_die;
          setMaterialProperty(ctx, *s, "YoungModulus2", yo2) | or_die;
          setMaterialProperty(ctx, *s, "YoungModulus3", yo3) | or_die;
          setMaterialProperty(ctx, *s, "PoissonRatio12", po12) | or_die;
          setMaterialProperty(ctx, *s, "PoissonRatio23", po23) | or_die;
          setMaterialProperty(ctx, *s, "PoissonRatio13", po13) | or_die;
          setMaterialProperty(ctx, *s, "ShearModulus12", sm12) | or_die;
          setMaterialProperty(ctx, *s, "ShearModulus23", sm23) | or_die;
          setMaterialProperty(ctx, *s, "ShearModulus13", sm13) | or_die;
        }
      };

  auto set_temperature = [&ctx, &or_die, &p](auto &m) {
    setExternalStateVariable(ctx, m.s0, "Temperature", p.temperature) | or_die;
    setExternalStateVariable(ctx, m.s1, "Temperature", p.temperature) | or_die;
  };

  for (int i = 0; i < nMat; i++) {
    auto &mat = problem.getMaterial(ctx, i + 1, 0) | or_die;
    set_properties(mat, young1, young2, young3,      // young modulus
                   poisson12, poisson23, poisson13,  // poisson ration
                   shear12, shear23, shear13         // shear modulus
    );
    set_temperature(mat);
  }

  //
  const auto vectors = mm_opera_hpc::readVectorsFromFile(p.vect_file);
  const auto nGrains = p.single_material ? grains.number : nMat;
  if (vectors.size() != 2 * static_cast<std::size_t>(nGrains)) {
    throw std::invalid_argument(
        "setupMaterials : incorrect number of vectors in vector file");
  }
  if (p.single_material) {
    // orientations of the grains given point by point, and grain of each
    // point
    auto &mat = problem.getMaterial(ctx, 1, 0) | or_die;
    if (mat.b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      return;
    }
    const auto &qspace = mat.getPartialQuadratureSpace();
    auto a1 = std::make_shared<mfem_mgis::PartialQuadratureFunction>(qspace, 3);
    auto a2 = std::make_shared<mfem_mgis::PartialQuadratureFunction>(qspace, 3);
    auto v1 = a1->getValues();
    auto v2 = a2->getValues();
    const auto &elements = qspace.getElements();
    const auto &offsets = qspace.getElementOffsets();
    const auto *const pe = elements.HostRead();
    const auto *const po = offsets.HostRead();
    grains.of_points.assign(
        static_cast<std::size_t>(qspace.getNumberOfIntegrationPoints()), 0);
    for (int i = 0; i != elements.Size(); ++i) {
      const auto g = grains.of_elements[static_cast<std::size_t>(pe[i])] - 1;
      const auto &r1 = vectors[2 * static_cast<std::size_t>(g)];
      const auto &r2 = vectors[2 * static_cast<std::size_t>(g) + 1];
      for (auto o = static_cast<std::size_t>(po[i]);
           o != static_cast<std::size_t>(po[i + 1]); ++o) {
        grains.of_points[o] = g;
        for (std::size_t k = 0; k != 3; ++k) {
          v1[3 * o + k] = r1[k];
          v2[3 * o + k] = r2[k];
        }
      }
    }
    mat.setRotationMatrix(mfem_mgis::RotationMatrix3D{
        std::array<mfem_mgis::MaterialAxis3D, 2u>{a1, a2}});
    return;
  }

  std::array<mfem_mgis::MaterialAxis3D, 2u> r;
  for (int i = 0; i < nMat; i++) {
    auto &mat = problem.getMaterial(ctx, i + 1, 0) | or_die;
    if (mat.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {
      r[0] = vectors[2 * i];
      r[1] = vectors[2 * i + 1];
      mat.setRotationMatrix(mfem_mgis::RotationMatrix3D{r});
    }
  }
}

/*!
 * \return parameters given as `name=value` separated by commas: integers,
 * reals or strings
 * \param[in] o: parameters
 */
static mfem_mgis::Parameters parseNameValueParameters(const std::string& o) {
  auto params = mfem_mgis::Parameters{};
  auto start = std::string::size_type{};
  while (start < o.size()) {
    const auto end = std::min(o.find(',', start), o.size());
    const auto option = o.substr(start, end - start);
    const auto equal = option.find('=');
    if ((equal == std::string::npos) || (equal == 0) ||
        (equal + 1 == option.size())) {
      mfem_mgis::Profiler::Utils::Message("error: invalid parameter '" +
                                          option + "'");
      std::exit(EXIT_FAILURE);
    }
    const auto name = option.substr(0, equal);
    const auto value = option.substr(equal + 1);
    auto n = std::size_t{};
    try {
      const auto i = std::stoi(value, &n);
      if (n == value.size()) {
        params.replaceOrInsert(name, i);
        start = end + 1;
        continue;
      }
      const auto r = std::stod(value, &n);
      if (n == value.size()) {
        params.replaceOrInsert(name, r);
        start = end + 1;
        continue;
      }
    } catch (std::exception&) {
    }
    params.replaceOrInsert(name, value);
    start = end + 1;
  }
  return params;
}  // end of parseNameValueParameters

/*!
 * \return the linear solver of the elastic operator of `NewtonK0`
 * \param[in] p: parameters
 */
static mfem_mgis::Parameter getElasticOperatorLinearSolver(
    const NumericalParameters& p) {
  const auto s = std::string_view{p.elastic_operator_linear_solver};
  if (s == "HypreBoomerAMG") {
    return std::string{s};
  }
  if (s != "AmgX") {
    mfem_mgis::Profiler::Utils::Message(
        "error: invalid linear solver of the elastic operator (expected "
        "HypreBoomerAMG or AmgX)");
    std::exit(EXIT_FAILURE);
  }
  auto options = mfem_mgis::Parameters{};
  if (!std::string_view{p.amgx_parameters}.empty()) {
    options.replaceOrInsert("AmgXParameters",
                            parseNameValueParameters(p.amgx_parameters));
  }
  if (p.amgx_single_precision) {
    options.replaceOrInsert("Precision", std::string{"Single"});
  }
  return mfem_mgis::Parameters{{"AmgX", options}};
}  // end of getElasticOperatorLinearSolver

static void setLinearSolver(mfem_mgis::attributes::MayAbort,
                            mfem_mgis::Context &ctx,
                            mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
                            const TestParameters &params) {
  CatchTimeSection(ctx, "set_linear_solver");
  auto or_die = ctx.getFatalFailureHandler();
  // direct solvers: no tolerance nor preconditioner, ReuseAnalysis only given
  // when requested (for the versions of mfem-mgis without this parameter)
  const auto solver = std::string_view{params.linear_solver};
  if ((solver == "MUMPSSolver") || (solver == "CuDSSSolver")) {
    auto direct_solver_parameters = mfem_mgis::Parameters{};
    if (params.reuse_analysis) {
      direct_solver_parameters.insert(
          mfem_mgis::may_throw, mfem_mgis::Parameters{{"ReuseAnalysis", true}});
    }
    p.setLinearSolver(ctx, params.linear_solver, direct_solver_parameters) |
        or_die;
    return;
  }
  // pilote
  constexpr int defaultMaxNumOfIt = 5000;  // MaximumNumberOfIterations
  auto solverParameters = mfem_mgis::Parameters{};
  solverParameters.insert(
      mfem_mgis::may_throw,
      mfem_mgis::Parameters{{"VerbosityLevel", params.verbosity_level}});
  solverParameters.insert(
      mfem_mgis::may_throw,
      mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
  // solverParameters.insert(mfem_mgis::Parameters{{"AbsoluteTolerance", Tol}});
  // solverParameters.insert(mfem_mgis::Parameters{{"RelativeTolerance", Tol}});
  solverParameters.insert(
      mfem_mgis::may_throw,
      mfem_mgis::Parameters{{"Tolerance", params.linear_solver_tolerance}});

  // preconditioner
  auto options =
      mfem_mgis::Parameters{{"VerbosityLevel", params.verbosity_level}};
  // auto preconditioner = mfem_mgis::Parameters{{"Name","HypreDiagScale"},
  // {"Options",options}};
  if (std::string_view{params.linear_solver_preconditioner} != "none") {
    auto preconditioner = mfem_mgis::Parameters{
        {"Name", params.linear_solver_preconditioner}, {"Options", options}};
    solverParameters.insert(
        mfem_mgis::may_throw,
        mfem_mgis::Parameters{{"Preconditioner", preconditioner}});
  }
  p.setLinearSolver(ctx, params.linear_solver, solverParameters) | or_die;
}

static void addOptionalPostProcessings(
    mfem_mgis::attributes::MayAbort,
    mfem_mgis::Context &ctx,
    mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
    const PostProcessingParameters &params,
    mfem::ParaViewDataCollection *&exporter) {
  auto or_die = ctx.getFatalFailureHandler();
  // setup paraview
  exporter = new mfem::ParaViewDataCollection("PolycrystalParaviewOutput");
  exporter->SetDataFormat(mfem::VTKFormat::BINARY);
  exporter->SetMesh(&(p.getFiniteElementDiscretization().getMesh<true>()));

  if (params.export_von_Mises_stress || params.export_first_eigen_stress) {
    p.getImplementation<true>().addPostProcessing(ctx, [&ctx, &p, &params,
                                                        &exporter, &or_die](
                                                           mfem_mgis::real t,
                                                           mfem_mgis::real dt) {
      static int count = 0;
      exporter->SetCycle(count++);
      // the fields are the ones of the end of the time step
      exporter->SetTime(t + dt);
      const int nMat =
          getMaterialsAttributes(p.getFiniteElementDiscretization())
              .Max();
      // fes = FirstEigenStress
      // vmes = VonMisesEquivalentStress
      std::vector<mfem_mgis::ImmutablePartialQuadratureFunctionView> fesView{};
      std::vector<mfem_mgis::ImmutablePartialQuadratureFunctionView> vmesView{};
      // the functions must outlive the views used by the projections
      std::vector<mfem_mgis::PartialQuadratureFunction> fes_pqfs;
      std::vector<mfem_mgis::PartialQuadratureFunction> vmes_pqfs;
      fes_pqfs.reserve(nMat);
      vmes_pqfs.reserve(nMat);
      auto ets = mfem_mgis::Material::END_OF_TIME_STEP;
      // For each material
      for (int i = 1; i <= nMat; i++) {
        auto &m = p.getMaterial(ctx, i, 0) | or_die;
        if (params.export_first_eigen_stress) {
          auto &fes_pqf =
              fes_pqfs.emplace_back(m.getPartialQuadratureSpace(), 1);
          mfem_mgis::computeFirstEigenStress(ctx, fes_pqf, m, ets);
          fesView.push_back(fes_pqf);
        }
        if (params.export_von_Mises_stress) {
          auto &vmes_pqf =
              vmes_pqfs.emplace_back(m.getPartialQuadratureSpace(), 1);
          bool success =
              mfem_mgis::computeVonMisesEquivalentStress(ctx, vmes_pqf, m, ets);
          if (!success) {
            std::cout << "ERROR: computeVonMisesEquivalentStress returns"
                      << " false in addOptionalPostProcessings" << std::endl;
          }
          vmesView.push_back(vmes_pqf);
        }
      }

      // Used by computeL2Projection
      auto &lsf = mfem_mgis::LinearSolverFactory<true>::getFactory();
      auto &fespace =
          p.getFiniteElementDiscretization().getFiniteElementSpace<true>();
      auto linear_solver = lsf.generate(
          ctx, "HyprePCG", fespace,
          mfem_mgis::Parameters{
              {"Preconditioner", mfem_mgis::Parameters{{"Name", "HypreILU"}}}});

      // project solutions
      std::optional<mfem_mgis::L2ProjectionResult<true>> fes_projection;
      std::optional<mfem_mgis::L2ProjectionResult<true>> vmes_projection;
      if (params.export_first_eigen_stress) {
        fes_projection =
            mfem_mgis::computeL2Projection<true>(ctx, linear_solver, fesView);
        exporter->RegisterField("FirstEigenStress",
                                fes_projection->result.get());
      }
      if (params.export_von_Mises_stress) {
        vmes_projection =
            mfem_mgis::computeL2Projection<true>(ctx, linear_solver, vmesView);
        exporter->RegisterField("vonMisesStress",
                                vmes_projection->result.get());
      }

      // save paraview files
      exporter->Save();
    }) | or_die;
  }
}

static void addPostProcessings(mfem_mgis::attributes::MayAbort,
                               mfem_mgis::Context &ctx,
                               mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
                               const std::string &msg,
                               const bool paraview_export,
                               const Grains &grains) {
  auto or_die = ctx.getFatalFailureHandler();
  if (paraview_export) {
    p.addPostProcessing(ctx, "ParaviewExportResults",
                        {{"OutputFileName", msg}}) |
        or_die;
  }
  if (grains.number != 0) {
    addGrainsMeanStresses(mfem_mgis::may_abort, ctx, p, grains,
                          "avgStressPolycrystal");
    return;
  }
  p.addPostProcessing(ctx, "MeanThermodynamicForces",
                      {{"OutputFileName", "avgStressPolycrystal"}}) |
      or_die;
}  // end timer add_postprocessing_and_outputs

/*!
 * \brief gives all the elements the attribute 1: one material for all the
 * grains (option --single-material)
 * \return the grains, the grain of each element being its former attribute
 * \param[in, out] fed: finite element discretization, before the creation
 * of the problem
 */
static Grains mergeGrains(mfem_mgis::FiniteElementDiscretization &fed) {
  auto grains = Grains{};
  // attributes of all the processes
  grains.number = static_cast<int>(getMaterialsAttributes(fed).Max());
  auto &mesh = fed.getMesh<true>();
  grains.of_elements.resize(static_cast<std::size_t>(mesh.GetNE()));
  for (int e = 0; e != mesh.GetNE(); ++e) {
    grains.of_elements[static_cast<std::size_t>(e)] = mesh.GetAttribute(e);
    mesh.SetAttribute(e, 1);
  }
  mesh.SetAttributes();
  mfem_mgis::Profiler::Utils::Message("Nombre de grains : ", grains.number);
  return grains;
}  // end of mergeGrains

/*!
 * \brief averages of the stresses (in the global frame) of each grain when
 * all the grains are one material, written at each step as by the
 * post-processing MeanThermodynamicForces with one material per grain (same
 * file, same sums in the order of the elements of each grain)
 * \param[in, out] ctx: execution context
 * \param[in, out] p: problem
 * \param[in] grains: grains of the points of the material
 * \param[in] f: output file
 */
static void addGrainsMeanStresses(mfem_mgis::attributes::MayAbort,
                                  mfem_mgis::Context &ctx,
                                  mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
                                  const Grains &grains,
                                  const std::string &f) {
  auto or_die = ctx.getFatalFailureHandler();
  auto out = std::make_shared<std::ofstream>();
  const auto &m = p.getMaterial(ctx, 1, 0) | or_die;
  const auto nc =
      static_cast<std::size_t>(m.s1.thermodynamic_forces_stride);
  if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) {
    out->open(f);
    if (!(*out)) {
      mfem_mgis::abort(("can't open output file '" + f + "'").c_str());
    }
    *out << "# first column: time\n";
    auto c = std::size_t{2};
    for (int g = 1; g <= grains.number; ++g) {
      for (std::size_t k = 0; k != nc; ++k, ++c) {
        *out << "# " << c << " column: " << k
             << " component of the mean thermodynamic force for material "
             << g << '\n';
      }
    }
  }
  p.getImplementation<true>().addPostProcessing(
      ctx, [&ctx, &p, &grains, out, nc, or_die](const mfem_mgis::real t,
                                                const mfem_mgis::real dt) {
        auto &impl = p.getImplementation<true>();
        const auto &mat = p.getMaterial(ctx, 1, 0) | or_die;
        const auto &w = impl.getIntegrationPointsWeights(ctx, 1, 0) | or_die;
        const auto ng = static_cast<std::size_t>(grains.number);
        // integrals of the stresses and volume of each grain
        auto sums = std::vector<mfem_mgis::real>(ng * (nc + 1), 0);
        auto add = [&](const mgis::FieldView<const mfem_mgis::real> &values) {
          for (std::size_t o = 0; o != w.size(); ++o) {
            const auto g = static_cast<std::size_t>(grains.of_points[o]);
            auto *const s = sums.data() + g * (nc + 1);
            for (std::size_t k = 0; k != nc; ++k) {
              s[k] += w[o] * values(static_cast<mfem_mgis::size_type>(o),
                                    static_cast<mfem_mgis::size_type>(k));
            }
            s[nc] += w[o];
          }
        };
        if (mat.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {
          auto rotated = mfem::Vector{};
          mat.computeThermodynamicForcesInGlobalFrame(rotated);
          add(mgis::makePointByPointFieldView(
              std::span<const mfem_mgis::real>(
                  rotated.HostRead(), static_cast<std::size_t>(rotated.Size())),
              static_cast<mfem_mgis::size_type>(w.size()),
              static_cast<mfem_mgis::size_type>(nc),
              static_cast<mfem_mgis::size_type>(nc)));
        } else {
          mat.synchronizeOnHostForReading(
              mfem_mgis::MaterialArrays::StateArray::THERMODYNAMIC_FORCES,
              mfem_mgis::TimeStepStage::END_OF_TIME_STEP);
          add(mgis::behaviour::getThermodynamicForcesView(mat.s1));
        }
        auto totals = std::vector<mfem_mgis::real>(sums.size(), 0);
        MPI_Reduce(sums.data(), totals.data(), static_cast<int>(sums.size()),
                   mfem_mgis::mpi_type<mfem_mgis::real>, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) {
          *out << t + dt;
          for (std::size_t g = 0; g != ng; ++g) {
            const auto *const s = totals.data() + g * (nc + 1);
            for (std::size_t k = 0; k != nc; ++k) {
              *out << " " << s[k] / s[nc];
            }
          }
          *out << std::endl;
        }
      }) |
      or_die;
}  // end of addGrainsMeanStresses
