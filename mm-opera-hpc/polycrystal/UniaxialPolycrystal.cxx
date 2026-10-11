/*!
 * \file   Mono_UO2_CosH_Jaco.cxx
 * \brief
 * This example is modelling ...
 * \author Raphael Prat, Maxence Wangermez
 * \date   06/2023
 */

#include <memory>
#include <string>
#include <algorithm>
#include <exception>
#include <csignal>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <functional>
#include <string_view>

#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
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
   * \brief relative tolerance of the linear solver of the non linear solver
   * adapted to the decrease of the residual: `None`, `EisenstatWalker1` or
   * `EisenstatWalker2` (see `AdaptiveLinearSolverTolerance` in mfem-mgis)
   */
  const char *adaptive_linear_solver_tolerance = "None";
  //! \brief initial adaptive tolerance, the default one if not given
  mfem_mgis::real adaptive_linear_solver_initial_tolerance = -1;
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
  const char *library = "src/libBehaviour.so";
  //! \brief default mechanical behaviour
  const char *behaviour = "Mono_UO2_Cosh_Jaco3";
  //! \brief default file name describing orientation vectors
  const char *vect_file = "mesh/vectors_5crystals.txt";
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
                               const std::string &);
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
    const TestParameters &);
static void setLinearSolver(mfem_mgis::attributes::MayAbort,
                            mfem_mgis::Context &,
                            mfem_mgis::PeriodicNonLinearEvolutionProblem &,
                            const TestParameters &);

int main(int argc, char *argv[]) {
  auto ctx = mfem_mgis::Context{};
  ctx.enableProfiling(true);
  // mpi initialization here
  mfem_mgis::initialize(argc, argv);

  // get parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  parseCommandLineArguments(args, p);

  // creating the finite element workspace
  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      ctx, mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
                                 {"MeshReadMode", p.mesh_mode},
                                 {"FiniteElementFamily", "H1"},
                                 {"FiniteElementOrder", p.order},
                                 {"UnknownsSize", mfem_mgis::size_type{3}},
                                 {"NumberOfUniformRefinements", p.refinement},
                                 {"Parallel", true}});
  auto or_die = ctx.getFatalFailureHandler();
  auto problem_parameters = mfem_mgis::Parameters{};
  // options of the non linear solver: adaptive tolerance of its linear
  // solver, if required
  auto non_linear_solver_options = mfem_mgis::Parameters{};
  const auto initial_tolerance_given =
      p.adaptive_linear_solver_initial_tolerance != -1;
  if ((std::string_view{p.adaptive_linear_solver_tolerance} == "None") &&
      (initial_tolerance_given)) {
    mfem_mgis::Profiler::Utils::Message(
        "error: an initial tolerance is given without adaptive tolerance");
    std::exit(EXIT_FAILURE);
  }
  if (std::string_view{p.adaptive_linear_solver_tolerance} != "None") {
    auto adaptive_options = mfem_mgis::Parameters{};
    if (initial_tolerance_given) {
      adaptive_options.insert(ctx, "InitialTolerance",
                              p.adaptive_linear_solver_initial_tolerance) |
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
  // the constructor of the problem throws on invalid options (of the non
  // linear solver for instance)
  auto problem_ptr = [&ctx, &fed, &problem_parameters] {
    try {
      return std::make_unique<mfem_mgis::PeriodicNonLinearEvolutionProblem>(
          ctx, fed, mfem_mgis::FIX_XMIN, problem_parameters);
    } catch (std::exception &e) {
      mfem_mgis::Profiler::Utils::Message(std::string{"error: "} + e.what());
      std::exit(EXIT_FAILURE);
    }
  }();
  auto &problem = *problem_ptr;
  mm_opera_hpc::printMeshInformation(ctx, problem);

  // set problem
  mm_opera_hpc::MacroscropicElasticMaterialProperties mp;
  setupMaterials(mfem_mgis::may_abort, ctx, problem, mp, p);
  setLinearSolver(mfem_mgis::may_abort, ctx, problem, p);
  // add post processings
  mfem::ParaViewDataCollection *paraview_exporter = nullptr;
  if (p.post_processings) {
    addPostProcessings(mfem_mgis::may_abort, ctx, problem,
                       "OutputFile-Uniaxial-polycrystal");
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
  // evolution of the axial component of the deformation gradient
  const auto Fzz = [&p](const mfem_mgis::real ets) {
    return 1 + p.linear_strain_rate * ets;
  };
  // main output file
  std::ofstream out(p.output_file);
  if (!out) {
    std::cerr << "can't open output file '" << p.output_file << "'\n";
  }
  out.precision(14);
  // setting the simulation and running it
  const auto np = mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation::
      NumericalParameters{.macroscopic_elastic_material_properties = mp};
  mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation s(
      problem, Fzz, np, p.post_processings);
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
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.behaviour, "-b", "--behaviour", "Mechanical behaviour.");
  args.AddOption(&p.linear_solver, "", "--linear-solver",
                 "linear solver to be used");
  args.AddOption(&p.linear_solver_preconditioner, "",
                 "--linear-solver-preconditioner",
                 "preconditioner of the linear solver to be used");
  args.AddOption(&p.non_linear_solver, "", "--non-linear-solver",
                 "non linear solver: Newton (default) or NewtonK0");
  args.AddOption(&p.adaptive_linear_solver_tolerance, "",
                 "--adaptive-linear-solver-tolerance",
                 "relative tolerance of the linear solver adapted to the "
                 "decrease of the residual: None (default), "
                 "EisenstatWalker1 or EisenstatWalker2");
  args.AddOption(&p.adaptive_linear_solver_initial_tolerance, "",
                 "--adaptive-linear-solver-initial-tolerance",
                 "initial adaptive tolerance, the default one if negative");
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
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "choose the verbosity level");
  args.AddOption(&p.duration, "-d", "--duration",
                 "choose the duration (default = 5)");
  args.AddOption(&p.nstep, "-n", "--nstep",
                 "choose the number of steps (default = 40)");

  args.Parse();

  if (!args.Good()) {
    if (mfem_mgis::getMPIrank() == 0) args.PrintUsage(std::cout);
    mfem_mgis::finalize();
    exit(0);
  }
  if (p.mesh_file == nullptr) {
    if (mfem_mgis::getMPIrank() == 0)
      std::cout << "ERROR: Mesh file missing" << std::endl;
    args.PrintUsage(std::cout);
    mfem_mgis::abort(EXIT_FAILURE);
  }
  if (mfem_mgis::getMPIrank() == 0) args.PrintOptions(std::cout);
}

static void setupMaterials(
    mfem_mgis::attributes::MayAbort,
    mfem_mgis::Context &ctx,
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
    mm_opera_hpc::MacroscropicElasticMaterialProperties &mp,
    const TestParameters &p) {
  using namespace mgis::behaviour;
  auto or_die = ctx.getFatalFailureHandler();
  using real = mfem_mgis::real;

  CatchTimeSection(ctx, "set_mgis_stuff");

  // const int nMat = 8;
  const int nMat =
      getMaterialsAttributes(*(problem.getFiniteElementDiscretizationPointer()))
          .Max();
  mfem_mgis::Profiler::Utils::Message("Nombre de matériaux : ", nMat);

  for (int i = 0; i < nMat; i++) {
    problem.addBehaviourIntegrator(ctx, "Mechanics", i + 1, p.library,
                                   p.behaviour) |
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
  if (vectors.size() != 2 * nMat) {
    throw std::invalid_argument(
        "setupMaterials : incorrect number of vectors in vector file");
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
  // pilote
  constexpr int defaultMaxNumOfIt = 5000;  // MaximumNumberOfIterations
  auto solverParameters = mfem_mgis::Parameters{};
  solverParameters.insert(
      mfem_mgis::throwing,
      mfem_mgis::Parameters{{"VerbosityLevel", params.verbosity_level}});
  solverParameters.insert(
      mfem_mgis::throwing,
      mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
  // solverParameters.insert(mfem_mgis::Parameters{{"AbsoluteTolerance", Tol}});
  // solverParameters.insert(mfem_mgis::Parameters{{"RelativeTolerance", Tol}});
  solverParameters.insert(
      mfem_mgis::throwing,
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
        mfem_mgis::throwing,
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
      exporter->SetTime(t);
      const int nMat =
          getMaterialsAttributes(*(p.getFiniteElementDiscretizationPointer()))
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
              fes_pqfs.emplace_back(m.getPartialQuadratureSpacePointer(), 1);
          mfem_mgis::computeFirstEigenStress(ctx, fes_pqf, m, ets);
          fesView.push_back(fes_pqf);
        }
        if (params.export_von_Mises_stress) {
          auto &vmes_pqf =
              vmes_pqfs.emplace_back(m.getPartialQuadratureSpacePointer(), 1);
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
                               const std::string &msg) {
  auto or_die = ctx.getFatalFailureHandler();
  p.addPostProcessing(ctx, "ParaviewExportResults", {{"OutputFileName", msg}}) |
      or_die;
  p.addPostProcessing(ctx, "MeanThermodynamicForces",
                      {{"OutputFileName", "avgStressPolycrystal"}}) |
      or_die;
}  // end timer add_postprocessing_and_outputs
