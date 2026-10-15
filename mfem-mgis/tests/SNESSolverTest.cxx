/*!
 * \file   SNESSolverTest.cxx
 * \brief  Test of the non linear solver based on the SNES of PETSc (`SNES`
 * non linear solver), in parallel: uniaxial tension of a cube made of a
 * plastic material with a linear hardening, in two time steps:
 *
 * - Newton's method of SNES with its backtracking line search (default),
 *   without line search, with a trust region, with the linear solver of the
 *   problem (MUMPS) applied once or preconditioning a Krylov solver of
 *   PETSc, and with the Krylov solver and the preconditioner of PETSc,
 *   given by the option `LinearSolver` or as the linear solver of the
 *   problem, gives the solution of the Newton solver of mfem-mgis;
 * - Newton's method of SNES without line search, and the trust region
 *   method, make the iterations of the Newton solver (same criterion),
 *   within two iterations, and with the line search at most as many plus
 *   two;
 * - with the tangent operator applied without assembly, SNES gives the
 *   solution of the Newton solver with the same operator, the linear
 *   solvers of PETSc and the trust region method being refused (they
 *   require an assembled operator);
 * - with a relative tolerance below the floor of the rounding errors, the
 *   resolution is stopped by the stagnation of the residual and accepted;
 * - a failure of the integration of the behaviours at the Newton step is
 *   reported to SNES, whose line search reduces the step: a traction of
 *   a cube made of a linear elastic material whose integration fails above
 *   a trace of the strain which the solution exceeds, the initial guess
 *   being below it, fails after at least one iteration;
 * - a failure of the linear solver stops the resolution at the first
 *   iteration;
 * - the invalid options are refused, and SNES in sequential.
 *
 * With a GPU, SNES and the solvers of PETSc compute on the host, the
 * vectors being copied from and to the vectors of mfem-mgis on the device.
 * \date   15/10/2026
 */

#include <array>
#include <cmath>
#include <tuple>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx"
#include "Tolerance.hxx"

//! \brief number of failed checks
static int failures = 0;
//! \brief rank of the process
static int rank = 0;
/*!
 * \brief tolerance on the relative difference between the solutions of
 * SNES and the one of the Newton solver
 */
static constexpr auto solution_tolerance =
    mfem_mgis::unit_tests::tolerance(1e-8, 3e-4);

static void check(const bool b, const std::string& msg) {
  if (!b) {
    if (rank == 0) {
      std::cerr << "check failed: " << msg << '\n';
    }
    ++failures;
  }
}  // end of check

//! \brief results of a computation
struct Results {
  //! \brief status
  bool status = false;
  //! \brief unknowns at the end of the computation (local part)
  std::vector<mfem_mgis::real> u;
  //! \brief total number of iterations
  mfem_mgis::size_type iterations = 0;
  //! \brief number of iterations of the time step which failed, if any
  mfem_mgis::size_type failed_iterations = 0;
};

//! \brief settings of a computation
struct Settings {
  //! \brief value of the `NonLinearSolver` parameter, the Newton solver
  //! being used if empty
  mfem_mgis::Parameters solver;
  //! \brief linear solver
  const char* linear_solver = "MUMPSSolver";
  //! \brief options of the linear solver
  mfem_mgis::Parameters linear_solver_options;
  //! \brief relative tolerance of the non linear solver
  mfem_mgis::real relative_tolerance = 1e-12;
  //! \brief if true, the tangent operator is applied without assembly
  bool matrix_free = false;
};

/*!
 * \return the results of the computation
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] s: settings
 */
static Results compute(const char* const mesh,
                       const char* const library,
                       const Settings& s) {
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
      {"Parallel", true}};
  if (s.matrix_free) {
    params.insert(throwing, {{"BehaviourIntegration", "Batch"},
                             {"GradientsComputation", "Kernels"},
                             {"InnerForcesComputation", "Kernels"},
                             {"TangentOperatorComputation", "MatrixFree"}});
  }
  if (!s.solver.empty()) {
    params.insert(throwing, {{"NonLinearSolver", s.solver}});
  }
  auto problem = construct<NonLinearEvolutionProblem>(ctx, params) | or_die;
  problem.addBehaviourIntegrator(ctx, "Mechanics", "Cube", library,
                                 "IsotropicLinearHardeningPlasticity") |
      or_die;
  auto& m = problem.getMaterial(ctx, "Cube", 0) | or_die;
  for (auto* st : {&m.s0, &m.s1}) {
    mgis::behaviour::setExternalStateVariable(ctx, *st, "Temperature", 293.15) |
        or_die;
  }
  for (const auto& [boundary, component] :
       {std::pair{"Xmin", 0}, std::pair{"Ymin", 1}, std::pair{"Zmin", 2}}) {
    problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", boundary}, {"Component", component}}) |
        or_die;
  }
  problem.addUniformDirichletBoundaryCondition(
      ctx,
      {{"Boundary", "Xmax"},
       {"Component", 0},
       {"LoadingEvolution", [](const real t) noexcept { return 1e-2 * t; }}}) |
      or_die;
  problem.setLinearSolver(ctx, s.linear_solver, s.linear_solver_options) |
      or_die;
  problem.setSolverParameters(
      ctx, Parameters{{"VerbosityLevel", 0},
                      {"RelativeTolerance", s.relative_tolerance},
                      {"AbsoluteTolerance", 0.},
                      {"MaximumNumberOfIterations", 100}}) |
      or_die;
  auto r = Results{};
  constexpr auto nsteps = 2;
  const auto dt = real{1} / nsteps;
  auto t = real{0};
  for (int i = 0; i != nsteps; ++i) {
    auto c = Context{};
    const auto o = problem.solve(c, t, dt);
    if (!o.status) {
      r.failed_iterations = o.iterations;
      return r;
    }
    r.iterations += o.iterations;
    problem.update(ctx) | or_die;
    t += dt;
  }
  const auto& u = problem.getUnknowns(ets);
  r.u.assign(u.HostRead(), u.HostRead() + u.Size());
  r.status = true;
  return r;
}  // end of compute

/*!
 * \return the results of a traction of a cube made of the
 * `FailingIsotropicLinearElasticity` behaviour (whose integration fails
 * above a trace of the strain of 1e-3), in one time step: the solution has
 * a trace of the strain of 1.9e-3, the initial guess (no displacement) a
 * null one
 * \param[in] mesh: mesh of a cube
 * \param[in] library: behaviour library
 * \param[in] solver: value of the `NonLinearSolver` parameter
 */
static Results computeFailingIntegration(const char* const mesh,
                                         const char* const library,
                                         const mfem_mgis::Parameters& solver) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  auto problem =
      construct<NonLinearEvolutionProblem>(
          ctx,
          Parameters{
              {"MeshFileName", mesh},
              {"FiniteElementFamily", "H1"},
              {"FiniteElementOrder", 1},
              {"UnknownsSize", 3},
              {"NumberOfUniformRefinements", 1},
              {"Materials", Parameters{{"Cube", 1}}},
              {"Boundaries",
               Parameters{{"Ymin", 1}, {"Zmin", 2}, {"Xmax", 3}, {"Xmin", 5}}},
              {"Hypothesis", "Tridimensional"},
              {"NonLinearSolver", solver},
              {"Parallel", true}}) |
      or_die;
  problem.addBehaviourIntegrator(ctx, "Mechanics", "Cube", library,
                                 "FailingIsotropicLinearElasticity") |
      or_die;
  auto& m = problem.getMaterial(ctx, "Cube", 0) | or_die;
  for (auto* st : {&m.s0, &m.s1}) {
    mgis::behaviour::setExternalStateVariable(ctx, *st, "Temperature", 293.15) |
        or_die;
  }
  for (const auto& [boundary, component] :
       {std::pair{"Xmin", 0}, std::pair{"Ymin", 1}, std::pair{"Zmin", 2}}) {
    problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", boundary}, {"Component", component}}) |
        or_die;
  }
  // uniaxial stress of 1 GPa: trace of the strain 1e9 (1 - 2 nu) / E,
  // nu = 0.3, E = 208 GPa
  problem.addBoundaryCondition(
      ctx, std::make_unique<UniformImposedPressureBoundaryCondition>(
               problem.getFiniteElementDiscretizationPointer(), "Xmax",
               [](const real t) noexcept { return -1e9 * t; })) |
      or_die;
  problem.setLinearSolver(ctx, "MUMPSSolver", Parameters{}) | or_die;
  problem.setSolverParameters(ctx,
                              Parameters{{"VerbosityLevel", 0},
                                         {"RelativeTolerance", 1e-12},
                                         {"AbsoluteTolerance", 0.},
                                         {"MaximumNumberOfIterations", 20}}) |
      or_die;
  auto r = Results{};
  auto c = Context{};
  const auto o = problem.solve(c, 0, 1);
  r.status = o.status;
  r.iterations = o.iterations;
  r.failed_iterations = o.status ? 0 : o.iterations;
  return r;
}  // end of computeFailingIntegration

/*!
 * \return the maximum difference between two arrays, relative to the
 * maximum absolute value of the first one, on all the processes
 */
static double getRelativeDifference(const std::vector<mfem_mgis::real>& a,
                                    const std::vector<mfem_mgis::real>& b) {
  auto values = std::array<double, 2>{0, 0};
  if (a.size() != b.size()) {
    values[0] = 1;
  } else {
    for (std::size_t i = 0; i != a.size(); ++i) {
      values[0] =
          std::max(values[0], static_cast<double>(std::abs(a[i] - b[i])));
      values[1] = std::max(values[1], static_cast<double>(std::abs(a[i])));
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, values.data(), 2, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  return values[1] <= 0 ? values[0] : values[0] / values[1];
}  // end of getRelativeDifference

/*!
 * \brief an invalid option of SNES is refused
 * \param[in] mesh: mesh
 * \param[in] options: options of SNES
 * \param[in] message: part of the expected error message
 * \param[in] parallel: parallel computation
 */
static void checkRefusal(const char* const mesh,
                         const mfem_mgis::Parameters& options,
                         const std::string& message,
                         const bool parallel = true) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  const auto o = construct<NonLinearEvolutionProblem>(
      ctx, Parameters{{"MeshFileName", mesh},
                      {"FiniteElementFamily", "H1"},
                      {"FiniteElementOrder", 1},
                      {"UnknownsSize", 3},
                      {"Hypothesis", "Tridimensional"},
                      {"NonLinearSolver", Parameters{{"SNES", options}}},
                      {"Parallel", parallel}});
  check(isInvalid(o), message + ": not refused");
  // reading the error message clears it
  const auto e = ctx.getErrorMessage();
  check(e.find(message) != std::string::npos,
        message + ": unexpected error message (" + e + ")");
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
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  auto print = [](const std::string& msg) {
    if (rank == 0) {
      std::cout << msg << '\n';
    }
  };
  const auto reference = compute(mesh, library, Settings{});
  check(reference.status, "Newton: resolution failed");
  print("Newton: " + std::to_string(reference.iterations) + " iterations");
  // settings of a computation with SNES, given its options
  const auto snes = [](const Parameters& options) {
    auto s = Settings{};
    s.solver = Parameters{{"SNES", options}};
    return s;
  };
  // a computation with SNES, compared to a reference
  auto run = [mesh, library, &print](const std::string& name,
                                     const Results& ref, const Settings& s) {
    const auto r = compute(mesh, library, s);
    const auto d = getRelativeDifference(ref.u, r.u);
    print(name + ": " + std::to_string(r.iterations) +
          " iterations, relative difference with Newton " + std::to_string(d));
    check(r.status, name + ": resolution failed");
    check(d < solution_tolerance, name + ": solution differs from Newton's");
    return r;
  };
  const auto line_search = run("SNES, line search", reference, snes({}));
  check(line_search.iterations <= reference.iterations + 2,
        "SNES, line search: more iterations than the Newton solver");
  const auto newton =
      run("SNES, no line search", reference,
          snes({{"Options", Parameters{{"snes_linesearch_type", "none"}}}}));
  check((newton.iterations + 2 >= reference.iterations) &&
            (newton.iterations <= reference.iterations + 2),
        "SNES, no line search: not the iterations of the Newton solver");
  const auto trust_region =
      run("SNES, trust region", reference,
          snes({{"Options", Parameters{{"snes_type", "newtontr"}}}}));
  check((trust_region.iterations + 2 >= reference.iterations) &&
            (trust_region.iterations <= reference.iterations + 2),
        "SNES, trust region: not the iterations of the Newton solver");
  // the linear solver of the problem preconditioning a Krylov solver of
  // PETSc
  std::ignore = run("SNES, GMRES of PETSc preconditioned by MUMPS", reference,
                    snes({{"Options", Parameters{{"ksp_type", "gmres"},
                                                 {"ksp_rtol", 1e-12}}}}));
  // the Krylov solver and the preconditioner of PETSc
  const auto petsc_options = Parameters{{"ksp_type", "gmres"},
                                        {"ksp_rtol", 1e-12},
                                        {"ksp_max_it", 1000},
                                        {"pc_type", "gamg"}};
  std::ignore =
      run("SNES, GMRES and GAMG of PETSc", reference,
          snes({{"LinearSolver", "PETSc"}, {"Options", petsc_options}}));
  // the same solvers of PETSc as the linear solver of the problem
  {
    auto s = snes({});
    s.linear_solver = "PETScLinearSolver";
    s.linear_solver_options = Parameters{{"Options", petsc_options}};
    std::ignore =
        run("SNES, GMRES and GAMG of PETSc, linear solver of the problem",
            reference, s);
  }
  // tangent operator applied without assembly, compared to the Newton
  // solver with the same operator (conjugate gradient and Jacobi smoother)
  {
    const auto matrix_free = [](Settings s) {
      s.linear_solver = "CGSolver";
      s.linear_solver_options = Parameters{
          {"VerbosityLevel", 0},
          {"RelativeTolerance", 1e-14},
          {"AbsoluteTolerance", 0.},
          {"MaximumNumberOfIterations", 10000},
          {"Preconditioner", Parameters{{"Name", "OperatorJacobiSmoother"}}}};
      s.matrix_free = true;
      return s;
    };
    const auto mf_reference = compute(mesh, library, matrix_free(Settings{}));
    check(mf_reference.status, "Newton, matrix free: resolution failed");
    std::ignore = run("SNES, matrix free", mf_reference, matrix_free(snes({})));
    {
      const auto r = compute(mesh, library,
                             matrix_free(snes({{"LinearSolver", "PETSc"}})));
      check((!r.status) && (r.failed_iterations == 0),
            "SNES, matrix free: the linear solvers of PETSc are not refused");
    }
    // the trust region methods need the transposed tangent operator
    const auto r = compute(
        mesh, library,
        matrix_free(
            snes({{"Options", Parameters{{"snes_type", "newtontr"}}}})));
    check((!r.status) && (r.failed_iterations == 0),
          "SNES, matrix free: the trust region method is not refused");
  }
  // a relative tolerance below the floor of the rounding errors: the
  // resolution is stopped by the stagnation of the residual, and accepted
  {
    auto s = snes({});
    s.relative_tolerance = 1e-20;
    const auto r = compute(mesh, library, s);
    const auto d = getRelativeDifference(reference.u, r.u);
    print("SNES, stagnation: " + std::to_string(r.iterations) +
          " iterations, relative difference with Newton " + std::to_string(d));
    check(r.status && (r.iterations < 50) && (d < solution_tolerance),
          "SNES: resolution not accepted at the floor of the rounding errors");
  }
  // a failure of the integration of the behaviours at the Newton step,
  // reported to SNES: the line search reduces the step, and the resolution
  // fails after at least one iteration
  {
    const auto r = computeFailingIntegration(
        mesh, library, Parameters{{"SNES", Parameters{}}});
    print(std::string{"SNES, failure of the integration: "} +
          (r.status ? "converged" : "failed") + " after " +
          std::to_string(r.failed_iterations) + " iteration(s)");
    check((!r.status) && (r.failed_iterations >= 1),
          "SNES: failure of the integration not handled by the line search");
  }
  // a failure of the linear solver (the conjugate gradient stopped after
  // two iterations) stops the resolution at the first iteration
  {
    auto s = snes({});
    s.linear_solver = "CGSolver";
    s.linear_solver_options = Parameters{{"VerbosityLevel", 0},
                                         {"RelativeTolerance", 1e-12},
                                         {"MaximumNumberOfIterations", 2}};
    const auto r = compute(mesh, library, s);
    print(std::string{"SNES, failure of the linear solver: "} +
          (r.status ? "converged" : "failed") + " after " +
          std::to_string(r.failed_iterations) + " iteration(s)");
    check((!r.status) && (r.failed_iterations <= 1),
          "SNES: failure of the linear solver not reported");
  }
  // refusals
  checkRefusal(mesh, Parameters{{"LinearSolver", "Ginkgo"}},
               "invalid linear solver");
  checkRefusal(mesh, Parameters{{"Options", Parameters{{"pc_type", "gamg"}}}},
               "is refused with the linear solver of the problem");
  checkRefusal(mesh,
               Parameters{{"Options", Parameters{{"-snes_type", "newtontr"}}}},
               "without the leading dash");
  checkRefusal(
      mesh,
      Parameters{{"Options", Parameters{{"snes_convergence_test", "skip"}}}},
      "the convergence test is the one of mfem-mgis");
  checkRefusal(mesh, Parameters{{"StagnationRelativeTolerance", 1.}},
               "StagnationRelativeTolerance");
  checkRefusal(mesh, Parameters{{"MaximumNumberOfIterations", 0}},
               "maximum number of iterations");
  checkRefusal(mesh, Parameters{{"MaximumNewtonStep", 1.}},
               "MaximumNewtonStep");
  checkRefusal(mesh, Parameters{}, "only available in parallel", false);
  print(failures == 0 ? "success" : "failure");
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
