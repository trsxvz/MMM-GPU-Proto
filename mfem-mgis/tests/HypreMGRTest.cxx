/*!
 * \file   HypreMGRTest.cxx
 * \brief  Test of the multigrid reduction of hypre (`HypreMGR`), on a
 * cantilever beam in linear elasticity, the unknowns being ordered by nodes
 * and by vector dimension:
 *
 * - with the default settings of hypre (Jacobi F-relaxation), the GMRES
 *   solver of MFEM reduces the residual (on the GPU, it does not converge
 *   in a thousand iterations);
 * - with an F-relaxation by algebraic multigrid, with one level of
 *   reduction (by default, the last component being eliminated) and with
 *   two levels, the GMRES solver of MFEM converges in less than a third of
 *   the iterations of the Jacobi preconditioner (when it converges);
 * - the numbers of iterations hardly depend on the ordering of the unknowns
 *   (the components are given to hypre by markers; only the coarse solver,
 *   BoomerAMG, depends on the ordering);
 * - the block interpolation of two components eliminated at one level
 *   (unknowns ordered by vector dimension, on the host), and the block
 *   interpolation with the Jacobi relaxation eliminating one component
 *   (which uses the Jacobi relaxation of the GPU when hypre computes on a
 *   GPU) reduce the residual;
 * - the GMRES solver of hypre is preconditioned with the options of hypre
 *   (values which hypre also computes on a GPU);
 * - a second operator, the first one times two, gives the same numbers of
 *   iterations, and operators of other sizes, one of them with a process
 *   owning no unknown, are solved by the same preconditioner;
 * - invalid components and options are refused, the options of hypre 3.2.0
 *   which give wrong results or stop the program, the block options
 *   grouping unknowns of different nodes, and, when hypre computes on a
 *   GPU, the values of the options that hypre only computes on the host.
 * \date   15/10/2026
 */

#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <algorithm>
#include <type_traits>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/PreconditionerFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "ElasticitySystem.hxx"

using mfem_mgis::tests::ElasticitySystem;

//! \brief number of failed checks
static int failures = 0;
//! \brief rank of the process
static int rank = 0;
//! \brief if true, the computations are in double precision
static constexpr bool double_precision =
    std::is_same_v<mfem_mgis::real, double>;
/*!
 * \brief relative tolerance of the Krylov solvers. The GMRES solvers stop on
 * the preconditioned residual, which, with two levels of reduction, is
 * about four orders of magnitude smaller than the residual
 */
static constexpr double tolerance = double_precision ? 1e-12 : 1e-5;
/*!
 * \brief largest relative residual of a converged resolution (in single
 * precision, the GMRES solver of MFEM stops with residuals close to 1e-3
 * with the Jacobi preconditioner)
 */
static constexpr double maximum_residual = double_precision ? 1e-6 : 1e-2;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    if (rank == 0) {
      std::cerr << "check failed: " << msg << '\n';
    }
    ++failures;
  }
}  // end of check

/*!
 * \return if two numbers of iterations are close: rounding errors, and the
 * coarse solver of the multigrid reduction, which depends on the ordering of
 * the unknowns, change them slightly
 */
static bool areClose(const std::optional<int>& a, const std::optional<int>& b) {
  if ((!a.has_value()) || (!b.has_value())) {
    return false;
  }
  return 4 * std::abs(*a - *b) <= std::max(*a, *b);
}  // end of areClose

//! \brief result of a resolution
struct Resolution {
  //! \brief if true, the solver reported its convergence
  bool converged = false;
  //! \brief number of iterations
  int iterations = -1;
  //! \brief relative residual
  double residual = 0;
  //! \return the number of iterations, if it converged with a small residual
  std::optional<int> success() const {
    if (this->converged && (this->residual < maximum_residual)) {
      return this->iterations;
    }
    return {};
  }
};

/*!
 * \return the result of a resolution
 * \param[in] s: system
 * \param[in] a: matrix
 * \param[in] solver: name of the Krylov solver
 * \param[in] p: parameters of the preconditioner
 * \param[in] name: name of the resolution
 */
static Resolution solve(ElasticitySystem<true>& s,
                        const mfem::HypreParMatrix& a,
                        const char* const solver,
                        const mfem_mgis::Parameters& p,
                        const std::string& name) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto parameters =
      Parameters{{"MaximumNumberOfIterations", 1000}, {"Preconditioner", p}};
  if (std::string{solver}.starts_with("Hypre")) {
    parameters.insert(throwing, {{"Tolerance", tolerance}});
  } else {
    parameters.insert(throwing, {{"RelativeTolerance", tolerance},
                                 {"AbsoluteTolerance", 0.}});
  }
  auto ls = LinearSolverFactory<true>::getFactory().generate(
      ctx, solver, s.fespace, parameters);
  if (isInvalid(ls)) {
    throw std::runtime_error(ctx.getErrorMessage());
  }
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.linear_solver->SetOperator(a);
  ls.linear_solver->Mult(s.B, x);
  const auto r = s.getRelativeResidual(a, x);
  const auto converged = hasConverged(*(ls.linear_solver));
  const auto iterations =
      getNumberOfIterationsAtConvergence(*(ls.linear_solver)).value_or(-1);
  if (rank == 0) {
    std::cout << name << ": " << iterations << " iterations, residual " << r
              << (converged ? "" : " (not converged)") << '\n';
  }
  return {.converged = converged, .iterations = iterations, .residual = r};
}  // end of solve

//! \return the parameters of the preconditioner
static mfem_mgis::Parameters makeMGR(const mfem_mgis::Parameters& options) {
  return mfem_mgis::Parameters{{"Name", "HypreMGR"}, {"Options", options}};
}  // end of makeMGR

//! \brief numbers of iterations with an F-relaxation by algebraic multigrid
struct Iterations {
  //! \brief one level of reduction
  std::optional<int> one;
  //! \brief two levels of reduction
  std::optional<int> two;
};

/*!
 * \return the numbers of iterations of the GMRES solver of MFEM with an
 * F-relaxation by algebraic multigrid, with one and two levels of reduction
 * \param[in] s: system
 * \param[in] ordering: name of the ordering of the unknowns
 */
static Iterations run(ElasticitySystem<true>& s, const std::string& ordering) {
  using namespace mfem_mgis;
  const auto jacobi =
      solve(s, s.A, "GMRESSolver", Parameters{{"Name", "HypreDiagScale"}},
            "GMRESSolver, HypreDiagScale, " + ordering)
          .success();
  // in single precision, the GMRES solver of MFEM preconditioned by Jacobi
  // may not converge in a thousand iterations
  check(jacobi.has_value() || (!double_precision),
        ordering + ": resolution with Jacobi");
  // default settings of hypre: the residual is reduced
  check(solve(s, s.A, "GMRESSolver", makeMGR(Parameters{}),
              "GMRESSolver, HypreMGR (default settings), " + ordering)
                .residual < maximum_residual,
        ordering + ": resolution with the default settings");
  // F-relaxation by algebraic multigrid
  const auto amg = Parameters{{"RelaxationType", 2}};
  const auto faster = [&jacobi](const std::optional<int>& i) {
    return i.has_value() && ((!jacobi.has_value()) || (3 * (*i) < *jacobi));
  };
  const auto one = solve(s, s.A, "GMRESSolver", makeMGR(amg),
                         "GMRESSolver, HypreMGR (one level), " + ordering)
                       .success();
  check(faster(one), ordering + ": resolution with one level");
  const auto levels =
      std::vector<Parameter>{Parameter{std::vector<Parameter>{0, 1}},
                             Parameter{std::vector<Parameter>{0}}};
  const auto two = solve(s, s.A, "GMRESSolver",
                         makeMGR(Parameters{{"RelaxationType", 2},
                                            {"CoarseComponents", levels}}),
                         "GMRESSolver, HypreMGR (two levels), " + ordering)
                       .success();
  check(faster(two), ordering + ": resolution with two levels");
  // block interpolation of the two components eliminated at one level:
  // hypre groups consecutive unknowns, which are the ones of a node with
  // the ordering by vector dimension only (refused with the ordering by
  // nodes, and on a GPU, where hypre accesses invalid memory, see
  // checkRefusals)
  if ((ordering == "byVDIM") && (!isHypreUsingGPU())) {
    check(solve(s, s.A, "GMRESSolver",
                makeMGR(Parameters{
                    {"RelaxationType", 2},
                    {"InterpolationType", 12},
                    {"CoarseComponents", std::vector<Parameter>{0}}}),
                "GMRESSolver, HypreMGR (block interpolation), " + ordering)
              .success()
              .has_value(),
          ordering + ": resolution with the block interpolation");
  }
  // the block interpolation with the Jacobi relaxation, eliminating one
  // component: on a GPU, the Jacobi relaxation of the GPU is used (hypre
  // 3.2.0 would apply the one of the host)
  check(
      solve(s, s.A, "GMRESSolver",
            makeMGR(Parameters{{"InterpolationType", 12}}),
            "GMRESSolver, HypreMGR (block interpolation, Jacobi), " + ordering)
              .residual < maximum_residual,
      ordering +
          ": resolution with the block interpolation and the "
          "Jacobi relaxation");
  // the Krylov solver of hypre, and the options of hypre
  check(solve(s, s.A, "HypreGMRES",
              makeMGR(Parameters{{"RelaxationType", 2},
                                 {"NumberOfRelaxationSweeps", 2},
                                 {"InterpolationType", 12},
                                 {"RestrictionType", 2},
                                 {"VerbosityLevel", 0}}),
              "HypreGMRES, HypreMGR, " + ordering)
            .success()
            .has_value(),
        ordering + ": resolution with the GMRES solver of hypre");
  // a second operator, the first one times two: the preconditioned system
  // does not change
  auto a2 = mfem::HypreParMatrix(s.A);
  a2 *= 2.;
  check(areClose(solve(s, a2, "GMRESSolver", makeMGR(amg),
                       "GMRESSolver, HypreMGR, second operator, " + ordering)
                     .success(),
                 one),
        ordering + ": resolution with a second operator");
  return {one, two};
}  // end of run

/*!
 * \brief checks a preconditioner given operators of different sizes (the
 * markers of the components and the solver of hypre are built again), the
 * second system having a process owning no unknown
 * \param[in] s1: first system
 * \param[in] s2: second system
 */
static void checkOtherOperators(ElasticitySystem<true>& s1,
                                ElasticitySystem<true>& s2) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto ls = LinearSolverFactory<true>::getFactory().generate(
      ctx, "GMRESSolver", s1.fespace,
      Parameters{
          {"RelativeTolerance", tolerance},
          {"AbsoluteTolerance", 0.},
          {"MaximumNumberOfIterations", 1000},
          {"Preconditioner", makeMGR(Parameters{{"RelaxationType", 2}})}});
  if (isInvalid(ls)) {
    throw std::runtime_error(ctx.getErrorMessage());
  }
  auto i = 0;
  for (auto* const s : {&s1, &s2, &s1}) {
    auto x = mfem::Vector(s->B.Size());
    x.UseDevice(true);
    x = 0.;
    ls.linear_solver->SetOperator(s->A);
    ls.linear_solver->Mult(s->B, x);
    const auto r = s->getRelativeResidual(x);
    const auto iterations =
        getNumberOfIterationsAtConvergence(*(ls.linear_solver)).value_or(-1);
    if (rank == 0) {
      std::cout << "GMRESSolver, HypreMGR, operator " << i << ": " << iterations
                << " iterations, residual " << r << '\n';
    }
    check(hasConverged(*(ls.linear_solver)) && (r < maximum_residual),
          "resolution with the operator " + std::to_string(i));
    ++i;
  }
}  // end of checkOtherOperators

//! \brief checks the refusals of invalid options
static void checkRefusals(ElasticitySystem<true>& s,
                          ElasticitySystem<true>& s2) {
  using namespace mfem_mgis;
  const auto refused_by = [](ElasticitySystem<true>& system,
                             const Parameters& options,
                             const std::string& msg) {
    auto ctx = Context{};
    const auto p = PreconditionerFactory<true>::getFactory().generate(
        ctx, "HypreMGR", system.fespace, options);
    const auto error = std::string{ctx.getErrorMessage()};
    check((p == nullptr) && (error.find(msg) != std::string::npos),
          "invalid options not refused ('" + error + "')");
  };
  const auto refused = [&s, &refused_by](const Parameters& options,
                                         const std::string& msg) {
    refused_by(s, options, msg);
  };
  refused(Parameters{{"CoarseComponents", std::vector<Parameter>{0, 3}}},
          "invalid component");
  refused(Parameters{{"CoarseComponents", std::vector<Parameter>{0, 0}}},
          "invalid component");
  refused(Parameters{{"CoarseComponents", std::vector<Parameter>{0, 1, 2}}},
          "strict subset");
  refused(Parameters{{"CoarseComponents", std::vector<Parameter>{}}},
          "strict subset");
  refused(Parameters{{"CoarseComponents",
                      std::vector<Parameter>{
                          Parameter{std::vector<Parameter>{0, 1}},
                          Parameter{std::vector<Parameter>{2}}}}},
          "strict subset");
  refused(Parameters{{"CoarseComponents",
                      std::vector<Parameter>{
                          Parameter{std::vector<Parameter>{0, 1}}, 0}}},
          "invalid option 'CoarseComponents'");
  refused(Parameters{{"InterpolationType", 15}}, "invalid value");
  refused(Parameters{{"VerbosityLevel", 16}}, "invalid value");
  // block interpolation of two components with the ordering by nodes (the
  // system of this function)
  refused(Parameters{{"InterpolationType", 12},
                     {"CoarseComponents", std::vector<Parameter>{0}}},
          "ordered by vector dimension");
  // options of hypre 3.2.0 giving wrong results or stopping the program
  // (the second system has unknowns ordered by vector dimension)
  if (!isHypreUsingGPU()) {
    refused(Parameters{{"InterpolationType", 13}}, "wrong results");
    refused(Parameters{{"RestrictionType", 14}}, "invalid values");
  }
  refused_by(s2,
             Parameters{{"RestrictionType", 12},
                        {"CoarseComponents", std::vector<Parameter>{0}}},
             "stops the program");
  if (isHypreUsingGPU()) {
    refused_by(s2,
               Parameters{{"InterpolationType", 12},
                          {"CoarseComponents", std::vector<Parameter>{0}}},
               "accesses invalid memory");
  }
  refused(Parameters{{"RelaxationType", 10}}, "invalid value");
  if (isHypreUsingGPU()) {
    // values that hypre only computes on the host
    refused(Parameters{{"RelaxationType", 3}}, "computing on the GPU");
    refused(Parameters{{"InterpolationType", 3}}, "computing on the GPU");
    refused(Parameters{{"RestrictionType", 3}}, "computing on the GPU");
  }
  refused(Parameters{{"CoarsestLevelSolver", "BoomerAMG"}},
          "CoarsestLevelSolver");
}  // end of checkRefusals

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  auto mesh = mfem::Mesh::MakeCartesian3D(12, 3, 3, mfem::Element::HEXAHEDRON,
                                          4., 1., 1.);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    auto s1 = ElasticitySystem<true>(mesh, fec, mfem::Ordering::byNODES);
    auto s2 = ElasticitySystem<true>(mesh, fec, mfem::Ordering::byVDIM);
    const auto i1 = run(s1, "byNODES");
    const auto i2 = run(s2, "byVDIM");
    check(areClose(i1.one, i2.one) && areClose(i1.two, i2.two),
          "the ordering of the unknowns changes the number of iterations");
    checkRefusals(s1, s2);
    // a smaller mesh, all the elements on the first process
    auto mesh2 = mfem::Mesh::MakeCartesian3D(8, 2, 2, mfem::Element::HEXAHEDRON,
                                             4., 1., 1.);
    const auto partitioning =
        std::vector<int>(static_cast<std::size_t>(mesh2.GetNE()), 0);
    auto s3 = ElasticitySystem<true>(mesh2, fec, mfem::Ordering::byNODES,
                                     partitioning.data());
    checkOtherOperators(s1, s3);
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (rank == 0) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
