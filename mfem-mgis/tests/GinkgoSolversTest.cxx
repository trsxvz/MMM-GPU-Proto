/*!
 * \file   GinkgoSolversTest.cxx
 * \brief  Test of the linear solvers and of the preconditioners of Ginkgo
 * built by the factories, on the Laplace equation on a cube, in sequential
 * and in parallel (the generic test `LinearSolverFactoriesTest` solves it
 * with each entry and its default options):
 *
 * - a preconditioner of Ginkgo is given to a solver of Ginkgo by its factory
 *   (the local ones being wrapped in the Schwarz preconditioner of Ginkgo in
 *   parallel), and the algebraic multigrid of Ginkgo reduces the number of
 *   iterations compared to the Jacobi preconditioner;
 * - the options of mfem-mgis (dimension of the Krylov space) and the keys of
 *   the `Configuration` option reach Ginkgo;
 * - a preconditioner of another library is wrapped, and receives the
 *   operators of different sizes given to the solver;
 * - a solver of Ginkgo described by a whole configuration of Ginkgo,
 *   including its preconditioner and its stopping criteria, and a
 *   preconditioner of Ginkgo described by a whole configuration, given to a
 *   solver of MFEM;
 * - the relative tolerance can be read and changed (see
 *   `getRelativeTolerance` and `setRelativeTolerance`), and a resolution
 *   stopped by the maximum number of iterations is reported as not
 *   converged;
 * - the iterative refinement with a preconditioner of Ginkgo, the Krylov
 *   dimension of GCR and CB-GMRES, the storage precision of CB-GMRES, an
 *   operator applied without assembly (refused by BiCG);
 * - the direct solver of Ginkgo (sequential), with the LU and the Cholesky
 *   factorizations, and the multigrid preconditioner with two cycles (keys of
 *   the entries which the `Configuration` option can redefine);
 * - the L1 smoother of the Schwarz preconditioner (parallel);
 * - invalid configurations are refused: unknown keys, keys defined twice,
 *   stopping criteria or preconditioner given twice, generic solver of a
 *   type which is not an iterative solver, L1 smoother given by a generic
 *   configuration.
 * \date   15/10/2026
 */

#include <cmath>
#include <string>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <ginkgo/ginkgo.hpp>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "LaplaceSystem.hxx"

using mfem_mgis::tests::LaplaceSystem;

//! \brief number of failed checks
static int failures = 0;
//! \brief rank of the process
static int rank = 0;
//! \brief parallel computations
static bool parallel = false;
//! \brief tolerance of the resolutions
static constexpr double tolerance =
    std::is_same_v<mfem::real_t, double> ? 1e-10 : 1e-5;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

//! \brief outcome of a resolution
struct Resolution {
  //! \brief if true, the resolution converged
  bool converged = false;
  //! \brief number of iterations
  int iterations = 0;
  //! \brief norm of the residual relative to the right hand side
  double residual = 1;
};

/*!
 * \return true if a resolution converged in at least one iteration with a
 * relative residual lower than 100 times the tolerance
 * \param[in] r: resolution
 */
static bool isSolved(const Resolution& r) {
  return r.converged && (r.iterations > 0) && (r.residual < 100 * tolerance);
}  // end of isSolved

/*!
 * \return true if a preconditioner reduces the number of iterations to less
 * than two thirds of the ones of another resolution
 * \param[in] r: resolution
 * \param[in] r0: reference resolution
 */
static bool isFaster(const Resolution& r, const Resolution& r0) {
  return (r.iterations > 0) && (3 * r.iterations < 2 * r0.iterations);
}  // end of isFaster

/*!
 * \return the linear solver of the factory, throws on failure
 * \param[in] s: system
 * \param[in] name: name of the linear solver
 * \param[in] parameters: parameters of the linear solver
 */
template <bool is_parallel>
static mfem_mgis::LinearSolverHandler generate(
    LaplaceSystem<is_parallel>& s,
    const std::string& name,
    const mfem_mgis::Parameters& parameters) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto ls = LinearSolverFactory<is_parallel>::getFactory().generate(
      ctx, name, s.fespace, parameters);
  if (isInvalid(ls)) {
    throw std::runtime_error(name + ": " + ctx.getErrorMessage());
  }
  return ls;
}  // end of generate

/*!
 * \return the error message of the generation of a linear solver which
 * shall fail, or an empty string if the generation succeeds
 * \param[in] s: system
 * \param[in] name: name of the linear solver
 * \param[in] parameters: parameters of the linear solver
 */
template <bool is_parallel>
static std::string getGenerationError(LaplaceSystem<is_parallel>& s,
                                      const std::string& name,
                                      const mfem_mgis::Parameters& parameters) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  const auto ls = LinearSolverFactory<is_parallel>::getFactory().generate(
      ctx, name, s.fespace, parameters);
  return isInvalid(ls) ? ctx.getErrorMessage() : "";
}  // end of getGenerationError

//! \return the resolution of the system by a linear solver
template <bool is_parallel>
static Resolution solve(LaplaceSystem<is_parallel>& s,
                        mfem::Solver& ls,
                        const std::string& name) {
  using namespace mfem_mgis;
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.SetOperator(s.A);
  ls.Mult(s.B, x);
  auto r = Resolution{};
  r.converged = hasConverged(ls);
  r.iterations = getNumberOfIterationsAtConvergence(ls).value_or(-1);
  r.residual = s.getRelativeResidual(x);
  if (rank == 0) {
    std::cout << name << ": " << r.iterations << " iterations, residual "
              << r.residual << (r.converged ? "" : " (not converged)") << '\n';
  }
  return r;
}  // end of solve

//! \brief operator applied without assembly, wrapping a matrix
template <typename Matrix>
struct MatrixFreeOperator final : public mfem::Operator {
  //! \param[in] m: matrix
  explicit MatrixFreeOperator(const Matrix& m)
      : mfem::Operator(m.Height(), m.Width()), matrix(m) {}
  void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
    this->matrix.Mult(x, y);
  }

 private:
  //! \brief matrix
  const Matrix& matrix;
};

/*!
 * \return the resolution of the system by a linear solver, for an operator
 * which is not the matrix of the system
 * \param[in] s: system
 * \param[in] op: operator
 * \param[in] ls: linear solver
 * \param[in] name: name of the resolution
 */
template <bool is_parallel>
static Resolution solve(LaplaceSystem<is_parallel>& s,
                        const mfem::Operator& op,
                        mfem::Solver& ls,
                        const std::string& name) {
  using namespace mfem_mgis;
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.SetOperator(op);
  ls.Mult(s.B, x);
  auto r = Resolution{};
  r.converged = hasConverged(ls);
  r.iterations = getNumberOfIterationsAtConvergence(ls).value_or(-1);
  r.residual = s.getRelativeResidual(x);
  if (rank == 0) {
    std::cout << name << ": " << r.iterations << " iterations, residual "
              << r.residual << (r.converged ? "" : " (not converged)") << '\n';
  }
  return r;
}  // end of solve

//! \return the factory of the solver of Ginkgo of a linear solver
static std::shared_ptr<gko::LinOpFactory> getFactory(
    const mfem_mgis::LinearSolverHandler& ls) {
  const auto* const s =
      dynamic_cast<const mfem::Ginkgo::GinkgoIterativeSolver*>(
          ls.linear_solver.get());
  if (s == nullptr) {
    throw std::runtime_error("not a linear solver of Ginkgo");
  }
  return s->GetFactory();
}  // end of getFactory

//! \return the parameters of the conjugate gradient
template <typename Solver>
static auto getParameters(const mfem_mgis::LinearSolverHandler& ls) {
  return gko::as<typename Solver::Factory>(getFactory(ls))->get_parameters();
}  // end of getParameters

/*!
 * \return the options of the solvers of this test
 * \param[in] n: name of the preconditioner, if any
 * \param[in] o: options of the preconditioner
 */
static mfem_mgis::Parameters getOptions(const std::string& n = "",
                                        const mfem_mgis::Parameters& o = {}) {
  using namespace mfem_mgis;
  auto p = Parameters{{"RelativeTolerance", tolerance},
                      {"MaximumNumberOfIterations", 1000}};
  if (!n.empty()) {
    p.replaceOrInsert("Preconditioner",
                      Parameters{{"Name", n}, {"Options", o}});
  }
  return p;
}  // end of getOptions

template <bool is_parallel>
static void run(LaplaceSystem<is_parallel>& s) {
  using namespace mfem_mgis;
  using real_t = mfem::real_t;
  // preconditioners of Ginkgo given by their factories
  auto jacobi = generate(s, "GinkgoCG", getOptions("GinkgoJacobi"));
  const auto rj = solve(s, *(jacobi.linear_solver), "GinkgoCG, GinkgoJacobi");
  check(isSolved(rj), "GinkgoCG, GinkgoJacobi: resolution");
  const auto pj = getParameters<gko::solver::Cg<real_t>>(jacobi).preconditioner;
  if constexpr (is_parallel) {
#if defined(MFEM_USE_MPI) && GINKGO_BUILD_MPI
    using Schwarz = gko::experimental::distributed::preconditioner::Schwarz<
        real_t, int, mfem::Ginkgo::gko_hypre_bigint>;
    check(std::dynamic_pointer_cast<const typename Schwarz::Factory>(pj) !=
              nullptr,
          "GinkgoCG, GinkgoJacobi: Schwarz preconditioner in parallel");
#endif /* defined(MFEM_USE_MPI) && GINKGO_BUILD_MPI */
  } else {
    check(std::dynamic_pointer_cast<
              const gko::preconditioner::Jacobi<real_t, int>::Factory>(pj) !=
              nullptr,
          "GinkgoCG, GinkgoJacobi: Jacobi preconditioner");
  }
  auto multigrid = generate(s, "GinkgoCG", getOptions("GinkgoMultigrid"));
  const auto rm =
      solve(s, *(multigrid.linear_solver), "GinkgoCG, GinkgoMultigrid");
  check(isSolved(rm) && isFaster(rm, rj), "GinkgoCG, GinkgoMultigrid");
  // options of mfem-mgis and keys of the configuration
  auto gmres =
      generate(s, "GinkgoGMRES",
               Parameters{{"KDim", 40},
                          {"RelativeTolerance", tolerance},
                          {"MaximumNumberOfIterations", 1000},
                          {"Configuration", Parameters{{"flexible", true}}}});
  check(isSolved(solve(s, *(gmres.linear_solver), "GinkgoGMRES")),
        "GinkgoGMRES: resolution");
  const auto pg = getParameters<gko::solver::Gmres<real_t>>(gmres);
  check((pg.krylov_dim == 40) && (pg.flexible),
        "GinkgoGMRES: parameters not given to Ginkgo");
  // the iterative refinement with a preconditioner of Ginkgo as inner
  // solver
  auto ir = generate(s, "GinkgoIR", getOptions("GinkgoJacobi"));
  check(isSolved(solve(s, *(ir.linear_solver), "GinkgoIR, GinkgoJacobi")) &&
            (getParameters<gko::solver::Ir<real_t>>(ir).solver != nullptr),
        "GinkgoIR, GinkgoJacobi: resolution");
  auto gcr = generate(s, "GinkgoGCR",
                      Parameters{{"KDim", 20},
                                 {"RelativeTolerance", tolerance},
                                 {"MaximumNumberOfIterations", 1000}});
  check(isSolved(solve(s, *(gcr.linear_solver), "GinkgoGCR")) &&
            (getParameters<gko::solver::Gcr<real_t>>(gcr).krylov_dim == 20),
        "GinkgoGCR: resolution and Krylov dimension");
  // an operator applied without assembly
  auto matrix_free = MatrixFreeOperator(s.A);
  auto mf = generate(s, "GinkgoCG",
                     Parameters{{"RelativeTolerance", tolerance},
                                {"MaximumNumberOfIterations", 1000}});
  check(isSolved(solve(s, matrix_free, *(mf.linear_solver),
                       "GinkgoCG, operator applied without assembly")),
        "GinkgoCG: resolution with an operator applied without assembly");
  // a preconditioner of another library, operators of different sizes
  auto smoother = generate(
      s, "GinkgoCG",
      getOptions("OperatorJacobiSmoother", Parameters{{"Damping", 1.}}));
  check(isSolved(solve(s, *(smoother.linear_solver),
                       "GinkgoCG, OperatorJacobiSmoother")),
        "GinkgoCG, OperatorJacobiSmoother: resolution");
  auto mesh = mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON);
  auto fec = mfem::H1_FECollection(1, 3);
  auto s2 = LaplaceSystem<is_parallel>(mesh, fec);
  check(isSolved(solve(s2, *(smoother.linear_solver),
                       "GinkgoCG, OperatorJacobiSmoother, smaller mesh")),
        "GinkgoCG, OperatorJacobiSmoother: resolution on a smaller mesh");
  if constexpr (is_parallel) {
    // the L1 smoother of the Schwarz preconditioner
    auto l1 =
        generate(s, "GinkgoCG",
                 getOptions("GinkgoJacobi", Parameters{{"L1Smoother", true}}));
    check(isSolved(solve(s, *(l1.linear_solver), "GinkgoCG, L1 smoother")),
          "GinkgoCG, GinkgoJacobi with the L1 smoother: resolution");
    auto amg = generate(
        s, "GinkgoCG",
        getOptions("HypreBoomerAMG", Parameters{{"VerbosityLevel", 0}}));
    const auto ra = solve(s, *(amg.linear_solver), "GinkgoCG, HypreBoomerAMG");
    check(isSolved(ra) && isFaster(ra, rj), "GinkgoCG, HypreBoomerAMG");
  }
  // a solver described by a whole configuration of Ginkgo
  const auto whole = Parameters{
      {"type", "solver::Cg"},
      {"criteria",
       Parameters{{"iteration", 1000}, {"relative_residual_norm", tolerance}}},
      {"preconditioner",
       is_parallel ? Parameters{{"type", "preconditioner::Schwarz"},
                                {"local_solver",
                                 Parameters{{"type", "preconditioner::Ic"}}}}
                   : Parameters{{"type", "preconditioner::Ic"}}}};
  auto configured =
      generate(s, "GinkgoSolver", Parameters{{"Configuration", whole}});
  const auto rc =
      solve(s, *(configured.linear_solver), "GinkgoSolver, IC preconditioner");
  check(isSolved(rc), "GinkgoSolver: resolution");
  check(!getRelativeTolerance(*(configured.linear_solver)).has_value(),
        "GinkgoSolver: relative tolerance given by the configuration");
  // a preconditioner described by a whole configuration of Ginkgo, given to
  // a solver of MFEM
  const auto local =
      Parameters{{"type", "preconditioner::Jacobi"}, {"max_block_size", 1}};
  const auto preconditioner =
      is_parallel ? Parameters{{"type", "preconditioner::Schwarz"},
                               {"local_solver", local}}
                  : local;
  auto mfem_cg =
      generate(s, "CGSolver",
               getOptions("GinkgoPreconditioner",
                          Parameters{{"Configuration", preconditioner}}));
  check(isSolved(solve(s, *(mfem_cg.linear_solver),
                       "CGSolver, GinkgoPreconditioner")),
        "CGSolver, GinkgoPreconditioner: resolution");
  // relative tolerance
  auto cg = generate(s, "GinkgoCG", getOptions("GinkgoJacobi"));
  const auto r1 = solve(s, *(cg.linear_solver), "GinkgoCG");
  // the tolerance given is bounded below in single precision
  const auto expected = static_cast<double>(boundRelativeTolerance(tolerance));
  const auto otol = getRelativeTolerance(*(cg.linear_solver));
  check(otol.has_value() && (std::abs(*otol - expected) < 1e-3 * expected),
        "GinkgoCG: relative tolerance not read");
  check(setRelativeTolerance(*(cg.linear_solver), 1e-3),
        "GinkgoCG: relative tolerance not set");
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  cg.linear_solver->Mult(s.B, x);
  const auto n2 =
      getNumberOfIterationsAtConvergence(*(cg.linear_solver)).value_or(-1);
  if (rank == 0) {
    std::cout << "GinkgoCG, relative tolerance 1e-3: " << n2 << " iterations\n";
  }
  check(hasConverged(*(cg.linear_solver)) && (n2 > 0) && (n2 < r1.iterations) &&
            (s.getRelativeResidual(x) < 1e-2),
        "GinkgoCG: relative tolerance not taken into account");
  // the factory rebuilt with the new tolerance, used for a new operator
  const auto r2 = solve(s2, *(cg.linear_solver), "GinkgoCG, smaller mesh");
  check(r2.converged && (r2.iterations > 0) && (r2.residual < 1e-2),
        "GinkgoCG: resolution after a change of the relative tolerance");
  // a resolution stopped by the maximum number of iterations
  auto stopped = generate(s, "GinkgoCG",
                          Parameters{{"RelativeTolerance", tolerance},
                                     {"MaximumNumberOfIterations", 2}});
  const auto rs = solve(s, *(stopped.linear_solver), "GinkgoCG, 2 iterations");
  check((!rs.converged) && (rs.iterations == 2),
        "resolution stopped by the maximum number of iterations");
  // the solvers which are only available in sequential: the direct solver,
  // CB-GMRES and BiCG
  if constexpr (!is_parallel) {
    for (const auto& f : {"factorization::Lu", "factorization::Cholesky"}) {
      auto direct = generate(
          s, "GinkgoDirect",
          Parameters{{"Configuration",
                      Parameters{{"factorization", Parameters{{"type", f}}}}}});
      auto xd = mfem::Vector(s.B.Size());
      xd.UseDevice(true);
      xd = 0.;
      direct.linear_solver->SetOperator(s.A);
      direct.linear_solver->Mult(s.B, xd);
      const auto rd = s.getRelativeResidual(xd);
      std::cout << "GinkgoDirect, " << f << ": residual " << rd << '\n';
      check(rd < 100 * tolerance, std::string{"GinkgoDirect, "} + f);
    }
    auto cbgmres = generate(
        s, "GinkgoCBGMRES",
        Parameters{
            {"KDim", 30},
            {"RelativeTolerance", tolerance},
            {"MaximumNumberOfIterations", 1000},
            {"Configuration", Parameters{{"storage_precision", "reduce2"}}}});
    const auto pcb = getParameters<gko::solver::CbGmres<real_t>>(cbgmres);
    check(isSolved(solve(s, *(cbgmres.linear_solver), "GinkgoCBGMRES")) &&
              (pcb.krylov_dim == 30) &&
              (pcb.storage_precision ==
               gko::solver::cb_gmres::storage_precision::reduce2),
          "GinkgoCBGMRES: resolution and parameters");
    // BiCG requires an assembled operator
    auto bicg = generate(s, "GinkgoBiCG",
                         Parameters{{"RelativeTolerance", tolerance},
                                    {"MaximumNumberOfIterations", 1000}});
    auto raised = false;
    try {
      bicg.linear_solver->SetOperator(matrix_free);
    } catch (std::exception&) {
      raised = true;
    }
    check(raised, "GinkgoBiCG: operator applied without assembly not refused");
  } else {
    check(getGenerationError(s, "GinkgoBiCG", Parameters{})
                  .find("only available in sequential") != std::string::npos,
          "GinkgoBiCG: not reported as only available in sequential");
  }
  // invalid configurations
  const auto refused = [&s](const std::string& n, const Parameters& p,
                            const std::string& msg) {
    return getGenerationError(s, n, p).find(msg) != std::string::npos;
  };
  check(refused("GinkgoCG",
                Parameters{{"Configuration", Parameters{{"unknown", 1}}}},
                "invalid configuration of Ginkgo"),
        "unknown key not refused");
  check(refused("GinkgoCG",
                Parameters{
                    {"Configuration", Parameters{{"type", "solver::Gmres"}}}},
                "is already defined"),
        "type redefined not refused");
  check(refused("GinkgoSolver", Parameters{}, "no configuration of Ginkgo"),
        "generic solver without configuration not refused");
  check(
      refused("GinkgoSolver",
              Parameters{{"Configuration", whole}, {"RelativeTolerance", 1e-3}},
              "the stopping criteria are given both"),
      "stopping criteria given twice not refused");
  check(refused("GinkgoSolver",
                Parameters{
                    {"Configuration", whole},
                    {"Preconditioner", Parameters{{"Name", "GinkgoJacobi"}}}},
                "the preconditioner is given both"),
        "preconditioner given twice not refused");
  check(refused("GinkgoSolver",
                Parameters{{"Configuration", whole},
                           {"Preconditioner",
                            Parameters{{"Name", "OperatorJacobiSmoother"}}}},
                "the preconditioner is given both"),
        "preconditioner given twice (other library) not refused");
  check(refused("GinkgoSolver",
                Parameters{{"Configuration",
                            Parameters{{"type", "preconditioner::Jacobi"}}}},
                "does not describe an iterative solver"),
        "generic solver of a non iterative type not refused");
  check(refused("CGSolver",
                getOptions(
                    "GinkgoPreconditioner",
                    Parameters{{"Configuration",
                                Parameters{{"type", "preconditioner::Schwarz"},
                                           {"l1_smoother", true},
                                           {"local_solver", local}}}}),
                "L1 smoother"),
        "L1 smoother in a generic configuration not refused");
  // errors in a configuration which Ginkgo does not read entirely
  check(refused("GinkgoGMRES",
                Parameters{
                    {"KDim", 10},
                    {"Configuration", Parameters{{"ortho_method", "unknown"}}}},
                "invalid configuration of Ginkgo"),
        "invalid value of a key with a Krylov dimension not refused");
  check(refused(
            "GinkgoCG",
            getOptions("GinkgoILU",
                       Parameters{{"Configuration",
                                   Parameters{{"factorization",
                                               Parameters{{"type",
                                                           "factorization::"
                                                           "ParIlu"},
                                                          {"unknown", 1}}}}}}),
            "invalid configuration of Ginkgo"),
        "invalid nested configuration of a preconditioner not refused");
  if constexpr (!is_parallel) {
    check(refused("GinkgoBiCG", getOptions("GinkgoSOR"), "can be transposed"),
          "GinkgoBiCG: preconditioner which can't be transposed not refused");
  }
  // keys of the entries which can be redefined
  auto cycles = generate(
      s, "GinkgoCG",
      getOptions("GinkgoMultigrid",
                 Parameters{{"Configuration",
                             Parameters{{"criteria",
                                         Parameters{{"type", "Iteration"},
                                                    {"max_iters", 2}}}}}}));
  check(isSolved(solve(s, *(cycles.linear_solver),
                       "GinkgoCG, GinkgoMultigrid, two cycles")),
        "GinkgoCG, GinkgoMultigrid with two cycles: resolution");
}  // end of run

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "parallel computations");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
#ifdef MFEM_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif /* MFEM_USE_MPI */
  auto mesh =
      mfem::Mesh::MakeCartesian3D(12, 12, 12, mfem::Element::HEXAHEDRON);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    if (parallel) {
#if defined(MFEM_USE_MPI) && GINKGO_BUILD_MPI
      auto s = LaplaceSystem<true>(mesh, fec);
      run<true>(s);
#elif defined(MFEM_USE_MPI)
      // the solvers of Ginkgo are only available in sequential
      if (rank == 0) {
        std::cout << "skipped: Ginkgo is built without MPI\n";
      }
      return EXIT_SUCCESS;
#else  /* MFEM_USE_MPI */
      throw std::runtime_error("parallel computations are not supported");
#endif /* MFEM_USE_MPI */
    } else {
      auto s = LaplaceSystem<false>(mesh, fec);
      run<false>(s);
    }
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (rank == 0) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
