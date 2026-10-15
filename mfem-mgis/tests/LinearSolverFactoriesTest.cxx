/*!
 * \file   LinearSolverFactoriesTest.cxx
 * \brief  Generic test of the factories of linear solvers and of
 * preconditioners: each registered linear solver, and each registered
 * preconditioner given to the GMRES solver of MFEM, solves the Laplace
 * equation on a cube, the essential degrees of freedom being eliminated
 * with a diagonal set to one. The outcome is checked according to the kind
 * of solver (see `LinearSolverDescription`): convergence up to the
 * tolerance for the iterative and direct solvers, decrease of the residual
 * for the approximate solvers. On the linear elasticity, the checks of the
 * iterative and approximate solvers are weaker (see `checkResolution`): a
 * converged resolution has a residual below `elasticity_tolerance`, the
 * entries which do not converge on this problem are listed (in each
 * precision), and the approximate solvers only give a finite non zero
 * solution; the direct solvers are checked as on the Laplace equation. The
 * test also
 * checks the descriptions, the preconditioners refused by the linear solvers,
 * and, in sequential, the message of the solvers only available in parallel.
 *
 * The entries which only compute on the GPU are skipped when MFEM does not
 * compute on a CUDA GPU, and the ones which only compute on the host when
 * hypre computes on a GPU. With the `--parallel` option, the parallel
 * factories are checked. With the `--elasticity` option, the system solved
 * is a cantilever beam in linear elasticity (three components per node)
 * instead of the Laplace equation.
 * \date   15/10/2026
 */

#include <map>
#include <set>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>
#include <algorithm>
#include <type_traits>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/PreconditionerFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "LaplaceSystem.hxx"
#include "ElasticitySystem.hxx"

using mfem_mgis::tests::ElasticitySystem;
using mfem_mgis::tests::LaplaceSystem;

//! \brief number of failed checks
static int failures = 0;
/*!
 * \brief relative tolerance of the iterative solvers, the residuals being
 * checked up to 100 times this tolerance: the single precision of MFEM
 * bounds the tolerances
 */
static constexpr double tolerance =
    std::is_same_v<mfem_mgis::real, float> ? 1e-5 : 1e-10;
/*!
 * \brief bound of the residuals of the converged resolutions of the linear
 * elasticity: the GMRES solver of MFEM, preconditioned on the left, stops on
 * the preconditioned residual, which is much smaller than the residual for
 * the preconditioners of this problem (by a few hundred for the incomplete
 * factorizations). In single precision, the residuals of the converged
 * resolutions reach 1e-2 (`MINRESSolver`, `GMRESSolver` preconditioned by
 * `HypreEuclid`)
 */
static constexpr double elasticity_tolerance =
    std::is_same_v<mfem_mgis::real, float> ? 2e-2 : 1e-6;
/*!
 * \brief maximum number of iterations of the iterative solvers: 5000 for
 * the Laplace equation, 2000 for the linear elasticity (the slowest entry
 * converging on it needs 1500 iterations)
 */
static int maximum_number_of_iterations = 5000;
//! \brief if true, the computations are parallel
static bool parallel = false;
//! \brief rank of the process
static int rank = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \return true if the solver or preconditioner can compute in this run:
 * the entries which only compute on the GPU require MFEM computing on a
 * CUDA GPU (or on a HIP GPU for MAGMA, which supports both), the ones
 * which only compute on the host require hypre computing
 * on the host, except the ones of MFEM, SuiteSparse and oneMKL, which copy
 * the matrices and the vectors on the host
 * \param[in] d: description
 */
static bool isAvailable(const mfem_mgis::LinearSolverDescription& d) {
  using mfem_mgis::LinearSolverPlacement;
  if (d.placement == LinearSolverPlacement::DEVICE) {
    if (d.library == "MAGMA") {
      return mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                                  mfem::Backend::HIP_MASK);
    }
    return mfem::Device::Allows(mfem::Backend::CUDA_MASK);
  }
  if (d.placement == LinearSolverPlacement::HOST) {
    return (d.library == "MFEM") || (d.library == "SuiteSparse") ||
           (d.library == "oneMKL") || (!mfem_mgis::isHypreUsingGPU());
  }
  return true;
}  // end of isAvailable

//! \return if the given option is accepted
static bool accepts(const mfem_mgis::LinearSolverDescription& d,
                    const std::string& o) {
  return std::find(d.options.begin(), d.options.end(), o) != d.options.end();
}  // end of accepts

/*!
 * \return the options required by the entries which describe a whole
 * configuration of their library: a conjugate gradient and a Jacobi
 * preconditioner (wrapped in the Schwarz preconditioner in parallel) for
 * the generic entries of Ginkgo
 * \param[in] n: name of the entry
 */
template <bool is_parallel>
static mfem_mgis::Parameters getRequiredOptions(const std::string& n) {
  using namespace mfem_mgis;
  if (n == "GinkgoSolver") {
    return Parameters{{"Configuration", Parameters{{"type", "solver::Cg"}}}};
  }
  if (n == "GinkgoPreconditioner") {
    const auto jacobi = Parameters{{"type", "preconditioner::Jacobi"}};
    if constexpr (is_parallel) {
      return Parameters{
          {"Configuration", Parameters{{"type", "preconditioner::Schwarz"},
                                       {"local_solver", jacobi}}}};
    } else {
      return Parameters{{"Configuration", jacobi}};
    }
  }
  return Parameters{};
}  // end of getRequiredOptions

/*!
 * \return the parameters of a linear solver for this test: tolerances and
 * maximum number of iterations for the iterative solvers, and a
 * preconditioner for the solvers which accept one (a damped Jacobi
 * smoother, or a preconditioner of the same library), and the options
 * required by the entry (see `getRequiredOptions`)
 * \param[in] n: name of the linear solver
 * \param[in] d: description of the linear solver
 */
template <bool is_parallel>
static std::optional<mfem_mgis::Parameters> getTestParameters(
    const std::string& n, const mfem_mgis::LinearSolverDescription& d) {
  using namespace mfem_mgis;
  auto p = getRequiredOptions<is_parallel>(n);
  if (accepts(d, "RelativeTolerance")) {
    p.replaceOrInsert("RelativeTolerance", tolerance);
  }
  if (accepts(d, "Tolerance")) {
    p.replaceOrInsert("Tolerance", tolerance);
  }
  if (accepts(d, "AbsoluteTolerance")) {
    p.replaceOrInsert("AbsoluteTolerance", 0.);
  }
  if (accepts(d, "MaximumNumberOfIterations")) {
    p.replaceOrInsert("MaximumNumberOfIterations",
                      maximum_number_of_iterations);
  }
  if (d.preconditioners == AcceptedPreconditioners::ANY) {
    p.replaceOrInsert("Preconditioner",
                      Parameters{{"Name", "OperatorJacobiSmoother"},
                                 {"Options", Parameters{{"Damping", 0.5}}}});
  } else if (d.preconditioners == AcceptedPreconditioners::SAME_LIBRARY) {
    // a Jacobi preconditioner of the library if any (the simplest one, which
    // can be transposed), or its first available preconditioner
    const auto& f = PreconditionerFactory<is_parallel>::getFactory();
    auto preconditioner = std::optional<std::string>{};
    for (const auto& pn : f.getNames()) {
      const auto* const pd = f.getDescription(pn);
      if ((pd->library != d.library) || (!isAvailable(*pd))) {
        continue;
      }
      if (pn.find("Jacobi") != std::string::npos) {
        preconditioner = pn;
        break;
      }
      if (!preconditioner.has_value()) {
        preconditioner = pn;
      }
    }
    if (!preconditioner.has_value()) {
      return {};
    }
    p.replaceOrInsert("Preconditioner", Parameters{{"Name", *preconditioner}});
  }
  return p;
}  // end of getTestParameters

/*!
 * \return the options of a preconditioner for this test: a Jacobi smoother
 * for the wrappers of another preconditioner, an F-relaxation by algebraic
 * multigrid for the multigrid reduction of hypre (the Jacobi relaxation of
 * hypre, by default, is slower than the Jacobi preconditioner on the linear
 * elasticity), and the options required by the entry (see
 * `getRequiredOptions`)
 * \param[in] n: name of the preconditioner
 * \param[in] d: description of the preconditioner
 */
template <bool is_parallel>
static mfem_mgis::Parameters getTestOptions(
    const std::string& n, const mfem_mgis::LinearSolverDescription& d) {
  using namespace mfem_mgis;
  auto o = getRequiredOptions<is_parallel>(n);
  if (accepts(d, "Preconditioner")) {
    o.replaceOrInsert("Preconditioner",
                      Parameters{{"Name", "OperatorJacobiSmoother"}});
  }
  if (n == "HypreMGR") {
    o.replaceOrInsert("RelaxationType", 2);
  }
  return o;
}  // end of getTestOptions

//! \return the outcome of a resolution, in a short text
static std::string describe(const bool converged,
                            const std::optional<int> iterations,
                            const double residual) {
  auto s = std::string{};
  if (iterations.has_value()) {
    s += std::to_string(*iterations) + " iterations, ";
  }
  auto r = std::ostringstream{};
  r << std::scientific << std::setprecision(2) << residual;
  s += "residual " + r.str();
  if (!converged) {
    s += " (not converged)";
  }
  return s;
}  // end of describe

/*!
 * \brief checks the descriptions of the registered entries
 * \param[in] solvers: if true, the entries are linear solvers
 * \param[in] n: name
 * \param[in] d: description
 */
static void checkDescription(const bool solvers,
                             const std::string& n,
                             const mfem_mgis::LinearSolverDescription& d) {
  using namespace mfem_mgis;
  check(!d.library.empty(), n + ": no library");
  check(!d.description.empty(), n + ": no description");
  if (solvers) {
    check(d.kind != LinearSolverKind::PRECONDITIONER,
          n + ": linear solver described as a preconditioner");
    check(d.placement.has_value(), n + ": linear solver without placement");
    check((d.preconditioners == AcceptedPreconditioners::NONE) ||
              accepts(d, "Preconditioner"),
          n + ": preconditioners accepted without the option "
              "'Preconditioner'");
  } else {
    check(d.kind == LinearSolverKind::PRECONDITIONER,
          n + ": preconditioner not described as a preconditioner");
  }
}  // end of checkDescription

/*!
 * \brief checks the placements of the entries registered before the
 * descriptions, whose locations were deduced from their names
 * \param[in] n: name
 * \param[in] d: description
 */
static void checkPlacement(const std::string& n,
                           const mfem_mgis::LinearSolverDescription& d) {
  using mfem_mgis::LinearSolverPlacement;
  using Placement = std::optional<LinearSolverPlacement>;
  static const auto placements = std::map<std::string, Placement>{
      {"AmgX", LinearSolverPlacement::DEVICE},
      {"AmgXPCG", LinearSolverPlacement::DEVICE},
      {"AssembledOperator", Placement{}},
      {"BiCGSTABSolver", LinearSolverPlacement::MFEM},
      {"BlockILU", LinearSolverPlacement::HOST},
      {"CGSolver", LinearSolverPlacement::MFEM},
      {"CuDSSSolver", LinearSolverPlacement::DEVICE},
      {"DSmoother", LinearSolverPlacement::HOST},
      {"FGMRESSolver", LinearSolverPlacement::MFEM},
      {"GMRESSolver", LinearSolverPlacement::MFEM},
      {"GSSmoother", LinearSolverPlacement::HOST},
      {"HypreBoomerAMG", LinearSolverPlacement::HYPRE},
      {"HypreDiagScale", LinearSolverPlacement::HYPRE},
      {"HypreEuclid", LinearSolverPlacement::HOST},
      {"HypreFGMRES", LinearSolverPlacement::HYPRE},
      {"HypreGMRES", LinearSolverPlacement::HYPRE},
      {"HypreILU", LinearSolverPlacement::HYPRE},
      {"HypreParaSails", LinearSolverPlacement::HOST},
      {"HyprePCG", LinearSolverPlacement::HYPRE},
      {"KLU", LinearSolverPlacement::HOST},
      {"KLUSolver", LinearSolverPlacement::HOST},
      {"MAGMABiCGSTAB", LinearSolverPlacement::DEVICE},
      {"MAGMACG", LinearSolverPlacement::DEVICE},
      {"MAGMACGS", LinearSolverPlacement::DEVICE},
      {"MAGMAGMRES", LinearSolverPlacement::DEVICE},
      {"MAGMAIC", LinearSolverPlacement::DEVICE},
      {"MAGMAIDR", LinearSolverPlacement::DEVICE},
      {"MAGMAILU", LinearSolverPlacement::DEVICE},
      {"MAGMAJacobi", LinearSolverPlacement::DEVICE},
      {"MAGMAQMR", LinearSolverPlacement::DEVICE},
      {"MAGMATFQMR", LinearSolverPlacement::DEVICE},
      {"MKLClusterSparseSolver", LinearSolverPlacement::HOST},
      {"MKLPardisoSolver", LinearSolverPlacement::HOST},
      {"MINRESSolver", LinearSolverPlacement::MFEM},
      {"MUMPSSolver", LinearSolverPlacement::HOST},
      {"OperatorChebyshevSmoother", LinearSolverPlacement::MFEM},
      {"OperatorJacobiSmoother", LinearSolverPlacement::MFEM},
      {"OrthoSolver", LinearSolverPlacement::HOST},
      {"SLISolver", LinearSolverPlacement::MFEM},
      {"UMFPack", LinearSolverPlacement::HOST},
      {"UMFPackSolver", LinearSolverPlacement::HOST}};
  const auto p = placements.find(n);
  if (p != placements.end()) {
    check(d.placement == p->second, n + ": unexpected placement");
  }
}  // end of checkPlacement

/*!
 * \return true for the entries which, in single precision only, do not
 * converge on the linear elasticity with the options of this test: the
 * BiCGSTAB solver of MFEM in parallel (divergence), and the BiCGSTAB,
 * FGMRES and LGMRES solvers of hypre (stagnation above the tolerance)
 * \param[in] n: name of the entry
 */
static bool mayNotConvergeOnElasticityInSinglePrecision(const std::string& n) {
  static const auto names =
      std::set<std::string>{"HypreBiCGSTAB", "HypreFGMRES", "HypreLGMRES"};
  return names.contains(n) || ((n == "BiCGSTABSolver") && parallel);
}  // end of mayNotConvergeOnElasticityInSinglePrecision

/*!
 * \return true for the entries which do not converge on the linear
 * elasticity with the options of this test: only their report of the
 * convergence is checked. These are the conjugate gradient squared with a
 * Jacobi preconditioner (divergence for Ginkgo, breakdown for MAGMA), the
 * stationary iterations with a damped Jacobi preconditioner (`SLISolver`,
 * `GinkgoIR`: too slow), the pipelined conjugate gradient of Ginkgo on a
 * GPU (stagnation above the tolerance; it converges on the host), the
 * TFQMR of MAGMA (breakdown) and the multigrid of Ginkgo (aggregation
 * without the rigid body modes), and, in single precision, the entries of
 * `mayNotConvergeOnElasticityInSinglePrecision`.
 * \param[in] n: name of the entry
 */
static bool mayNotConvergeOnElasticity(const std::string& n) {
  static const auto names =
      std::set<std::string>{"GinkgoCGS", "MAGMACGS",   "SLISolver",
                            "GinkgoIR",  "MAGMATFQMR", "GinkgoMultigrid"};
  if (names.contains(n)) {
    return true;
  }
  if ((n == "GinkgoPipeCG") &&
      (mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                            mfem::Backend::HIP_MASK))) {
    return true;
  }
  return std::is_same_v<mfem_mgis::real, float> &&
         mayNotConvergeOnElasticityInSinglePrecision(n);
}  // end of mayNotConvergeOnElasticity

/*!
 * \brief checks the outcome of a resolution
 * \param[in] n: name of the resolution
 * \param[in] entry: name of the entry tested
 * \param[in] kind: kind of the entry tested (the direct solvers are
 * checked as on the Laplace equation)
 * \param[in] converged: convergence reported
 * \param[in] r: relative residual
 * \param[in] elasticity: if true, the system is the linear elasticity
 */
static void checkResolution(const std::string& n,
                            const std::string& entry,
                            const mfem_mgis::LinearSolverKind kind,
                            const bool converged,
                            const double r,
                            const bool elasticity) {
  if ((!elasticity) || (kind == mfem_mgis::LinearSolverKind::DIRECT)) {
    check(converged && (r < 100 * tolerance), n + ": resolution");
    return;
  }
  check((!converged) || (r < elasticity_tolerance),
        n + ": convergence reported with a large residual");
  check(converged || mayNotConvergeOnElasticity(entry), n + ": resolution");
}  // end of checkResolution

/*!
 * \brief checks all the registered linear solvers and preconditioners
 * \param[in, out] s: system (Laplace equation or linear elasticity)
 * \param[in] elasticity: if true, the system is the linear elasticity
 */
template <bool is_parallel, template <bool> class System>
static void run(System<is_parallel>& s, const bool elasticity) {
  using namespace mfem_mgis;
  const auto& lf = LinearSolverFactory<is_parallel>::getFactory();
  const auto& pf = PreconditionerFactory<is_parallel>::getFactory();
  auto print = [](const std::string& msg) {
    if (rank == 0) {
      std::cout << msg << '\n';
    }
  };
  // linear solvers
  for (const auto& n : lf.getNames()) {
    const auto* const d = lf.getDescription(n);
    check(d != nullptr, n + ": no description");
    if (d == nullptr) {
      continue;
    }
    checkDescription(true, n, *d);
    checkPlacement(n, *d);
    if (!isAvailable(*d)) {
      print(n + ": skipped (not available in this run)");
      continue;
    }
    const auto oparams = getTestParameters<is_parallel>(n, *d);
    if (!oparams.has_value()) {
      print(n + ": skipped (no preconditioner of " + d->library +
            " available)");
      continue;
    }
    auto ctx = Context{};
    auto ls = lf.generate(ctx, n, s.fespace, *oparams);
    if (isInvalid(ls)) {
      check(false, n + ": generation failed: " + ctx.getErrorMessage());
      continue;
    }
    auto x = mfem::Vector(s.B.Size());
    x.UseDevice(true);
    x = 0.;
    ls.linear_solver->SetOperator(s.A);
    ls.linear_solver->Mult(s.B, x);
    const auto r = s.getRelativeResidual(x);
    const auto converged = hasConverged(*(ls.linear_solver));
    const auto iterations =
        getNumberOfIterationsAtConvergence(*(ls.linear_solver));
    check((d->preconditioners == AcceptedPreconditioners::NONE) ||
              (ls.preconditioner != nullptr),
          n + ": the preconditioner is not built");
    print(n + ": " + describe(converged, iterations, r));
    if (d->kind == LinearSolverKind::APPROXIMATE) {
      // one cycle of a scalar algebraic multigrid (default options) does not
      // reduce the residual of the linear elasticity: the solution is only
      // checked to be finite and not null
      const auto xnorm = [&s, &x] {
        if constexpr (is_parallel) {
          // not `mfem::ParNormlp`, which also computes the norm from the
          // values of the host without synchronizing them (MFEM 4.10): a
          // solution written on the GPU had a null norm
          return std::sqrt(mfem::InnerProduct(s.fespace.GetComm(), x, x));
        } else {
          static_cast<void>(s);
          return x.Norml2();
        }
      }();
      check(elasticity ? std::isfinite(r) && (xnorm > 0) : r < 1,
            n + ": the residual does not decrease");
    } else {
      checkResolution(n, n, d->kind, converged, r, elasticity);
    }
  }
  // preconditioners, given to the GMRES solver of MFEM
  for (const auto& n : pf.getNames()) {
    const auto* const d = pf.getDescription(n);
    check(d != nullptr, n + ": no description");
    if (d == nullptr) {
      continue;
    }
    checkDescription(false, n, *d);
    checkPlacement(n, *d);
    if ((d->placement.has_value()) && (!isAvailable(*d))) {
      print(n + ": skipped (not available in this run)");
      continue;
    }
    // the orthogonalizing wrapper is meant for the operators whose kernel
    // is the constant vector, which is not the case of this system: it is
    // tested by MFEMSolversTest
    if (n == "OrthoSolver") {
      print(n + ": skipped (operators whose kernel is the constant vector)");
      continue;
    }
    // the multigrid reduction splits the unknowns by their components: it
    // is tested on the linear elasticity, and by HypreMGRTest
    if ((n == "HypreMGR") && (s.fespace.GetVDim() < 2)) {
      print(n + ": skipped (unknowns with one component)");
      continue;
    }
    auto ctx = Context{};
    auto ls = lf.generate(
        ctx, "GMRESSolver", s.fespace,
        Parameters{
            {"RelativeTolerance", tolerance},
            {"AbsoluteTolerance", 0.},
            {"MaximumNumberOfIterations", maximum_number_of_iterations},
            {"KDim", 50},
            {"Preconditioner",
             Parameters{{"Name", n},
                        {"Options", getTestOptions<is_parallel>(n, *d)}}}});
    if (isInvalid(ls)) {
      check(false, n + ": generation failed: " + ctx.getErrorMessage());
      continue;
    }
    auto x = mfem::Vector(s.B.Size());
    x.UseDevice(true);
    x = 0.;
    ls.linear_solver->SetOperator(s.A);
    ls.linear_solver->Mult(s.B, x);
    const auto r = s.getRelativeResidual(x);
    const auto converged = hasConverged(*(ls.linear_solver));
    print("GMRESSolver, " + n + ": " +
          describe(converged,
                   getNumberOfIterationsAtConvergence(*(ls.linear_solver)), r));
    check(ls.preconditioner != nullptr,
          "GMRESSolver, " + n + ": the preconditioner is not built");
    checkResolution("GMRESSolver, " + n, n, LinearSolverKind::ITERATIVE,
                    converged, r, elasticity);
  }
  // preconditioners refused by the descriptions: a preconditioner given to
  // a linear solver which accepts none, and a preconditioner of another
  // library given to a linear solver which only accepts the ones of its
  // library (the Jacobi smoother of MFEM is always available)
  for (const auto& n : lf.getNames()) {
    const auto* const d = lf.getDescription(n);
    if ((d == nullptr) || (!isAvailable(*d)) ||
        (d->preconditioners == AcceptedPreconditioners::ANY) ||
        (d->library == "MFEM")) {
      continue;
    }
    const auto message =
        d->preconditioners == AcceptedPreconditioners::NONE
            ? std::string{"accepts no preconditioner"}
            : "only accepts the preconditioners of " + d->library;
    auto ctx = Context{};
    auto ls = lf.generate(
        ctx, n, s.fespace,
        Parameters{{"Preconditioner",
                    Parameters{{"Name", "OperatorJacobiSmoother"}}}});
    check(isInvalid(ls) &&
              (ctx.getErrorMessage().find(message) != std::string::npos),
          n + ": preconditioner 'OperatorJacobiSmoother' not refused");
  }
  // in sequential, the solvers and the preconditioners only available in
  // parallel are reported as such
  if constexpr (!is_parallel) {
#ifdef MFEM_USE_MPI
    auto ctx = Context{};
    auto ls = lf.generate(ctx, "HyprePCG", s.fespace, Parameters{});
    check(isInvalid(ls) &&
              (ctx.getErrorMessage().find("only available in parallel") !=
               std::string::npos),
          "HyprePCG: not reported as only available in parallel");
    auto ctx2 = Context{};
    auto p = pf.generate(ctx2, "HypreBoomerAMG", s.fespace, Parameters{});
    check((p == nullptr) &&
              (ctx2.getErrorMessage().find("only available in parallel") !=
               std::string::npos),
          "HypreBoomerAMG: not reported as only available in parallel");
#endif /* MFEM_USE_MPI */
  }
}  // end of run

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "parallel computations");
  auto elasticity = false;
  args.AddOption(&elasticity, "-e", "--elasticity", "-no-e", "--no-elasticity",
                 "linear elasticity instead of the Laplace equation");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
#ifdef MFEM_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif /* MFEM_USE_MPI */
  if (elasticity) {
    maximum_number_of_iterations = 2000;
  }
  // a cube, or a beam for the linear elasticity
  auto mesh = elasticity ? mfem::Mesh::MakeCartesian3D(
                               12, 3, 3, mfem::Element::HEXAHEDRON, 4., 1., 1.)
                         : mfem::Mesh::MakeCartesian3D(
                               8, 8, 8, mfem::Element::HEXAHEDRON);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    const auto test = [&mesh, &fec, elasticity]<bool is_parallel>() {
      if (elasticity) {
        auto s = ElasticitySystem<is_parallel>(mesh, fec);
        run<is_parallel>(s, elasticity);
      } else {
        auto s = LaplaceSystem<is_parallel>(mesh, fec);
        run<is_parallel>(s, elasticity);
      }
    };
    if (parallel) {
#ifdef MFEM_USE_MPI
      test.template operator()<true>();
#else  /* MFEM_USE_MPI */
      throw std::runtime_error("parallel computations are not supported");
#endif /* MFEM_USE_MPI */
    } else {
      test.template operator()<false>();
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
