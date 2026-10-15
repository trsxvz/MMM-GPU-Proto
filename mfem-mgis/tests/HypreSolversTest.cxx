/*!
 * \file   HypreSolversTest.cxx
 * \brief  Test of the solvers and of the preconditioners of hypre which are
 * not wrapped by MFEM (BiCGSTAB, COGMRES, LGMRES, hybrid solver,
 * factorized sparse approximate inverse), of the smoothers of hypre
 * (`HypreSmoother`) and of the options of BoomerAMG, on the Laplace
 * equation on a cube, in parallel (the generic test
 * `LinearSolverFactoriesTest` solves it with each entry and its default
 * options):
 *
 * - the Krylov solvers converge with BoomerAMG, report their number of
 *   iterations, give access to their relative tolerance, report a
 *   resolution with a null right hand side as converged and a resolution
 *   stopped by the maximum number of iterations as not converged;
 * - the hybrid solver converges, with each of its Krylov solvers, and
 *   switches to BoomerAMG when the convergence of the diagonal scaling is
 *   too slow;
 * - the factorized sparse approximate inverse reduces the number of
 *   iterations compared to the diagonal scaling on one process;
 * - each type of smoother of hypre (but the lumped Jacobi smoother, meant
 *   for mass matrices) preconditions the GMRES solver of MFEM;
 * - a second operator, of the same pattern, is given to the Krylov solvers
 *   (with the factorized sparse approximate inverse, built again), to the
 *   hybrid solver, and to BoomerAMG, whose options are applied again;
 * - the options of BoomerAMG are accepted, and invalid values refused.
 * \date   15/10/2026
 */

#include <cmath>
#include <memory>
#include <string>
#include <string_view>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
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
static mfem_mgis::LinearSolverHandler generate(
    LaplaceSystem<true>& s,
    const std::string& name,
    const mfem_mgis::Parameters& parameters) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto ls = LinearSolverFactory<true>::getFactory().generate(
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
static std::string getGenerationError(LaplaceSystem<true>& s,
                                      const std::string& name,
                                      const mfem_mgis::Parameters& parameters) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  const auto ls = LinearSolverFactory<true>::getFactory().generate(
      ctx, name, s.fespace, parameters);
  return isInvalid(ls) ? ctx.getErrorMessage() : "";
}  // end of getGenerationError

//! \return the resolution of the system by a linear solver
static Resolution solve(LaplaceSystem<true>& s,
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

/*!
 * \return the parameters of a Krylov solver of hypre
 * \param[in] p: preconditioner
 * \param[in] o: options of the preconditioner
 */
static mfem_mgis::Parameters getKrylovParameters(
    const std::string& p, const mfem_mgis::Parameters& o = {}) {
  using namespace mfem_mgis;
  return Parameters{
      {"Tolerance", tolerance},
      {"MaximumNumberOfIterations", 1000},
      {"Preconditioner", Parameters{{"Name", p}, {"Options", o}}}};
}  // end of getKrylovParameters

/*!
 * \return the parameters of the GMRES solver of MFEM
 * \param[in] p: preconditioner
 * \param[in] o: options of the preconditioner
 */
static mfem_mgis::Parameters getGMRESParameters(
    const std::string& p, const mfem_mgis::Parameters& o = {}) {
  using namespace mfem_mgis;
  return Parameters{
      {"RelativeTolerance", tolerance},
      {"MaximumNumberOfIterations", 2000},
      {"KDim", 100},
      {"Preconditioner", Parameters{{"Name", p}, {"Options", o}}}};
}  // end of getGMRESParameters

/*!
 * \return the resolution of the system with another operator: the matrix of
 * the system times two, whose solution is the half of the one of the system
 * (a new operator of the same pattern, as in the iterations of a Newton
 * solver)
 * \param[in] s: system
 * \param[in] a2: matrix of the system times two
 * \param[in] ls: linear solver
 * \param[in] name: name of the resolution
 */
static Resolution solveDoubled(LaplaceSystem<true>& s,
                               const mfem::HypreParMatrix& a2,
                               mfem::Solver& ls,
                               const std::string& name) {
  using namespace mfem_mgis;
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.SetOperator(a2);
  ls.Mult(s.B, x);
  auto res = Resolution{};
  res.converged = hasConverged(ls);
  res.iterations = getNumberOfIterationsAtConvergence(ls).value_or(-1);
  res.residual = s.getRelativeResidual(a2, x);
  if (rank == 0) {
    std::cout << name << ": " << res.iterations << " iterations, residual "
              << res.residual << (res.converged ? "" : " (not converged)")
              << '\n';
  }
  return res;
}  // end of solveDoubled

static void run(LaplaceSystem<true>& s) {
  using namespace mfem_mgis;
  const auto amg = Parameters{{"VerbosityLevel", 0}};
  // Krylov solvers
  // COGMRES only computes on the host
  const auto cogmres = !isHypreUsingGPU();
  if (!cogmres) {
    check(
        getGenerationError(s, "HypreCOGMRES",
                           getKrylovParameters("HypreBoomerAMG", amg))
                .find("on the host, which hypre does not") != std::string::npos,
        "HypreCOGMRES: not refused when hypre computes on a GPU");
  }
  for (const auto& n : {"HypreBiCGSTAB", "HypreCOGMRES", "HypreLGMRES"}) {
    if ((!cogmres) && (std::string_view{n} == "HypreCOGMRES")) {
      continue;
    }
    auto ls = generate(s, n, getKrylovParameters("HypreBoomerAMG", amg));
    check(isSolved(solve(s, *(ls.linear_solver), std::string{n} + ", AMG")),
          std::string{n} + ": resolution");
    // the tolerance given is bounded below in single precision
    const auto expected =
        static_cast<double>(boundRelativeTolerance(tolerance));
    const auto otol = getRelativeTolerance(*(ls.linear_solver));
    check(otol.has_value() && (std::abs(static_cast<double>(*otol) -
                                        expected) <= 1e-3 * expected),
          std::string{n} + ": relative tolerance not read");
    check(setRelativeTolerance(*(ls.linear_solver), 1e-3) &&
              (getRelativeTolerance(*(ls.linear_solver)) ==
               static_cast<real>(1e-3)),
          std::string{n} + ": relative tolerance not set");
  }
  for (const auto& n : {"HypreCOGMRES", "HypreLGMRES"}) {
    if ((!cogmres) && (std::string_view{n} == "HypreCOGMRES")) {
      continue;
    }
    auto p = getKrylovParameters("HypreDiagScale");
    p.replaceOrInsert("KDim", 5);
    auto ls = generate(s, n, p);
    check(isSolved(solve(s, *(ls.linear_solver), std::string{n} + ", KDim 5")),
          std::string{n} + ": resolution with a Krylov space of dimension 5");
  }
  // a null right hand side: the initial residual is null
  for (const auto& n : {"HypreBiCGSTAB", "HypreCOGMRES", "HypreLGMRES"}) {
    if ((!cogmres) && (std::string_view{n} == "HypreCOGMRES")) {
      continue;
    }
    auto ls = generate(s, n, getKrylovParameters("HypreDiagScale"));
    auto x = mfem::Vector(s.B.Size());
    auto zero = mfem::Vector(s.B.Size());
    x.UseDevice(true);
    zero.UseDevice(true);
    x = 1.;
    zero = 0.;
    ls.linear_solver->SetOperator(s.A);
    ls.linear_solver->Mult(zero, x);
    check(hasConverged(*(ls.linear_solver)) && (x.Normlinf() <= 0),
          std::string{n} + ": resolution with a null right hand side");
  }
  auto stopped = generate(
      s, "HypreBiCGSTAB",
      Parameters{{"Tolerance", tolerance}, {"MaximumNumberOfIterations", 2}});
  const auto rs =
      solve(s, *(stopped.linear_solver), "HypreBiCGSTAB, 2 iterations");
  check((!rs.converged) && (rs.iterations == 2),
        "resolution stopped by the maximum number of iterations");
  // a second operator of the same pattern: the preconditioners and the
  // hybrid solver of hypre are built again, and given to the Krylov solvers
  // the matrix is copied and scaled: on a GPU, `mfem::Add` sorts the
  // columns of each row, while the diagonal scaling of hypre reads the
  // diagonal as the first entry of each row
  auto a2 = std::make_unique<mfem::HypreParMatrix>(s.A);
  *a2 *= 2.;
  for (const auto& n : {"HypreBiCGSTAB", "HypreLGMRES", "HyprePCG"}) {
    auto ls = generate(s, n, getKrylovParameters("HypreFSAI"));
    const auto name = std::string{n} + ", HypreFSAI";
    check(isSolved(solve(s, *(ls.linear_solver), name)) &&
              isSolved(solveDoubled(s, *a2, *(ls.linear_solver),
                                    name + ", second operator")),
          name + ": resolution with a second operator");
  }
  {
    auto hybrid = generate(s, "HypreHybrid",
                           Parameters{{"Tolerance", tolerance},
                                      {"MaximumNumberOfIterations", 1000}});
    check(isSolved(solve(s, *(hybrid.linear_solver), "HypreHybrid")) &&
              isSolved(solveDoubled(s, *a2, *(hybrid.linear_solver),
                                    "HypreHybrid, second operator")),
          "HypreHybrid: resolution with a second operator");
  }
  // hybrid solver
  for (const auto& t : {"PCG", "GMRES", "BiCGSTAB"}) {
    auto hybrid = generate(s, "HypreHybrid",
                           Parameters{{"Tolerance", tolerance},
                                      {"MaximumNumberOfIterations", 1000},
                                      {"SolverType", t}});
    check(isSolved(solve(s, *(hybrid.linear_solver),
                         std::string{"HypreHybrid, "} + t)),
          std::string{"HypreHybrid, "} + t + ": resolution");
  }
  // a resolution stopped by the maximum number of iterations
  auto hybrid_stopped = generate(
      s, "HypreHybrid",
      Parameters{{"Tolerance", tolerance}, {"MaximumNumberOfIterations", 2}});
  const auto rhs =
      solve(s, *(hybrid_stopped.linear_solver), "HypreHybrid, 2 iterations");
  check(!rhs.converged,
        "HypreHybrid: resolution stopped by the maximum number of iterations "
        "reported as converged");
  // hypre switches to BoomerAMG when the convergence factor of the diagonal
  // scaling exceeds the given rate: a rate close to zero switches at once
  auto switched = generate(s, "HypreHybrid",
                           Parameters{{"Tolerance", tolerance},
                                      {"MaximumNumberOfIterations", 1000},
                                      {"ConvergenceRate", 1e-3}});
  auto ds = generate(s, "HyprePCG", getKrylovParameters("HypreDiagScale"));
  const auto rds = solve(s, *(ds.linear_solver), "HyprePCG, HypreDiagScale");
  const auto rsw =
      solve(s, *(switched.linear_solver), "HypreHybrid, switched to AMG");
  check(isSolved(rsw) && isFaster(rsw, rds),
        "HypreHybrid: no switch to BoomerAMG");
  // factorized sparse approximate inverse
  // the factorized approximate inverse of hypre is computed on the
  // diagonal block of each process: in parallel, it does not reduce the
  // number of iterations of the diagonal scaling on this small mesh, and
  // only the convergence is checked
  auto fsai = generate(s, "HyprePCG", getKrylovParameters("HypreFSAI"));
  const auto rf = solve(s, *(fsai.linear_solver), "HyprePCG, HypreFSAI");
  auto size = int{};
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  check(isSolved(rf) && ((size > 1) || (rf.iterations < rds.iterations)),
        "HyprePCG, HypreFSAI: no reduction of the number of iterations");
  auto mfem_fsai = generate(
      s, "CGSolver",
      Parameters{
          {"RelativeTolerance", tolerance},
          {"MaximumNumberOfIterations", 1000},
          {"Preconditioner",
           Parameters{{"Name", "HypreFSAI"},
                      {"Options", Parameters{{"MaximumSteps", 10},
                                             {"MaximumStepSize", 2}}}}}});
  check(isSolved(solve(s, *(mfem_fsai.linear_solver), "CGSolver, HypreFSAI")),
        "CGSolver, HypreFSAI: resolution");
  // smoothers of hypre (the lumped Jacobi smoother, which divides by the
  // sums of the rows, does not apply to a stiffness matrix)
  // the on-processor forward solve, the Taubin and FIR smoothers only
  // compute on the host
  for (const auto& t : {"Jacobi", "l1Jacobi", "l1GS", "l1GStr", "GS", "OPFS",
                        "Chebyshev", "Taubin", "FIR"}) {
    const auto host_only = (std::string_view{t} == "OPFS") ||
                           (std::string_view{t} == "Taubin") ||
                           (std::string_view{t} == "FIR");
    if ((host_only) && (isHypreUsingGPU())) {
      check(getGenerationError(
                s, "GMRESSolver",
                getGMRESParameters("HypreSmoother", Parameters{{"Type", t}}))
                    .find("on the host, which hypre does not") !=
                std::string::npos,
            std::string{"HypreSmoother: "} + t +
                " not refused when hypre computes on a GPU");
      continue;
    }
    auto smoother =
        generate(s, "GMRESSolver",
                 getGMRESParameters("HypreSmoother", Parameters{{"Type", t}}));
    check(isSolved(solve(s, *(smoother.linear_solver),
                         std::string{"GMRESSolver, HypreSmoother "} + t)),
          std::string{"GMRESSolver, HypreSmoother "} + t + ": resolution");
  }
  // options of BoomerAMG
  auto options = amg;
  options.replaceOrInsert("TruncationFactor", 0.1);
  options.replaceOrInsert("MaximumElementsPerRow", 4);
  options.replaceOrInsert("AggressiveCoarseningLevels", 1);
  options.replaceOrInsert("AggressiveInterpolationType", 4);
  options.replaceOrInsert("NumberOfPaths", 2);
  options.replaceOrInsert("MaximumCoarseSize", 20);
  options.replaceOrInsert("MinimumCoarseSize", 2);
  options.replaceOrInsert("RelaxationWeight", 0.9);
  options.replaceOrInsert("KeepTranspose", true);
  auto tuned =
      generate(s, "HyprePCG", getKrylovParameters("HypreBoomerAMG", options));
  check(isSolved(solve(s, *(tuned.linear_solver), "HyprePCG, tuned AMG")) &&
            isSolved(solveDoubled(s, *a2, *(tuned.linear_solver),
                                  "HyprePCG, tuned AMG, second operator")),
        "HyprePCG, BoomerAMG with options: resolution");
  // the options given to hypre are applied again for each operator (MFEM
  // creates a new BoomerAMG of hypre)
  {
    const auto* const p =
        dynamic_cast<const mfem::HypreBoomerAMG*>(tuned.preconditioner.get());
    auto truncation = HYPRE_Real{};
    auto elements = HYPRE_Int{};
    auto maximum_coarse_size = HYPRE_Int{};
    auto minimum_coarse_size = HYPRE_Int{};
    if (p != nullptr) {
      HYPRE_BoomerAMGGetTruncFactor(*p, &truncation);
      HYPRE_BoomerAMGGetPMaxElmts(*p, &elements);
      HYPRE_BoomerAMGGetMaxCoarseSize(*p, &maximum_coarse_size);
      HYPRE_BoomerAMGGetMinCoarseSize(*p, &minimum_coarse_size);
    }
    check((p != nullptr) &&
              (std::abs(static_cast<double>(truncation) - 0.1) < 1e-6) &&
              (elements == 4) && (maximum_coarse_size == 20) &&
              (minimum_coarse_size == 2),
          "HypreBoomerAMG: options not applied to the second operator");
  }
  auto chebyshev = amg;
  chebyshev.replaceOrInsert("RelaxationType", 16);
  chebyshev.replaceOrInsert("ChebyshevOrder", 3);
  auto rc =
      generate(s, "HyprePCG", getKrylovParameters("HypreBoomerAMG", chebyshev));
  check(isSolved(solve(s, *(rc.linear_solver), "HyprePCG, AMG, Chebyshev")),
        "HyprePCG, BoomerAMG with Chebyshev smoothers: resolution");
  // invalid values
  const auto refused = [&s](const std::string& n, const Parameters& p,
                            const std::string& msg) {
    return getGenerationError(s, n, p).find(msg) != std::string::npos;
  };
  check(
      refused("GMRESSolver",
              getGMRESParameters("HypreSmoother", Parameters{{"Type", "SOR"}}),
              "invalid type"),
      "invalid type of smoother not refused");
  check(refused("HyprePCG", getKrylovParameters("HypreSmoother"),
                "applied by MFEM"),
        "smoother of hypre given to a Krylov solver of hypre not refused");
  check(refused("HyprePCG",
                getKrylovParameters("HypreBoomerAMG",
                                    Parameters{{"TruncationFactor", -1.}}),
                "invalid value"),
        "invalid truncation factor not refused");
  check(
      refused("HypreHybrid", Parameters{{"SolverType", "CG"}}, "invalid type"),
      "invalid Krylov solver of the hybrid solver not refused");
  check(refused("HypreHybrid", Parameters{{"KDim", 0}}, "invalid value"),
        "invalid Krylov dimension of the hybrid solver not refused");
  check(refused("HypreHybrid", Parameters{{"Tolerance", 2.}},
                "invalid tolerance"),
        "invalid tolerance of the hybrid solver not refused");
  check(refused(
            "HyprePCG",
            getKrylovParameters("HypreFSAI", Parameters{{"AlgorithmType", 0}}),
            "invalid value"),
        "invalid algorithm of the factorized approximate inverse not refused");
}  // end of run

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
  auto mesh =
      mfem::Mesh::MakeCartesian3D(12, 12, 12, mfem::Element::HEXAHEDRON);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    auto s = LaplaceSystem<true>(mesh, fec);
    run(s);
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (rank == 0) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
