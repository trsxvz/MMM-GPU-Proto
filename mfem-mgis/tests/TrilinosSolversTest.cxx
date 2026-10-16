/*!
 * \file   TrilinosSolversTest.cxx
 * \brief  Test of the solvers and of the preconditioners of Trilinos built
 * by the factories, in parallel (the generic test
 * `LinearSolverFactoriesTest` solves the Laplace equation and the linear
 * elasticity with their default options):
 *
 * - the direct solver of Amesos2 (KLU2) solves the Laplace equation and the
 *   linear elasticity, the unknowns being ordered by nodes and by vector
 *   dimension, up to the round-off errors: the matrices of hypre are
 *   converted exactly, the unknowns ordered by nodes being permuted; Tacho
 *   solves the linear elasticity, if Trilinos is built with it;
 * - the Krylov solvers of Belos (CG, GMRES) converge, preconditioned by
 *   Ifpack2 (Jacobi relaxation, Chebyshev polynomial, incomplete
 *   factorization), and the algebraic multigrid of MueLu and the two level
 *   Schwarz preconditioner of FROSch reduce the number of iterations
 *   compared to the Jacobi relaxation, MueLu also with its factories
 *   written with Kokkos (`use kokkos refactor`);
 * - one V-cycle of MueLu, used as an approximate linear solver, reduces
 *   the residual, and takes the initial guess into account in the
 *   iterative mode (the stationary iteration converges);
 * - on the linear elasticity, MueLu needs close numbers of iterations for
 *   both orderings of the unknowns, and fewer iterations with the rigid
 *   body modes than with an aggregation of the unknowns one by one, and
 *   FROSch converges for both orderings;
 * - the fast incomplete factorizations of Ifpack2 (FAST_ILU, FAST_IC),
 *   applied to the diagonal block of each process, need close numbers of
 *   iterations to the incomplete factorization RILUK, if Trilinos is built
 *   with ShyLU FastILU;
 * - the block preconditioners of Teko precondition the linear elasticity,
 *   if mfem-mgis uses Teko: the block Jacobi preconditioner with MueLu for
 *   both orderings of the unknowns, and the block Gauss-Seidel
 *   preconditioner with Amesos2 on blocks of two and one components;
 *   invalid definitions of the inverses are refused with an explanation,
 *   and a process owning no unknown is refused;
 * - the preconditioners of Trilinos precondition the Krylov solvers of
 *   MFEM, and the preconditioners of mfem-mgis the Krylov solvers of Belos;
 * - the options of mfem-mgis set the tolerance and the maximum number of
 *   iterations, the parameters of Belos overriding them, and the relative
 *   tolerance can be read and changed (see `getRelativeTolerance` and
 *   `setRelativeTolerance`); a resolution stopped by the maximum number of
 *   iterations is reported as not converged;
 * - a second operator, the first one times two, is solved by the same
 *   solver, and a system with a process owning no unknown is solved
 *   (refused by FROSch and Teko);
 * - invalid parameters are refused.
 * \date   15/10/2026
 */

#include <cmath>
#include <tuple>
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "LaplaceSystem.hxx"
#include "ElasticitySystem.hxx"

using mfem_mgis::tests::ElasticitySystem;
using mfem_mgis::tests::LaplaceSystem;

//! \brief number of failed checks
static int failures = 0;
//! \brief rank of the process
static int rank = 0;

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
 * \return true if a resolution converged with a relative residual lower
 * than the given bound
 * \param[in] r: resolution
 * \param[in] bound: bound of the relative residual
 */
static bool isSolved(const Resolution& r, const double bound = 1e-8) {
  return r.converged && (r.residual < bound);
}  // end of isSolved

/*!
 * \return true if a preconditioner reduces the number of iterations to
 * less than two thirds of the one of another preconditioner
 * \param[in] r: resolution with the preconditioner
 * \param[in] r2: resolution with the other preconditioner
 */
static bool isFaster(const Resolution& r, const Resolution& r2) {
  return (r.iterations > 0) && (3 * r.iterations < 2 * r2.iterations);
}  // end of isFaster

/*!
 * \return the linear solver of the factory, throws on failure
 * \param[in] s: system
 * \param[in] name: name of the linear solver
 * \param[in] parameters: parameters of the linear solver
 */
template <typename System>
static mfem_mgis::LinearSolverHandler generate(
    System& s,
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
template <typename System>
static std::string getGenerationError(System& s,
                                      const std::string& name,
                                      const mfem_mgis::Parameters& parameters) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  const auto ls = LinearSolverFactory<true>::getFactory().generate(
      ctx, name, s.fespace, parameters);
  return isInvalid(ls) ? ctx.getErrorMessage() : "";
}  // end of getGenerationError

/*!
 * \return the resolution of a system by a linear solver
 * \param[in] s: system
 * \param[in] m: matrix of the system
 * \param[in] ls: linear solver
 * \param[in] name: name of the resolution, printed with its outcome
 */
template <typename System>
static Resolution solve(System& s,
                        const mfem::HypreParMatrix& m,
                        mfem::Solver& ls,
                        const std::string& name) {
  using namespace mfem_mgis;
  // vectors on the device, if any (the right hand side of the systems is
  // not), so that the solvers of Trilinos receive the vectors of the GPU
  auto b = mfem::Vector(s.B.Size());
  b.UseDevice(true);
  b = s.B;
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.SetOperator(m);
  ls.Mult(b, x);
  auto r = Resolution{};
  r.converged = hasConverged(ls);
  r.iterations = getNumberOfIterationsAtConvergence(ls).value_or(0);
  r.residual = s.getRelativeResidual(m, x);
  if (rank == 0) {
    std::cout << name << ": " << r.iterations << " iterations, residual "
              << r.residual << (r.converged ? "" : " (not converged)") << '\n';
  }
  return r;
}  // end of solve

//! \return the resolution of a system by a linear solver
template <typename System>
static Resolution solve(System& s, mfem::Solver& ls, const std::string& name) {
  return solve(s, s.A, ls, name);
}  // end of solve

/*!
 * \return the relative residual after some iterations of a linear solver M
 * in the iterative mode, x <- x + M (b - A x) if the linear solver takes the
 * initial guess into account
 * \param[in] s: system
 * \param[in] ls: linear solver
 * \param[in] n: number of iterations
 * \param[in] name: name of the resolution, printed with its outcome
 */
template <typename System>
static double iterate(System& s,
                      mfem::Solver& ls,
                      const int n,
                      const std::string& name) {
  auto b = mfem::Vector(s.B.Size());
  b.UseDevice(true);
  b = s.B;
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.SetOperator(s.A);
  ls.iterative_mode = true;
  for (int i = 0; i != n; ++i) {
    ls.Mult(b, x);
  }
  ls.iterative_mode = false;
  const auto r = s.getRelativeResidual(s.A, x);
  if (rank == 0) {
    std::cout << name << ": residual " << r << " after " << n
              << " iterations\n";
  }
  return r;
}  // end of iterate

/*!
 * \return the parameters of a Krylov solver of Belos
 * \param[in] solver: name of the Krylov solver of Belos
 * \param[in] preconditioner: name of the preconditioner
 * \param[in] options: options of the preconditioner
 */
static mfem_mgis::Parameters getBelosParameters(
    const std::string& solver,
    const std::string& preconditioner,
    const mfem_mgis::Parameters& options = {}) {
  using namespace mfem_mgis;
  return Parameters{{"Solver", solver},
                    {"RelativeTolerance", 1e-10},
                    {"MaximumNumberOfIterations", 2000},
                    {"Preconditioner", Parameters{{"Name", preconditioner},
                                                  {"Options", options}}}};
}  // end of getBelosParameters

/*!
 * \return the options of MueLu on the small systems of this test: MueLu
 * solves directly the systems smaller than its parameter `coarse: max size`
 * (2000 by default), which is lowered so that the multigrid has several
 * levels
 * \param[in] p: other parameters of MueLu
 */
static mfem_mgis::Parameters getMueLuOptions(mfem_mgis::Parameters p = {}) {
  p.insert(mfem_mgis::may_throw, "coarse: max size", 20);
  return mfem_mgis::Parameters{{"Parameters", p}};
}  // end of getMueLuOptions

/*!
 * \brief checks the direct solvers of Amesos2
 * \param[in] s: Laplace equation
 * \param[in] e1: linear elasticity, unknowns ordered by nodes
 * \param[in] e2: linear elasticity, unknowns ordered by vector dimension
 */
static void checkDirectSolvers(LaplaceSystem<true>& s,
                               ElasticitySystem<true>& e1,
                               ElasticitySystem<true>& e2) {
  using namespace mfem_mgis;
  auto ls = generate(s, "Amesos2Solver", {});
  check(
      isSolved(solve(s, *(ls.linear_solver), "Amesos2Solver, Laplace"), 1e-10),
      "Amesos2Solver: resolution of the Laplace equation");
  auto ls1 = generate(e1, "Amesos2Solver", {});
  check(isSolved(solve(e1, *(ls1.linear_solver),
                       "Amesos2Solver, elasticity, byNODES"),
                 1e-10),
        "Amesos2Solver: resolution of the linear elasticity (byNODES)");
  auto ls2 = generate(e2, "Amesos2Solver", Parameters{{"Solver", "KLU2"}});
  check(isSolved(solve(e2, *(ls2.linear_solver),
                       "Amesos2Solver, elasticity, byVDIM"),
                 1e-10),
        "Amesos2Solver: resolution of the linear elasticity (byVDIM)");
  // Tacho (Cholesky factorization, on the GPU when Tpetra computes there),
  // if Trilinos is built with it
  const auto tacho = Parameters{{"Solver", "Tacho"}};
  if (getGenerationError(e1, "Amesos2Solver", tacho).empty()) {
    auto ls3 = generate(e1, "Amesos2Solver", tacho);
    check(isSolved(solve(e1, *(ls3.linear_solver),
                         "Amesos2Solver, Tacho, elasticity, byNODES"),
                   1e-10),
          "Amesos2Solver: resolution of the linear elasticity by Tacho");
  } else if (rank == 0) {
    std::cout << "Amesos2Solver, Tacho: skipped (not available)\n";
  }
  // the factorization as preconditioner of the conjugate gradient of MFEM
  auto cg = generate(
      s, "CGSolver",
      Parameters{
          {"RelativeTolerance", 1e-10},
          {"MaximumNumberOfIterations", 10},
          {"Preconditioner", Parameters{{"Name", "Amesos2Preconditioner"}}}});
  const auto r =
      solve(s, *(cg.linear_solver), "CGSolver, Amesos2Preconditioner");
  check(isSolved(r) && (r.iterations <= 2),
        "CGSolver, Amesos2Preconditioner: resolution");
}  // end of checkDirectSolvers

/*!
 * \brief checks the Krylov solvers of Belos and the preconditioners of
 * Trilinos on the Laplace equation
 * \param[in] s: system
 */
static void checkKrylovSolvers(LaplaceSystem<true>& s) {
  using namespace mfem_mgis;
  auto jacobi = generate(s, "BelosSolver", getBelosParameters("CG", "Ifpack2"));
  const auto rj = solve(s, *(jacobi.linear_solver), "Belos CG, Ifpack2 Jacobi");
  check(isSolved(rj), "Belos CG, Ifpack2 Jacobi: resolution");
  auto chebyshev = generate(
      s, "BelosSolver",
      getBelosParameters("CG", "Ifpack2", Parameters{{"Type", "CHEBYSHEV"}}));
  check(isSolved(solve(s, *(chebyshev.linear_solver),
                       "Belos CG, Ifpack2 Chebyshev")),
        "Belos CG, Ifpack2 Chebyshev: resolution");
  auto ilu = generate(
      s, "BelosSolver",
      getBelosParameters("GMRES", "Ifpack2", Parameters{{"Type", "RILUK"}}));
  check(isSolved(solve(s, *(ilu.linear_solver), "Belos GMRES, Ifpack2 RILUK")),
        "Belos GMRES, Ifpack2 RILUK: resolution");
  auto muelu = generate(s, "BelosSolver", getBelosParameters("CG", "MueLu"));
  const auto rm = solve(s, *(muelu.linear_solver), "Belos CG, MueLu");
  check(isSolved(rm) && isFaster(rm, rj), "Belos CG, MueLu: resolution");
  // the factories of MueLu written with Kokkos (the default of MueLu on all
  // the nodes but the serial one, not of mfem-mgis)
  auto kokkos =
      generate(s, "BelosSolver",
               getBelosParameters(
                   "CG", "MueLu",
                   Parameters{{"Parameters",
                               Parameters{{"use kokkos refactor", true}}}}));
  const auto rk = solve(s, *(kokkos.linear_solver),
                        "Belos CG, MueLu, factories written with Kokkos");
  check(isSolved(rk) && isFaster(rk, rj),
        "Belos CG, MueLu with the factories written with Kokkos: resolution");
  auto frosch = generate(s, "BelosSolver", getBelosParameters("CG", "FROSch"));
  const auto rf = solve(s, *(frosch.linear_solver), "Belos CG, FROSch");
  // the iterations of the two level Schwarz method grow with the number of
  // subdomains on this small system (20 against 25 for the Jacobi
  // relaxation on three processes). With one process, the only subdomain
  // is factorized.
  auto size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  check(isSolved(rf) && (rf.iterations < rj.iterations) &&
            ((size != 1) || (rf.iterations <= 2)),
        "Belos CG, FROSch: resolution");
  // a preconditioner of Trilinos with a Krylov solver of MFEM
  auto cg =
      generate(s, "CGSolver",
               Parameters{{"RelativeTolerance", 1e-10},
                          {"MaximumNumberOfIterations", 2000},
                          {"Preconditioner", Parameters{{"Name", "MueLu"}}}});
  const auto rcg = solve(s, *(cg.linear_solver), "CGSolver, MueLu");
  check(isSolved(rcg) && isFaster(rcg, rj), "CGSolver, MueLu: resolution");
  // one V-cycle of MueLu as an approximate linear solver, which reports no
  // iterations
  auto vcycle = generate(s, "MueLu", getMueLuOptions());
  const auto rv = solve(s, *(vcycle.linear_solver), "MueLu, one V-cycle");
  check(rv.converged && (rv.iterations == 0) && (rv.residual < 0.5),
        "MueLu: one V-cycle");
  // the initial guess is taken into account in the iterative mode: without
  // it, each V-cycle would give the residual of the first one
  check(iterate(s, *(vcycle.linear_solver), 30, "MueLu, stationary iteration") <
            1e-6,
        "MueLu: stationary iteration");
  // a preconditioner of mfem-mgis with a Krylov solver of Belos
  auto amg =
      generate(s, "BelosSolver", getBelosParameters("CG", "HypreBoomerAMG"));
  const auto ramg = solve(s, *(amg.linear_solver), "Belos CG, HypreBoomerAMG");
  check(isSolved(ramg) && isFaster(ramg, rj),
        "Belos CG, HypreBoomerAMG: resolution");
  // default solver (GMRES) without preconditioner
  auto gmres = generate(s, "BelosSolver",
                        Parameters{{"RelativeTolerance", 1e-10},
                                   {"MaximumNumberOfIterations", 2000}});
  check(isSolved(solve(s, *(gmres.linear_solver), "Belos GMRES")),
        "Belos GMRES: resolution");
  // a second operator, the first one times two
  auto a2 = mfem::HypreParMatrix(s.A);
  a2 *= 2.;
  check(isSolved(solve(s, a2, *(muelu.linear_solver), "Belos CG, MueLu, 2 A")),
        "Belos CG, MueLu: resolution with a second operator");
  check(isSolved(solve(s, *(muelu.linear_solver), "Belos CG, MueLu, A")),
        "Belos CG, MueLu: resolution with the first operator");
  // the prolongators and the restrictions of the first operator kept for
  // the following ones: the coarse operators of 2 A are those of A times
  // two, hence the same iterations (up to the rounding errors of the GPU)
  auto reuse = generate(
      s, "BelosSolver",
      getBelosParameters(
          "CG", "MueLu",
          Parameters{{"Parameters", Parameters{{"reuse: type", "RP"}}}}));
  const auto rr = solve(s, *(reuse.linear_solver), "Belos CG, MueLu, RP, A");
  const auto rr2 =
      solve(s, a2, *(reuse.linear_solver), "Belos CG, MueLu, RP, 2 A");
  const auto rr3 = solve(s, *(reuse.linear_solver), "Belos CG, MueLu, RP, A");
  check(isSolved(rr) && isSolved(rr2) && isSolved(rr3) &&
            (std::abs(rr2.iterations - rr.iterations) <= 1) &&
            (std::abs(rr3.iterations - rr.iterations) <= 1) && isFaster(rr, rj),
        "Belos CG, MueLu: resolutions with the hierarchy reused");
  // the whole hierarchy of A kept for 2 A (its smoothers and its coarse
  // operators): the preconditioner is no longer the multigrid of 2 A, which
  // changes the iterations, showing that the hierarchy is reused
  auto full = generate(
      s, "BelosSolver",
      getBelosParameters(
          "CG", "MueLu",
          Parameters{{"Parameters", Parameters{{"reuse: type", "full"}}}}));
  const auto rfull =
      solve(s, *(full.linear_solver), "Belos CG, MueLu, full, A");
  const auto rfull2 =
      solve(s, a2, *(full.linear_solver), "Belos CG, MueLu, full, 2 A");
  check(isSolved(rfull) && (std::abs(rfull.iterations - rr.iterations) <= 1) &&
            (rfull2.iterations > rr2.iterations + 3),
        "Belos CG, MueLu: hierarchy not reused");
}  // end of checkKrylovSolvers

/*!
 * \brief checks the algebraic multigrid of MueLu on the linear elasticity
 * \param[in] e1: unknowns ordered by nodes
 * \param[in] e2: unknowns ordered by vector dimension
 */
static void checkElasticity(ElasticitySystem<true>& e1,
                            ElasticitySystem<true>& e2) {
  using namespace mfem_mgis;
  auto m1 = generate(e1, "BelosSolver",
                     getBelosParameters("CG", "MueLu", getMueLuOptions()));
  const auto r1 = solve(e1, *(m1.linear_solver), "Belos CG, MueLu, byNODES");
  auto m2 = generate(e2, "BelosSolver",
                     getBelosParameters("CG", "MueLu", getMueLuOptions()));
  const auto r2 = solve(e2, *(m2.linear_solver), "Belos CG, MueLu, byVDIM");
  // several iterations: the multigrid is not a direct solver on this system
  check(isSolved(r1) && isSolved(r2) && (r1.iterations > 1),
        "Belos CG, MueLu: resolution of the linear elasticity");
  // the matrices assembled by MFEM for both orderings differ by round-off
  // errors (order of the sums), which may change a few aggregates: the
  // numbers of iterations are only close
  check(
      std::abs(r1.iterations - r2.iterations) <= std::max(3, r2.iterations / 5),
      "Belos CG, MueLu: the ordering of the unknowns changes the number "
      "of iterations");
  // aggregation of the unknowns one by one, without the rigid body modes
  auto m3 =
      generate(e2, "BelosSolver",
               getBelosParameters(
                   "CG", "MueLu",
                   getMueLuOptions(Parameters{{"number of equations", 1}})));
  const auto r3 =
      solve(e2, *(m3.linear_solver), "Belos CG, MueLu, scalar aggregation");
  check(isSolved(r3), "Belos CG, MueLu, scalar aggregation: resolution");
  check(isFaster(r1, r3) && isFaster(r2, r3),
        "Belos CG, MueLu: the rigid body modes do not reduce the number of "
        "iterations");
  // the hierarchy, built with the rigid body modes, reused for a second
  // operator
  auto m4 = generate(
      e1, "BelosSolver",
      getBelosParameters("CG", "MueLu",
                         getMueLuOptions(Parameters{{"reuse: type", "RP"}})));
  const auto r4 = solve(e1, *(m4.linear_solver), "Belos CG, MueLu, RP, A");
  auto a2 = mfem::HypreParMatrix(e1.A);
  a2 *= 2.;
  const auto r5 =
      solve(e1, a2, *(m4.linear_solver), "Belos CG, MueLu, RP, 2 A");
  check(isSolved(r4) && isSolved(r5) && isFaster(r4, r3) &&
            (std::abs(r5.iterations - r4.iterations) <= 1),
        "Belos CG, MueLu: hierarchy with the rigid body modes reused");
  auto jacobi =
      generate(e1, "BelosSolver", getBelosParameters("GMRES", "Ifpack2"));
  const auto rj = solve(e1, *(jacobi.linear_solver),
                        "Belos GMRES, Ifpack2 Jacobi, byNODES");
  check(isSolved(rj),
        "Belos GMRES, Ifpack2 Jacobi: resolution of the linear elasticity");
  // two level Schwarz preconditioner, with the null space of the linear
  // elasticity (9 iterations against 103 for the Jacobi relaxation on two
  // processes)
  auto f1 = generate(e1, "BelosSolver", getBelosParameters("CG", "FROSch"));
  const auto rf1 = solve(e1, *(f1.linear_solver), "Belos CG, FROSch, byNODES");
  check(isSolved(rf1) && isFaster(rf1, rj),
        "Belos CG, FROSch: resolution of the linear elasticity (byNODES)");
  auto f2 = generate(e2, "BelosSolver", getBelosParameters("CG", "FROSch"));
  const auto rf2 = solve(e2, *(f2.linear_solver), "Belos CG, FROSch, byVDIM");
  check(isSolved(rf2) && isFaster(rf2, rj),
        "Belos CG, FROSch: resolution of the linear elasticity (byVDIM)");
}  // end of checkElasticity

/*!
 * \brief checks the fast incomplete factorizations of Ifpack2 (ShyLU
 * FastILU), if Trilinos is built with them, refused otherwise
 * \param[in] s: system
 */
static void checkFastILU(LaplaceSystem<true>& s) {
  using namespace mfem_mgis;
  const auto filu =
      getBelosParameters("GMRES", "Ifpack2", Parameters{{"Type", "FAST_ILU"}});
  const auto e = getGenerationError(s, "BelosSolver", filu);
  if (!e.empty()) {
    check(e.find("ShyLU_NodeFastILU") != std::string::npos,
          "Ifpack2: FAST_ILU refused without explanation (" + e + ")");
    if (rank == 0) {
      std::cout << "Ifpack2, FAST_ILU: skipped (not available)\n";
    }
    return;
  }
  // the incomplete factorization without fill-in of Ifpack2, applied to
  // the diagonal block of each process, as the fast factorizations
  auto riluk = generate(
      s, "BelosSolver",
      getBelosParameters("GMRES", "Ifpack2", Parameters{{"Type", "RILUK"}}));
  const auto rr =
      solve(s, *(riluk.linear_solver), "Belos GMRES, Ifpack2 RILUK");
  // the fast factorizations approximate it: close numbers of iterations
  const auto close = [&rr](const Resolution& r) {
    return isSolved(r) && (std::abs(r.iterations - rr.iterations) <=
                           std::max(3, rr.iterations / 5));
  };
  auto ilu = generate(s, "BelosSolver", filu);
  const auto r1 =
      solve(s, *(ilu.linear_solver), "Belos GMRES, Ifpack2 FAST_ILU");
  check(isSolved(rr) && close(r1), "Belos GMRES, Ifpack2 FAST_ILU: resolution");
  // parameters given to the factorization: one sweep and one iteration of
  // the triangular solves, a coarser approximation, and an invalid number
  // of sweeps, refused when the preconditioner is created
  auto ilu1 = generate(
      s, "BelosSolver",
      getBelosParameters(
          "GMRES", "Ifpack2",
          Parameters{
              {"Type", "FAST_ILU"},
              {"Parameters", Parameters{{"sweeps", 1},
                                        {"triangular solve iterations", 1}}}}));
  const auto r3 = solve(s, *(ilu1.linear_solver),
                        "Belos GMRES, Ifpack2 FAST_ILU, one sweep");
  check(isSolved(r3) && (r3.iterations != r1.iterations),
        "Belos GMRES, Ifpack2 FAST_ILU: parameters of the factorization");
  const auto e2 = getGenerationError(
      s, "BelosSolver",
      getBelosParameters(
          "GMRES", "Ifpack2",
          Parameters{{"Type", "FAST_ILU"},
                     {"Parameters", Parameters{{"sweeps", 0}}}}));
  check(e2.find("sweeps") != std::string::npos,
        "Ifpack2: invalid parameter of FAST_ILU not refused (" + e2 + ")");
  // symmetric factorization, for the conjugate gradient
  auto ic = generate(
      s, "BelosSolver",
      getBelosParameters("CG", "Ifpack2", Parameters{{"Type", "FAST_IC"}}));
  const auto r2 = solve(s, *(ic.linear_solver), "Belos CG, Ifpack2 FAST_IC");
  check(close(r2), "Belos CG, Ifpack2 FAST_IC: resolution");
}  // end of checkFastILU

/*!
 * \brief checks the block preconditioners of Teko, if mfem-mgis uses Teko,
 * refused otherwise
 * \param[in] s: Laplace equation
 * \param[in] e1: linear elasticity, unknowns ordered by nodes
 * \param[in] e2: linear elasticity, unknowns ordered by vector dimension
 */
static void checkTeko(LaplaceSystem<true>& s,
                      ElasticitySystem<true>& e1,
                      ElasticitySystem<true>& e2) {
  using namespace mfem_mgis;
#ifdef MFEM_MGIS_HAVE_TEKO
  // block Jacobi, one block per component of the displacement, each block
  // approximated by a multigrid of MueLu with several levels
  const auto amg = Parameters{
      {"Type", "MueLu"}, {"verbosity", "none"}, {"coarse: max size", 20}};
  const auto jacobi = Parameters{
      {"Parameters", Parameters{{"Jacobi", Parameters{{"Type", "Block Jacobi"},
                                                      {"Inverse Type", "AMG"}}},
                                {"AMG", amg}}},
      {"Inverse", "Jacobi"}};
  auto j1 =
      generate(e1, "BelosSolver", getBelosParameters("CG", "Teko", jacobi));
  const auto r1 =
      solve(e1, *(j1.linear_solver), "Belos CG, Teko block Jacobi, byNODES");
  // the only block preconditioner of the definitions, without 'Inverse'
  const auto jacobi2 = remove(jacobi, std::vector<std::string>{"Inverse"});
  auto j2 =
      generate(e2, "BelosSolver", getBelosParameters("CG", "Teko", jacobi2));
  const auto r2 =
      solve(e2, *(j2.linear_solver), "Belos CG, Teko block Jacobi, byVDIM");
  check(isSolved(r1) && isSolved(r2) && (r1.iterations > 1),
        "Belos CG, Teko block Jacobi: resolution of the linear elasticity");
  check(
      std::abs(r1.iterations - r2.iterations) <= std::max(3, r2.iterations / 5),
      "Belos CG, Teko block Jacobi: the ordering of the unknowns changes the "
      "number of iterations");
  // default parameters, on the linear elasticity and on the Laplace
  // equation (one block)
  auto d1 = generate(e1, "BelosSolver", getBelosParameters("CG", "Teko"));
  check(isSolved(solve(e1, *(d1.linear_solver), "Belos CG, Teko, default")),
        "Belos CG, Teko: resolution with the default parameters");
  auto d2 = generate(s, "BelosSolver", getBelosParameters("CG", "Teko"));
  check(isSolved(solve(s, *(d2.linear_solver), "Belos CG, Teko, Laplace")),
        "Belos CG, Teko: resolution of the Laplace equation");
  // block Gauss-Seidel, a block of two components and a block of one
  // component, factorized by Amesos2 (not symmetric: GMRES)
  const auto gs = Parameters{
      {"Parameters",
       Parameters{{"Gauss-Seidel", Parameters{{"Type", "Block Gauss-Seidel"},
                                              {"Inverse Type", "Amesos2"}}}}},
      {"Blocks", list{2, 1}}};
  auto g = generate(e1, "BelosSolver", getBelosParameters("GMRES", "Teko", gs));
  const auto rg =
      solve(e1, *(g.linear_solver), "Belos GMRES, Teko block Gauss-Seidel");
  check(isSolved(rg),
        "Belos GMRES, Teko block Gauss-Seidel: resolution of the linear "
        "elasticity");
  // invalid parameters
  const auto refused = [&e1](const Parameters& options,
                             const std::string& msg) {
    const auto e = getGenerationError(
        e1, "BelosSolver", getBelosParameters("CG", "Teko", options));
    check(e.find(msg) != std::string::npos,
          "Teko: invalid parameters not refused (expected '" + msg +
              "', got '" + e + "')");
  };
  refused(Parameters{{"Inverse", "Unknown"}}, "unknown inverse 'Unknown'");
  refused(Parameters{{"Blocks", list{1, 1}}},
          "is not the number of components of the unknowns");
  refused(Parameters{{"Blocks", list{3, 0}}},
          "the blocks must have at least one component");
  refused(Parameters{{"Parameters", Parameters{{"Jacobi", 1}}}},
          "is not the definition of an inverse");
  refused(Parameters{{"Parameters",
                      Parameters{{"Jacobi",
                                  Parameters{{"Inverse Type", "Amesos2"}}}}}},
          "has no type");
  refused(Parameters{{"Parameters",
                      Parameters{{"Jacobi", Parameters{{"Type", "Unknown"}}}}}},
          "unknown type 'Unknown' of the inverse 'Jacobi'");
  refused(Parameters{{"Parameters",
                      Parameters{{"Jacobi",
                                  Parameters{{"Type", "Block Jacobi"},
                                             {"Inverse Type", "Unknown"}}}}}},
          "unknown inverse 'Unknown' given by the parameter 'Inverse Type'");
  refused(Parameters{{"Parameters", Parameters{{"AMG", amg}}}},
          "defines no block preconditioner of Teko");
  refused(Parameters{{"Parameters",
                      Parameters{{"A", Parameters{{"Type", "Block Jacobi"}}},
                                 {"B", Parameters{{"Type", "Block Jacobi"}}}}}},
          "the parameter 'Inverse' is required");
  refused(Parameters{{"Inverse", "Amesos2"}},
          "is not a block preconditioner of Teko");
  refused(Parameters{{"Parameters",
                      Parameters{{"Jacobi",
                                  Parameters{{"Type", "Block Jacobi"},
                                             {"Inverse Type 0", "Amesos2"}}}}}},
          "invalid parameter 'Inverse Type 0'");
  refused(Parameters{{"Parameters",
                      Parameters{
                          {"Ifpack2", Parameters{{"Type", "Block Jacobi"}}}}}},
          "has the name of a solver or a preconditioner of Stratimikos");
  refused(
      Parameters{
          {"Parameters",
           Parameters{{"Jacobi", Parameters{{"Type", "Block Jacobi"},
                                            {"Inverse Type", "CG"}}},
                      {"CG", Parameters{{"Type", "Belos"},
                                        {"Use Preconditioner", "Amesos2"}}}}}},
      "is not a preconditioner");
#else  /* MFEM_MGIS_HAVE_TEKO */
  const auto e =
      getGenerationError(e1, "BelosSolver", getBelosParameters("CG", "Teko"));
  check(e.find("enable-teko") != std::string::npos,
        "Teko: refused without explanation (" + e + ")");
  static_cast<void>(s);
  static_cast<void>(e2);
  if (rank == 0) {
    std::cout << "Teko: skipped (not available)\n";
  }
#endif /* MFEM_MGIS_HAVE_TEKO */
}  // end of checkTeko

/*!
 * \brief checks the options of mfem-mgis and of Belos
 * \param[in] s: system
 * \param[in] elasticity: linear elasticity, a system of another size
 */
static void checkOptions(LaplaceSystem<true>& s,
                         ElasticitySystem<true>& elasticity) {
  using namespace mfem_mgis;
  const auto jacobi = Parameters{{"Name", "Ifpack2"}};
  auto ls = generate(s, "BelosSolver",
                     Parameters{{"Solver", "CG"},
                                {"RelativeTolerance", 1e-6},
                                {"AbsoluteTolerance", 0.},
                                {"Preconditioner", jacobi}});
  const auto rtol = getRelativeTolerance(*(ls.linear_solver));
  check(rtol.has_value() && (std::abs(*rtol - 1e-6) < 1e-12),
        "BelosSolver: relative tolerance");
  const auto r = solve(s, *(ls.linear_solver), "Belos CG, tolerance 1e-6");
  check(r.converged && (r.residual < 1e-5) && (r.residual > 1e-10),
        "Belos CG, tolerance 1e-6: resolution");
  check(setRelativeTolerance(*(ls.linear_solver), 1e-10) &&
            (std::abs(getRelativeTolerance(*(ls.linear_solver)).value_or(0) -
                      1e-10) < 1e-14),
        "BelosSolver: change of the relative tolerance");
  const auto r2 = solve(s, *(ls.linear_solver), "Belos CG, tolerance 1e-10");
  check(isSolved(r2) && (r2.iterations > r.iterations),
        "Belos CG, tolerance 1e-10: resolution");
  // the parameters of Belos override the options of mfem-mgis
  auto ls2 = generate(
      s, "BelosSolver",
      Parameters{{"Solver", "CG"},
                 {"RelativeTolerance", 1e-10},
                 {"Parameters", Parameters{{"Convergence Tolerance", 1e-4}}},
                 {"Preconditioner", jacobi}});
  check(std::abs(getRelativeTolerance(*(ls2.linear_solver)).value_or(0) -
                 1e-4) < 1e-12,
        "BelosSolver: relative tolerance given by the parameters of Belos");
  // maximum number of iterations
  auto ls3 = generate(s, "BelosSolver",
                      Parameters{{"Solver", "CG"},
                                 {"RelativeTolerance", 1e-10},
                                 {"MaximumNumberOfIterations", 3},
                                 {"Preconditioner", jacobi}});
  const auto r3 = solve(s, *(ls3.linear_solver), "Belos CG, 3 iterations");
  check((!r3.converged) && (r3.iterations == 3),
        "Belos CG: resolution stopped by the maximum number of iterations");
  // the hybrid GMRES (polynomial preconditioning) gives the convergence
  // parameters to its outer solver
  auto hybrid = generate(s, "BelosSolver",
                         Parameters{{"Solver", "HYBRID BLOCK GMRES"},
                                    {"RelativeTolerance", 1e-6},
                                    {"MaximumNumberOfIterations", 2000}});
  check(std::abs(getRelativeTolerance(*(hybrid.linear_solver)).value_or(0) -
                 1e-6) < 1e-12,
        "Belos hybrid GMRES: relative tolerance");
  const auto rh = solve(s, *(hybrid.linear_solver), "Belos hybrid GMRES, 1e-6");
  check(rh.converged && (rh.residual < 1e-5) && (rh.residual > 1e-10),
        "Belos hybrid GMRES, tolerance 1e-6: resolution");
  check(setRelativeTolerance(*(hybrid.linear_solver), 1e-10),
        "Belos hybrid GMRES: change of the relative tolerance");
  // a second operator of another size (the linear elasticity, the unknowns
  // being seen as scalar unknowns): the polynomial is built again
  check(isSolved(solve(elasticity, elasticity.A, *(hybrid.linear_solver),
                       "Belos hybrid GMRES, 1e-10, elasticity")),
        "Belos hybrid GMRES, tolerance 1e-10: resolution of a second "
        "system");
  auto a2 = mfem::HypreParMatrix(s.A);
  a2 *= 2.;
  // the aliases of the hybrid GMRES, and its default tolerance (the one of
  // its outer solver)
  auto poly = generate(
      s, "BelosSolver",
      Parameters{{"Solver", "GMRESPOLY"}, {"MaximumNumberOfIterations", 2000}});
  check(getRelativeTolerance(*(poly.linear_solver)).has_value(),
        "Belos hybrid GMRES: default relative tolerance");
  check(isSolved(solve(s, *(poly.linear_solver), "Belos GMRESPOLY"), 1e-6),
        "Belos GMRESPOLY: resolution");
  // GCRODR (recycling GMRES) and a native GMRES of Tpetra, with a second
  // operator; the native solvers of Tpetra ignore the unknown parameters
  // (the output frequency set by mfem-mgis is not given to them)
  for (const auto* const solver : {"GCRODR", "TPETRA GMRES"}) {
    const auto p = Parameters{{"Solver", solver},
                              {"RelativeTolerance", 1e-10},
                              {"MaximumNumberOfIterations", 2000}};
    auto verbose = p;
    verbose.insert(may_throw, "VerbosityLevel", 2);
    const auto error = getGenerationError(s, "BelosSolver", verbose);
    check(error.empty(), "Belos " + std::string{solver} +
                             ": verbosity refused ('" + error + "')");
    auto krylov = generate(s, "BelosSolver", p);
    check(isSolved(solve(s, *(krylov.linear_solver), solver)) &&
              isSolved(solve(s, a2, *(krylov.linear_solver), solver)),
          "Belos " + std::string{solver} + ": resolutions");
  }
  // parameters of Ifpack2 absent from `Ifpack2::getValidParameters`
  const auto accepted = [&s](const std::string& type, const Parameters& p) {
    const auto options = Parameters{{"Type", type}, {"Parameters", p}};
    const auto error = getGenerationError(
        s, "BelosSolver", getBelosParameters("CG", "Ifpack2", options));
    check(error.empty(),
          "Ifpack2, " + type + ": valid parameters refused ('" + error + "')");
  };
  accepted("RELAXATION",
           Parameters{{"relaxation: fix tiny diagonal entries", true},
                      {"relaxation: inner sweeps", 2}});
  accepted("CHEBYSHEV", Parameters{{"chebyshev: boost factor", 1.1}});
  accepted("RILUK", Parameters{{"fact: iluk overalloc", 2.}});
  // blocks of three rows (the block relaxations read the matrix when their
  // parameters are set)
  auto blocks = generate(
      s, "BelosSolver",
      getBelosParameters(
          "GMRES", "Ifpack2",
          Parameters{
              {"Type", "BLOCK RELAXATION"},
              {"Parameters", Parameters{{"partitioner: local parts", -3}}}}));
  check(isSolved(solve(s, *(blocks.linear_solver),
                       "Belos GMRES, Ifpack2 block relaxation")),
        "Belos GMRES, Ifpack2 block relaxation: resolution");
  auto schwarz = generate(
      s, "BelosSolver",
      getBelosParameters(
          "GMRES", "Ifpack2",
          Parameters{
              {"Type", "SCHWARZ"},
              {"Parameters",
               Parameters{{"schwarz: inner preconditioner name", "RILUK"}}}}));
  check(isSolved(solve(s, *(schwarz.linear_solver),
                       "Belos GMRES, Ifpack2 Schwarz, RILUK")),
        "Belos GMRES, Ifpack2 Schwarz: resolution");
  // all the preconditioners of Ifpack2 are created, their parameters being
  // checked by a preconditioner built without matrix
  const auto error = getGenerationError(
      s, "BelosSolver",
      getBelosParameters("CG", "Ifpack2", Parameters{{"Type", "NoSuchType"}}));
  auto count = 0;
  const auto b = error.find("(supported: ");
  const auto e = error.find(')', b);
  if ((b != std::string::npos) && (e != std::string::npos)) {
    // names in quotes, separated by commas
    auto names = std::istringstream{error.substr(b + 12, e - b - 12)};
    for (auto n = std::string{}; std::getline(names, n, '\'');) {
      if (!std::getline(names, n, '\'')) {
        break;
      }
      ++count;
      if (n == "HIPTMAIR") {
        // needs the matrices of its auxiliary space, which mfem-mgis does
        // not give
        continue;
      }
      // the factorization of Amesos2 reads the sublist `Amesos2` only
      accepted(n, n == "AMESOS2" ? Parameters{{"Amesos2", Parameters{}}}
                                 : Parameters{});
    }
  }
  check(count > 10, "Ifpack2: list of the preconditioners");
  // the parameters of ML, translated by MueLu, are not checked at the
  // creation (MueLu 17.2.1 translates them to an option of the aggregation
  // which its version on Kokkos does not support: no resolution)
  const auto ml = getGenerationError(
      s, "BelosSolver",
      getBelosParameters(
          "CG", "MueLu",
          Parameters{{"Parameters", Parameters{{"parameterlist: syntax", "ml"},
                                               {"ML output", 0}}}}));
  check(ml.empty(), "MueLu: parameters of ML refused ('" + ml + "')");
}  // end of checkOptions

/*!
 * \brief checks the resolution of a system with a process owning no
 * unknown
 * \param[in] s: system
 */
static void checkEmptyProcess(ElasticitySystem<true>& s) {
  using namespace mfem_mgis;
  auto direct = generate(s, "Amesos2Solver", {});
  check(isSolved(solve(s, *(direct.linear_solver),
                       "Amesos2Solver, process without unknowns"),
                 1e-10),
        "Amesos2Solver: resolution with a process owning no unknown");
  auto muelu = generate(s, "BelosSolver",
                        getBelosParameters("CG", "MueLu", getMueLuOptions()));
  check(isSolved(solve(s, *(muelu.linear_solver),
                       "Belos CG, MueLu, process without unknowns")),
        "Belos CG, MueLu: resolution with a process owning no unknown");
  auto ifpack2 =
      generate(s, "BelosSolver", getBelosParameters("CG", "Ifpack2"));
  check(isSolved(solve(s, *(ifpack2.linear_solver),
                       "Belos CG, Ifpack2, process without unknowns")),
        "Belos CG, Ifpack2: resolution with a process owning no unknown");
  // FROSch and Teko do not support empty subdomains: refused by all the
  // processes (in sequential, the only process owns all the unknowns)
  auto size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#ifdef MFEM_MGIS_HAVE_TEKO
  auto teko = generate(s, "BelosSolver", getBelosParameters("CG", "Teko"));
  if (size == 1) {
    check(isSolved(solve(s, *(teko.linear_solver), "Belos CG, Teko")),
          "Belos CG, Teko: resolution of the small system");
  } else {
    auto refused = false;
    try {
      static_cast<void>(solve(s, *(teko.linear_solver),
                              "Belos CG, Teko, process without unknowns"));
    } catch (std::exception& e) {
      refused = std::string{e.what()}.find("each process must own unknowns") !=
                std::string::npos;
    }
    check(refused, "Belos CG, Teko: process owning no unknown not refused");
  }
#endif /* MFEM_MGIS_HAVE_TEKO */
  auto frosch = generate(s, "BelosSolver", getBelosParameters("CG", "FROSch"));
  if (size == 1) {
    check(isSolved(solve(s, *(frosch.linear_solver), "Belos CG, FROSch")),
          "Belos CG, FROSch: resolution of the small system");
  } else {
    auto refused = false;
    try {
      static_cast<void>(solve(s, *(frosch.linear_solver),
                              "Belos CG, FROSch, process without unknowns"));
    } catch (std::exception& e) {
      refused = std::string{e.what()}.find("each process must own unknowns") !=
                std::string::npos;
    }
    check(refused, "Belos CG, FROSch: process owning no unknown not refused");
  }
}  // end of checkEmptyProcess

/*!
 * \brief checks that invalid parameters are refused
 * \param[in] s: system
 */
static void checkRefusals(LaplaceSystem<true>& s) {
  using namespace mfem_mgis;
  const auto refused = [&s](const std::string& name,
                            const Parameters& parameters,
                            const std::string& msg) {
    const auto error = getGenerationError(s, name, parameters);
    check(error.find(msg) != std::string::npos,
          name + ": invalid parameters not refused ('" + error + "')");
  };
  refused("BelosSolver", Parameters{{"Solver", "NoSuchSolver"}},
          "NoSuchSolver");
  refused("BelosSolver", Parameters{{"AbsoluteTolerance", 1e-3}},
          "only have a relative tolerance");
  refused(
      "BelosSolver",
      Parameters{{"Parameters", Parameters{{"Orthogonalization",
                                            std::vector<Parameter>{1, 2}}}}},
      "invalid value of the parameter 'Orthogonalization'");
  refused("Amesos2Solver", Parameters{{"Solver", "NoSuchSolver"}},
          "not available");
  refused("BelosSolver",
          getBelosParameters(
              "CG", "MueLu",
              Parameters{{"Parameters", Parameters{{"coarse: max sise", 20}}}}),
          "unknown parameter 'coarse: max sise'");
  refused("BelosSolver",
          getBelosParameters(
              "CG", "Ifpack2",
              Parameters{{"Parameters", Parameters{{"relaxation: sweep", 2}}}}),
          "relaxation: sweep");
  refused(
      "BelosSolver",
      Parameters{{"Solver", "CG"},
                 {"Parameters", Parameters{{"Convergence Tolerence", 1e-4}}}},
      "Convergence Tolerence");
  refused(
      "BelosSolver",
      Parameters{{"Solver", "HYBRID BLOCK GMRES"},
                 {"Parameters", Parameters{{"Outer Solver Params",
                                            Parameters{{"Num Blocs", 10}}}}}},
      "Num Blocs");
  // the outer solver of the hybrid GMRES is built by the generic factory
  // of Belos, without the native solvers of Tpetra
  refused(
      "BelosSolver",
      Parameters{{"Solver", "HYBRID BLOCK GMRES"},
                 {"Parameters", Parameters{{"Outer Solver", "TPETRA GMRES"}}}},
      "unknown outer solver");
  refused("BelosSolver",
          Parameters{{"Solver", "LSQR"}, {"RelativeTolerance", 1e-8}},
          "no relative tolerance");
  // LSQR and the native solvers of Tpetra read parameters absent from the
  // list of their valid parameters
  using UncheckedParameter = std::tuple<std::string, std::string, Parameter>;
  const auto unchecked = std::vector<UncheckedParameter>{
      {"LSQR", "lambda", 0.}, {"TPETRA GMRES", "Orthogonalization", "ICGS"}};
  for (const auto& [solver, parameter, value] : unchecked) {
    const auto error = getGenerationError(
        s, "BelosSolver",
        Parameters{{"Solver", solver},
                   {"Parameters", Parameters{{parameter, value}}}});
    check(error.empty(), "Belos " + solver + ": parameter '" + parameter +
                             "' refused ('" + error + "')");
  }
  refused("BelosSolver",
          getBelosParameters(
              "CG", "FROSch",
              Parameters{{"Parameters",
                          Parameters{{"DirichletBoundary Numerical Tolerance",
                                      1e-12}}}}),
          "DirichletBoundary Numerical Tolerance");
  refused(
      "BelosSolver",
      Parameters{{"Preconditioner",
                  Parameters{{"Name", "Ifpack2"},
                             {"Options", Parameters{{"Type", "NoSuchType"}}}}}},
      "Ifpack2");
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
  auto mesh =
      mfem::Mesh::MakeCartesian3D(16, 16, 16, mfem::Element::HEXAHEDRON);
  auto beam = mfem::Mesh::MakeCartesian3D(12, 3, 3, mfem::Element::HEXAHEDRON,
                                          4., 1., 1.);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    auto s = LaplaceSystem<true>(mesh, fec);
    auto e1 = ElasticitySystem<true>(beam, fec, mfem::Ordering::byNODES);
    auto e2 = ElasticitySystem<true>(beam, fec, mfem::Ordering::byVDIM);
    checkDirectSolvers(s, e1, e2);
    checkKrylovSolvers(s);
    checkElasticity(e1, e2);
    checkFastILU(s);
    checkTeko(s, e1, e2);
    checkOptions(s, e1);
    checkRefusals(s);
    // a smaller mesh, whose elements are on the first process (two
    // processes), or split between the first and the last processes (more
    // processes): the other processes own no unknown
    auto beam2 = mfem::Mesh::MakeCartesian3D(8, 2, 2, mfem::Element::HEXAHEDRON,
                                             4., 1., 1.);
    auto size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    auto partitioning =
        std::vector<int>(static_cast<std::size_t>(beam2.GetNE()), 0);
    if (size > 2) {
      for (auto i = partitioning.size() / 2; i != partitioning.size(); ++i) {
        partitioning[i] = size - 1;
      }
    }
    auto e3 = ElasticitySystem<true>(beam2, fec, mfem::Ordering::byNODES,
                                     partitioning.data());
    checkEmptyProcess(e3);
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (rank == 0) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
