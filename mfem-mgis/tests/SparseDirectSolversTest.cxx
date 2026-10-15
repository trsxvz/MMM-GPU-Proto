/*!
 * \file   SparseDirectSolversTest.cxx
 * \brief  Test of the direct solvers of MUMPS, SuperLU_DIST and STRUMPACK,
 * and of the preconditioner of STRUMPACK, built by the factories, on the
 * Laplace equation on a cube, in parallel (the generic test
 * `LinearSolverFactoriesTest` solves it with each entry and its default
 * options):
 *
 * - MUMPS accepts each reordering strategy and the tolerance of its block
 *   low-rank factorization, and solves the system (the strategy used is not
 *   reported by MFEM);
 * - SuperLU_DIST solves the system with its options, reuses its analysis
 *   for a new operator of the same pattern and is built again for an
 *   operator of another size; the ParMETIS ordering with the parallel
 *   symbolic factorization is refused on one process and with the reuse of
 *   the analysis, and accepted in both cases with the sequential symbolic
 *   factorization;
 * - STRUMPACK solves the system in double and in mixed precision, with the
 *   options of STRUMPACK given by the `Options` parameter, reuses its
 *   reordering for a new operator of the same pattern and is built again
 *   for an operator of another size;
 * - the Krylov solver of STRUMPACK, with the compression of the frontal
 *   matrices, reports its convergence, and its non convergence when it
 *   stops at the maximum number of iterations;
 * - the approximate factorization of STRUMPACK (block low-rank compression)
 *   preconditions the conjugate gradient of MFEM;
 * - invalid options are refused.
 * \date   15/10/2026
 */

#include <cmath>
#include <string>
#include <type_traits>
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
/*!
 * \brief tolerance on the relative residual of the direct solvers of MUMPS
 * (SuperLU_DIST and STRUMPACK are only available in double precision)
 */
static constexpr double direct_tolerance =
    std::is_same_v<mfem::real_t, double> ? 1e-12 : 1e-5;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

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

/*!
 * \return the relative residual of the resolution of the system of a
 * matrix, with the right hand side of a system
 * \param[in] s: system
 * \param[in] a: matrix
 * \param[in] ls: linear solver
 * \param[in] name: name of the resolution
 */
static double solve(LaplaceSystem<true>& s,
                    const mfem::HypreParMatrix& a,
                    mfem::Solver& ls,
                    const std::string& name) {
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.SetOperator(a);
  ls.Mult(s.B, x);
  const auto r = s.getRelativeResidual(a, x);
  if (rank == 0) {
    std::cout << name << ": residual " << r << '\n';
  }
  return r;
}  // end of solve

/*!
 * \return the relative residual of the resolution of a system
 * \param[in] s: system
 * \param[in] ls: linear solver
 * \param[in] name: name of the resolution
 */
static double solve(LaplaceSystem<true>& s,
                    mfem::Solver& ls,
                    const std::string& name) {
  return solve(s, s.A, ls, name);
}  // end of solve

static void run(LaplaceSystem<true>& s,
                [[maybe_unused]] LaplaceSystem<true>& s2) {
  using namespace mfem_mgis;
  // a second operator of the same pattern, as in the iterations of a Newton
  // solver: the matrix of the system times two
  auto a2 = mfem::HypreParMatrix(s.A);
  a2 *= 2.;
#ifdef MFEM_USE_MUMPS
  // MUMPS and its reordering strategies (the ones which are not available
  // are replaced by MUMPS with a warning)
  for (const auto& r : {"Automatic", "AMD", "PORD", "METIS"}) {
    auto mumps =
        generate(s, "MUMPSSolver",
                 Parameters{{"ReorderingStrategy", r}, {"VerbosityLevel", 0}});
    check(solve(s, *(mumps.linear_solver), std::string{"MUMPSSolver, "} + r) <
              direct_tolerance,
          std::string{"MUMPSSolver, "} + r + ": resolution");
  }
  // the block low-rank factorization is approximate; on this small mesh,
  // MUMPS activates it (ICNTL(35)) but compresses no front, which are too
  // small: the resolution is exact. The option is only checked to be
  // accepted, the use of the block low-rank factorization not being
  // reported by MFEM
  {
    auto blr = generate(s, "MUMPSSolver", Parameters{{"BLRTolerance", 1e-6}});
    check(solve(s, *(blr.linear_solver), "MUMPSSolver, BLR") < direct_tolerance,
          "MUMPSSolver, BLR: resolution");
  }
  check(getGenerationError(s, "MUMPSSolver",
                           Parameters{{"ReorderingStrategy", "RCM"}})
                .find("invalid reordering strategy") != std::string::npos,
        "MUMPSSolver: invalid reordering strategy not refused");
  check(getGenerationError(s, "MUMPSSolver", Parameters{{"BLRTolerance", -1.}})
                .find("invalid tolerance") != std::string::npos,
        "MUMPSSolver: negative tolerance not refused");
#endif /* MFEM_USE_MUMPS */
#ifdef MFEM_USE_SUPERLU
  // SuperLU_DIST, its options, the reuse of its analysis
  auto superlu = generate(s, "SuperLUSolver",
                          Parameters{{"ColumnPermutation", "METIS_AT_PLUS_A"},
                                     {"RowPermutation", "LargeDiag_MC64"},
                                     {"IterativeRefinement", "Double"},
                                     {"Equilibrate", true},
                                     {"ReuseAnalysis", true}});
  check(solve(s, *(superlu.linear_solver), "SuperLUSolver") < 1e-12,
        "SuperLUSolver: resolution");
  check(solve(s, a2, *(superlu.linear_solver),
              "SuperLUSolver, analysis reused, second operator") < 1e-12,
        "SuperLUSolver: resolution of a second operator with the analysis "
        "reused");
  check(solve(s2, *(superlu.linear_solver), "SuperLUSolver, smaller mesh") <
            1e-12,
        "SuperLUSolver: resolution on a smaller mesh");
  {
    auto size = int{};
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1) {
      check(getGenerationError(s, "SuperLUSolver",
                               Parameters{{"ColumnPermutation", "ParMETIS"}})
                    .find("at least two processes") != std::string::npos,
            "SuperLUSolver: ParMETIS on one process not refused");
    } else {
      auto parmetis = generate(s, "SuperLUSolver",
                               Parameters{{"ColumnPermutation", "ParMETIS"}});
      check(solve(s, *(parmetis.linear_solver), "SuperLUSolver, ParMETIS") <
                1e-12,
            "SuperLUSolver, ParMETIS: resolution");
      check(solve(s, *(parmetis.linear_solver),
                  "SuperLUSolver, ParMETIS, second factorization") < 1e-12,
            "SuperLUSolver, ParMETIS: second resolution");
    }
  }
  check(getGenerationError(s, "SuperLUSolver",
                           Parameters{{"ColumnPermutation", "ParMETIS"},
                                      {"ReuseAnalysis", true}})
                .find("can't be reused") != std::string::npos,
        "SuperLUSolver: reuse of the analysis with ParMETIS not refused");
  // the ParMETIS ordering with the sequential symbolic factorization, on
  // any number of processes, and with the reuse of the analysis
  {
    auto parmetis =
        generate(s, "SuperLUSolver",
                 Parameters{{"ColumnPermutation", "ParMETIS"},
                            {"ParallelSymbolicFactorization", false},
                            {"ReuseAnalysis", true}});
    const auto name = std::string{
        "SuperLUSolver, ParMETIS, sequential symbolic factorization"};
    check((solve(s, *(parmetis.linear_solver), name) < 1e-12) &&
              (solve(s, a2, *(parmetis.linear_solver),
                     name + ", second operator") < 1e-12),
          name + ": resolution");
  }
  check(getGenerationError(s, "SuperLUSolver",
                           Parameters{{"ParallelSymbolicFactorization", true},
                                      {"ReuseAnalysis", true}})
                .find("can't be reused") != std::string::npos,
        "SuperLUSolver: reuse of the analysis with the parallel symbolic "
        "factorization not refused");
  // operators set without resolution, and a solver destroyed before any
  // resolution
  {
    auto lazy = generate(s, "SuperLUSolver", Parameters{});
    lazy.linear_solver->SetOperator(s2.A);
    check(solve(s, *(lazy.linear_solver), "SuperLUSolver, operator set twice") <
              1e-12,
          "SuperLUSolver: resolution after two operators");
    auto unused = generate(s, "SuperLUSolver", Parameters{});
    unused.linear_solver->SetOperator(s.A);
  }
  check(getGenerationError(s, "SuperLUSolver",
                           Parameters{{"ColumnPermutation", "AMD"}})
                .find("invalid value") != std::string::npos,
        "SuperLUSolver: invalid column permutation not refused");
#endif /* MFEM_USE_SUPERLU */
#ifdef MFEM_USE_STRUMPACK
  // STRUMPACK, in double and mixed precisions
  for (const auto& n : {"STRUMPACKSolver", "STRUMPACKMixedPrecisionSolver"}) {
    auto strumpack = generate(
        s, n,
        Parameters{{"Options", Parameters{{"sp_reordering_method", "metis"}}}});
    check((solve(s, *(strumpack.linear_solver), n) < 1e-10) &&
              hasConverged(*(strumpack.linear_solver)),
          std::string{n} + ": resolution");
    // the reordering is reused for an operator of the same pattern, and
    // computed again for an operator of another size (the iterative
    // refinement of the mixed precision solver stops by default at the
    // relative tolerance of STRUMPACK, 1e-6)
    auto reused = generate(
        s, n,
        Parameters{{"ReuseAnalysis", true}, {"RelativeTolerance", 1e-12}});
    const auto name = std::string{n} + ", analysis reused";
    check((solve(s, *(reused.linear_solver), name) < 1e-10) &&
              (solve(s, a2, *(reused.linear_solver),
                     name + ", second operator") < 1e-10) &&
              (solve(s2, *(reused.linear_solver), name + ", smaller mesh") <
               1e-10),
          name + ": resolutions");
  }
  // compression with a Krylov solver of STRUMPACK
  auto blr = generate(
      s, "STRUMPACKSolver",
      Parameters{{"RelativeTolerance", 1e-10},
                 {"Options", Parameters{{"sp_compression", "blr"},
                                        {"sp_compression_min_sep_size", 16},
                                        {"blr_rel_tol", 1e-2},
                                        {"sp_Krylov_solver", "pgmres"}}}});
  check(
      (solve(s, *(blr.linear_solver), "STRUMPACKSolver, BLR, GMRES") < 1e-8) &&
          hasConverged(*(blr.linear_solver)) &&
          (getNumberOfIterationsAtConvergence(*(blr.linear_solver))
               .value_or(0) > 0),
      "STRUMPACKSolver, BLR: resolution");
  // a resolution stopped by the maximum number of iterations, which
  // STRUMPACK reports as a success
  auto stopped = generate(
      s, "STRUMPACKSolver",
      Parameters{{"RelativeTolerance", 1e-14},
                 {"MaximumNumberOfIterations", 2},
                 {"Options", Parameters{{"sp_compression", "blr"},
                                        {"sp_compression_min_sep_size", 16},
                                        {"blr_rel_tol", 1e-1},
                                        {"sp_Krylov_solver", "pgmres"}}}});
  check((solve(s, *(stopped.linear_solver),
               "STRUMPACKSolver, BLR, GMRES, 2 iterations") > 1e-8) &&
            (!hasConverged(*(stopped.linear_solver))),
        "STRUMPACKSolver: resolution stopped by the maximum number of "
        "iterations reported as converged");
  // the approximate factorization as a preconditioner of MFEM
  auto cg = generate(
      s, "CGSolver",
      Parameters{
          {"RelativeTolerance", 1e-10},
          {"MaximumNumberOfIterations", 200},
          {"Preconditioner",
           Parameters{
               {"Name", "STRUMPACKPreconditioner"},
               {"Options",
                Parameters{
                    {"Options", Parameters{{"sp_compression_min_sep_size", 16},
                                           {"blr_rel_tol", 1e-2}}}}}}}});
  check((solve(s, *(cg.linear_solver), "CGSolver, STRUMPACKPreconditioner") <
         1e-8) &&
            hasConverged(*(cg.linear_solver)),
        "CGSolver, STRUMPACKPreconditioner: resolution");
  check(getGenerationError(
            s, "STRUMPACKSolver",
            Parameters{{"Options", Parameters{{"--sp_compression", "blr"}}}})
                .find("leading dashes") != std::string::npos,
        "STRUMPACKSolver: option with leading dashes not refused");
#endif /* MFEM_USE_STRUMPACK */
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
  auto mesh2 = mfem::Mesh::MakeCartesian3D(8, 8, 8, mfem::Element::HEXAHEDRON);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    auto s = LaplaceSystem<true>(mesh, fec);
    auto s2 = LaplaceSystem<true>(mesh2, fec);
    run(s, s2);
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (rank == 0) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
