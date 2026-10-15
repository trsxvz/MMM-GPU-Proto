/*!
 * \file   MKLSolversTest.cxx
 * \brief  Test of the sparse direct solvers of oneMKL, PARDISO in
 * sequential and the Cluster Sparse Solver in parallel (`--parallel`), on
 * the Laplace equation on a cube:
 *
 * - the system is solved with the non symmetric, symmetric indefinite and
 *   symmetric positive definite types of matrix (the upper triangle being
 *   given to oneMKL for the symmetric types);
 * - the analysis is reused for a second operator of the same pattern (the
 *   matrix times two), and computed again for an operator of another size
 *   and, in sequential, for an operator of the same size and of another
 *   pattern (number of analyses checked);
 * - a singular matrix (a null row whose right hand side is not null, an
 *   inconsistent system) is reported, oneMKL perturbing the null pivot;
 * - invalid options are refused.
 * \date   15/10/2026
 */

#include <array>
#include <cmath>
#include <memory>
#include <algorithm>
#include <tuple>
#include <string>
#include <utility>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MPI.hxx"
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
//! \brief tolerance on the relative residual of the direct solvers
static constexpr double tolerance =
    std::is_same_v<mfem::real_t, double> ? 1e-12 : 1e-5;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

//! \return the name of the solver of oneMKL
static const char* getSolverName() {
  return parallel ? "MKLClusterSparseSolver" : "MKLPardisoSolver";
}  // end of getSolverName

/*!
 * \return the linear solver of the factory, throws on failure
 * \param[in] s: system
 * \param[in] parameters: parameters of the linear solver
 */
template <bool is_parallel>
static mfem_mgis::LinearSolverHandler generate(
    LaplaceSystem<is_parallel>& s, const mfem_mgis::Parameters& parameters) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto ls = LinearSolverFactory<is_parallel>::getFactory().generate(
      ctx, getSolverName(), s.fespace, parameters);
  if (isInvalid(ls)) {
    throw std::runtime_error(ctx.getErrorMessage());
  }
  return ls;
}  // end of generate

/*!
 * \return the relative residual of the resolution of the system of a
 * matrix, with the right hand side of a system, and if the solver reports
 * a success
 * \param[in] s: system
 * \param[in] a: matrix
 * \param[in] ls: linear solver
 * \param[in] name: name of the resolution
 */
template <bool is_parallel>
static bool solve(LaplaceSystem<is_parallel>& s,
                  const typename LaplaceSystem<is_parallel>::Matrix& a,
                  mfem::Solver& ls,
                  const std::string& name) {
  using namespace mfem_mgis;
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.SetOperator(a);
  ls.Mult(s.B, x);
  const auto r = s.getRelativeResidual(a, x);
  const auto converged = hasConverged(ls);
  if (rank == 0) {
    std::cout << name << ": residual " << r
              << (converged ? "" : " (failure reported)") << '\n';
  }
  return converged && (r < tolerance);
}  // end of solve

template <bool is_parallel>
static void run(LaplaceSystem<is_parallel>& s, LaplaceSystem<is_parallel>& s2) {
  using namespace mfem_mgis;
  using Matrix = typename LaplaceSystem<is_parallel>::Matrix;
  auto a2 = Matrix(s.A);
  a2 *= 2.;
  // types of matrix (the Laplace matrix, with the eliminated degrees of
  // freedom on the diagonal, is symmetric positive definite)
  const auto types = std::array<std::pair<const char*, Parameters>, 3>{
      {{"non symmetric", Parameters{}},
       {"symmetric indefinite", Parameters{{"Symmetric", true}}},
       {"symmetric positive definite",
        Parameters{{"Symmetric", true}, {"PositiveDefinite", true}}}}};
  for (const auto& [n, p] : types) {
    auto ls = generate(s, p);
    check(solve(s, s.A, *(ls.linear_solver),
                std::string{getSolverName()} + ", " + n),
          std::string{n} + ": resolution");
  }
  // reuse of the analysis, then an operator of another size
  {
    auto ls = generate(s, Parameters{{"ReuseAnalysis", true}});
    const auto name = std::string{getSolverName()} + ", analysis reused";
    check(solve(s, s.A, *(ls.linear_solver), name) &&
              solve(s, a2, *(ls.linear_solver), name + ", second operator"),
          "reuse of the analysis");
    check(getNumberOfAnalyses(*(ls.linear_solver)) == size_type{1},
          "the analysis is not reused for an operator of the same pattern");
    check(solve(s2, s2.A, *(ls.linear_solver), name + ", smaller mesh"),
          "resolution on a smaller mesh");
    check(getNumberOfAnalyses(*(ls.linear_solver)) == size_type{2},
          "the analysis is not computed again for an operator of another "
          "size");
    // an operator of the same size and of another pattern (two coupling
    // terms added between the first and the last unknowns)
    if constexpr (!is_parallel) {
      std::ignore = s.A.HostReadI();
      std::ignore = s.A.HostReadJ();
      std::ignore = s.A.HostReadData();
      const auto n = s.A.Height();
      auto extra = mfem::SparseMatrix(n, n);
      extra.Add(0, n - 1, 1e-6);
      extra.Add(n - 1, 0, 1e-6);
      extra.Finalize();
      const auto a3 =
          std::unique_ptr<mfem::SparseMatrix>(mfem::Add(s.A, extra));
      check(solve(s, *a3, *(ls.linear_solver), name + ", another pattern") &&
                (getNumberOfAnalyses(*(ls.linear_solver)) == size_type{3}),
            "the analysis is not computed again for an operator of another "
            "pattern");
    }
  }
  // a singular matrix, the system being inconsistent: a null row whose
  // right hand side is not null (the rows of the eliminated degrees of
  // freedom have a null right hand side: their system is consistent)
  if constexpr (!is_parallel) {
    const auto* const b = s.B.HostRead();
    const auto row = static_cast<int>(
        std::find_if(b, b + s.B.Size(), [](const auto v) { return v != 0; }) -
        b);
    auto singular = Matrix(s.A);
    std::ignore = singular.HostReadWriteI();
    std::ignore = singular.HostReadWriteJ();
    std::ignore = singular.HostReadWriteData();
    singular.EliminateRow(row, mfem::Operator::DIAG_ZERO);
    auto ls = generate(s, Parameters{});
    auto reported = false;
    try {
      ls.linear_solver->SetOperator(singular);
      auto x = mfem::Vector(s.B.Size());
      x = 0.;
      ls.linear_solver->Mult(s.B, x);
      reported = !hasConverged(*(ls.linear_solver));
    } catch (CollectiveFailure&) {
      reported = true;
    }
    if (rank == 0) {
      std::cout << getSolverName() << ", singular matrix: "
                << (reported ? "failure reported" : "not reported") << '\n';
    }
    check(reported, "singular matrix not reported");
  }
  // refusals
  auto ctx = Context{};
  const auto refused = LinearSolverFactory<is_parallel>::getFactory().generate(
      ctx, getSolverName(), s.fespace, Parameters{{"PositiveDefinite", true}});
  check(isInvalid(refused) &&
            (ctx.getErrorMessage().find("requires the option 'Symmetric'") !=
             std::string::npos),
        "PositiveDefinite without Symmetric not refused");
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
      mfem::Mesh::MakeCartesian3D(10, 10, 10, mfem::Element::HEXAHEDRON);
  auto mesh2 = mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    if (parallel) {
#if defined(MFEM_USE_MPI) && defined(MFEM_USE_MKL_CPARDISO)
      auto s = LaplaceSystem<true>(mesh, fec);
      auto s2 = LaplaceSystem<true>(mesh2, fec);
      run<true>(s, s2);
#else
      throw std::runtime_error("the Cluster Sparse Solver is not available");
#endif
    } else {
#ifdef MFEM_USE_MKL_PARDISO
      auto s = LaplaceSystem<false>(mesh, fec);
      auto s2 = LaplaceSystem<false>(mesh2, fec);
      run<false>(s, s2);
#else
      throw std::runtime_error("PARDISO is not available");
#endif
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
