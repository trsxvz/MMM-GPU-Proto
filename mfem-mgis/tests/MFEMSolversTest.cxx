/*!
 * \file   MFEMSolversTest.cxx
 * \brief  Test of the preconditioners of the sparse matrices of MFEM
 * (Gauss-Seidel and Jacobi smoothers, block ILU factorization), of the
 * direct solvers of SuiteSparse used as preconditioners or solvers, and of
 * the orthogonalizing wrapper of a preconditioner, on the Laplace equation
 * on a cube, in sequential and in parallel (the preconditioners being
 * applied to the diagonal block of each process), with the `--device`
 * option on a GPU (the matrices and the vectors being copied on the host):
 *
 * - each type of the smoothers, with one or two iterations, and the
 *   reorderings of the block ILU factorization precondition the GMRES
 *   solver of MFEM;
 * - the direct solvers of SuiteSparse give the solution in one iteration in
 *   sequential, and precondition the diagonal blocks in parallel;
 * - the KLU direct solver (sequential), which reports a singular matrix;
 * - each preconditioner is set up again for a second operator of the same
 *   pattern (the matrix times two);
 * - the orthogonalizing wrapper solves the Laplace equation with null flux
 *   boundary conditions, whose kernel is the constant vector, its solution
 *   being orthogonal to the constant vector;
 * - an operator applied without assembly and invalid options are refused.
 * \date   15/10/2026
 */

#include <cmath>
#include <string>
#include <tuple>
#include <type_traits>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
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
 * \return true if a resolution converged with a relative residual lower
 * than 100 times the tolerance
 * \param[in] r: resolution
 */
static bool isSolved(const Resolution& r) {
  return r.converged && (r.residual < 100 * tolerance);
}  // end of isSolved

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

/*!
 * \return the resolution of the system of a matrix, with the right hand
 * side of a system, by a linear solver
 * \param[in] s: system
 * \param[in] a: matrix
 * \param[in] ls: linear solver
 * \param[in] name: name of the resolution
 * \param[out] solution: solution, if not null
 */
template <bool is_parallel>
static Resolution solve(LaplaceSystem<is_parallel>& s,
                        const typename LaplaceSystem<is_parallel>::Matrix& a,
                        mfem::Solver& ls,
                        const std::string& name,
                        mfem::Vector* const solution = nullptr) {
  using namespace mfem_mgis;
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.SetOperator(a);
  ls.Mult(s.B, x);
  auto r = Resolution{};
  r.converged = hasConverged(ls);
  r.iterations = getNumberOfIterationsAtConvergence(ls).value_or(-1);
  r.residual = s.getRelativeResidual(a, x);
  if (rank == 0) {
    std::cout << name << ": " << r.iterations << " iterations, residual "
              << r.residual << (r.converged ? "" : " (not converged)") << '\n';
  }
  if (solution != nullptr) {
    *solution = x;
  }
  return r;
}  // end of solve

/*!
 * \return the resolution of a system by a linear solver
 * \param[in] s: system
 * \param[in] ls: linear solver
 * \param[in] name: name of the resolution
 */
template <bool is_parallel>
static Resolution solve(LaplaceSystem<is_parallel>& s,
                        mfem::Solver& ls,
                        const std::string& name) {
  return solve(s, s.A, ls, name);
}  // end of solve

/*!
 * \return true if a linear solver solves the system, then the system of a
 * second operator of the same pattern
 * \param[in] s: system
 * \param[in] a2: second operator
 * \param[in] ls: linear solver
 * \param[in] name: name of the resolution
 */
template <bool is_parallel>
static bool solveTwice(LaplaceSystem<is_parallel>& s,
                       const typename LaplaceSystem<is_parallel>::Matrix& a2,
                       mfem::Solver& ls,
                       const std::string& name) {
  const auto r = solve(s, ls, name);
  const auto r2 = solve(s, a2, ls, name + ", second operator");
  return isSolved(r) && isSolved(r2);
}  // end of solveTwice

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
      {"MaximumNumberOfIterations", 1000},
      {"KDim", 100},
      {"Preconditioner", Parameters{{"Name", p}, {"Options", o}}}};
}  // end of getGMRESParameters

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

template <bool is_parallel>
static void run(LaplaceSystem<is_parallel>& s,
                LaplaceSystem<is_parallel>& neumann) {
  using namespace mfem_mgis;
  using Matrix = typename LaplaceSystem<is_parallel>::Matrix;
  // a second operator of the same pattern: the matrix times two
  auto a2 = Matrix(s.A);
  a2 *= 2.;
  // smoothers and incomplete factorization of the sparse matrices
  for (const auto& t : {"Symmetric", "Forward", "Backward"}) {
    auto gs = generate(
        s, "GMRESSolver",
        getGMRESParameters("GSSmoother",
                           Parameters{{"Type", t}, {"NumberOfIterations", 2}}));
    check(solveTwice(s, a2, *(gs.linear_solver),
                     std::string{"GMRESSolver, GSSmoother "} + t),
          std::string{"GSSmoother "} + t);
  }
  for (const auto& t : {"Jacobi", "l1Jacobi"}) {
    for (const auto n : {1, 2}) {
      auto d =
          generate(s, "GMRESSolver",
                   getGMRESParameters("DSmoother",
                                      Parameters{{"Type", t},
                                                 {"Damping", 0.8},
                                                 {"NumberOfIterations", n}}));
      const auto name = std::string{"GMRESSolver, DSmoother "} + t + ", " +
                        std::to_string(n) + " iteration(s)";
      check(solveTwice(s, a2, *(d.linear_solver), name), name);
    }
  }
  for (const auto& r : {"MinimumDiscardedFill", "None"}) {
    auto ilu =
        generate(s, "GMRESSolver",
                 getGMRESParameters("BlockILU", Parameters{{"Reordering", r}}));
    check(solveTwice(s, a2, *(ilu.linear_solver),
                     std::string{"GMRESSolver, BlockILU "} + r),
          std::string{"BlockILU "} + r);
  }
#ifdef MFEM_MGIS_HAVE_SUITESPARSE
  // the direct solvers of SuiteSparse
  for (const auto& n : {"UMFPack", "KLU"}) {
    auto direct = generate(s, "GMRESSolver", getGMRESParameters(n));
    const auto name = std::string{"GMRESSolver, "} + n;
    const auto r = solve(s, *(direct.linear_solver), name);
    const auto r2 =
        solve(s, a2, *(direct.linear_solver), name + ", second operator");
    check(isSolved(r) && isSolved(r2) &&
              (is_parallel || ((r.iterations <= 1) && (r2.iterations <= 1))),
          name);
  }
  if constexpr (!is_parallel) {
    auto klu = generate(s, "KLUSolver", Parameters{});
    check((solve(s, *(klu.linear_solver), "KLUSolver").residual < 1e-12) &&
              (solve(s, a2, *(klu.linear_solver), "KLUSolver, second operator")
                   .residual < 1e-12),
          "KLUSolver: resolution");
    // a singular matrix (a null row) is reported by all the processes
    // `EliminateRow` modifies the arrays of the host, made valid before
    auto singular = Matrix(s.A);
    std::ignore = singular.HostReadWriteI();
    std::ignore = singular.HostReadWriteJ();
    std::ignore = singular.HostReadWriteData();
    singular.EliminateRow(0, mfem::Operator::DIAG_ZERO);
    auto reported = false;
    try {
      klu.linear_solver->SetOperator(singular);
    } catch (CollectiveFailure&) {
      reported = true;
    }
    check(reported, "KLUSolver: singular matrix not reported");
    check(getGenerationError(s, "KLUSolver", Parameters{{"Verbosity", 1}})
                  .find("Verbosity") != std::string::npos,
          "KLUSolver: unknown option not refused");
  }
#endif /* MFEM_MGIS_HAVE_SUITESPARSE */
  // the orthogonalizing wrapper, for an operator whose kernel is the
  // constant vector
  auto ortho = generate(
      neumann, "CGSolver",
      Parameters{
          {"RelativeTolerance", tolerance},
          {"MaximumNumberOfIterations", 1000},
          {"Preconditioner",
           Parameters{{"Name", "OrthoSolver"},
                      {"Options",
                       Parameters{{"Preconditioner",
                                   Parameters{{"Name", "GSSmoother"}}}}}}}});
  auto x = mfem::Vector{};
  const auto ro = solve(neumann, neumann.A, *(ortho.linear_solver),
                        "CGSolver, OrthoSolver, null flux", &x);
  const auto mean = neumann.getRelativeMean(x);
  // without the wrapper, the solution is not orthogonal to the constant
  // vector
  auto plain = generate(
      neumann, "CGSolver",
      Parameters{{"RelativeTolerance", tolerance},
                 {"MaximumNumberOfIterations", 1000},
                 {"Preconditioner", Parameters{{"Name", "GSSmoother"}}}});
  auto x_plain = mfem::Vector{};
  std::ignore = solve(neumann, neumann.A, *(plain.linear_solver),
                      "CGSolver, GSSmoother, null flux", &x_plain);
  const auto mean_plain = neumann.getRelativeMean(x_plain);
  if (rank == 0) {
    std::cout << "relative mean of the solution: " << mean
              << " with OrthoSolver, " << mean_plain << " without\n";
  }
  check(isSolved(ro) && (mean < 100 * tolerance),
        "CGSolver, OrthoSolver: resolution with null flux boundary "
        "conditions");
  // refusals
  auto gs = generate(s, "GMRESSolver", getGMRESParameters("GSSmoother"));
  auto matrix_free = MatrixFreeOperator(s.A);
  auto raised = false;
  try {
    gs.linear_solver->SetOperator(matrix_free);
  } catch (std::exception&) {
    raised = true;
  }
  check(raised, "GSSmoother: operator applied without assembly not refused");
  const auto refused = [&s](const std::string& p, const Parameters& o,
                            const std::string& msg) {
    return getGenerationError(s, "GMRESSolver", getGMRESParameters(p, o))
               .find(msg) != std::string::npos;
  };
  check(refused("GSSmoother", Parameters{{"Type", "Diagonal"}}, "invalid type"),
        "GSSmoother: invalid type not refused");
  check(refused("DSmoother", Parameters{{"Type", "SOR"}}, "invalid type"),
        "DSmoother: invalid type not refused");
  check(refused("BlockILU", Parameters{{"Reordering", "RCM"}},
                "invalid reordering"),
        "BlockILU: invalid reordering not refused");
  check(refused("DSmoother", Parameters{{"Damping", 0.}},
                "invalid damping factor"),
        "DSmoother: null damping factor not refused");
  check(
      refused("OrthoSolver", Parameters{}, "'Preconditioner' is not declared"),
      "OrthoSolver: missing preconditioner not refused");
  check(refused("OrthoSolver",
                Parameters{{"Preconditioner", Parameters{{"Name", "None"}}}},
                "requires a preconditioner"),
        "OrthoSolver: no preconditioner not refused");
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
  auto mesh = mfem::Mesh::MakeCartesian3D(8, 8, 8, mfem::Element::HEXAHEDRON);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    if (parallel) {
#ifdef MFEM_USE_MPI
      auto s = LaplaceSystem<true>(mesh, fec);
      auto neumann = LaplaceSystem<true>(mesh, fec, false);
      run<true>(s, neumann);
#else  /* MFEM_USE_MPI */
      throw std::runtime_error("parallel computations are not supported");
#endif /* MFEM_USE_MPI */
    } else {
      auto s = LaplaceSystem<false>(mesh, fec);
      auto neumann = LaplaceSystem<false>(mesh, fec, false);
      run<false>(s, neumann);
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
