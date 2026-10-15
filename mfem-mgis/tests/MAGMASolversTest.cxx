/*!
 * \file   MAGMASolversTest.cxx
 * \brief  Test of the sparse solvers of MAGMA on the Laplace equation on a
 * cube:
 *
 * - in sequential, each Krylov solver of MAGMA, without preconditioner and
 *   with each preconditioner of MAGMA, solves the system, then a second
 *   operator (the matrix times two) and an operator of another size; the
 *   incomplete Cholesky factorization reduces the number of iterations of
 *   the conjugate gradient. QMR is only run with the Jacobi and incomplete
 *   LU preconditioners: without look-ahead, it breaks down on this system
 *   without preconditioner, and it is refused with the incomplete Cholesky
 *   factorization, whose transpose MAGMA applies wrongly. The right hand
 *   side uses the memory of the device of MFEM (given to MAGMA when MFEM
 *   computes on the GPU), and the memory of the host in one case (copied);
 * - a resolution from the solution (initial guess kept by
 *   `iterative_mode`) needs no iteration; each solver stopped by the
 *   maximum number of iterations is reported as not converged; the
 *   dimensions of the Krylov space of GMRES and of the shadow space of IDR
 *   change the number of iterations; a null right hand side gives a null
 *   solution;
 * - the preconditioners of MAGMA are applied by the conjugate gradient of
 *   MFEM, in sequential and in parallel (`--parallel`, diagonal block of
 *   each process): the Jacobi preconditioner needs the number of
 *   iterations of the Jacobi smoother of MFEM, the incomplete
 *   factorizations the number of iterations of the incomplete LU
 *   factorization of MFEM;
 * - a null diagonal entry (on the first process only in parallel) makes
 *   the setup of the preconditioners fail, on all the processes;
 * - invalid options are refused.
 * \date   15/10/2026
 */

#include <array>
#include <cmath>
#include <tuple>
#include <memory>
#include <string>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/PreconditionerFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "LaplaceSystem.hxx"

using mfem_mgis::tests::LaplaceSystem;

//! \brief number of failed checks
static int failures = 0;
//! \brief rank of the process
static int rank = 0;
//! \brief parallel computations
static bool parallel = false;
//! \brief bound of the relative residuals (relative tolerance of 1e-10)
static constexpr double tolerance = 1e-9;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \return the linear solver of the factory, throws on failure
 * \param[in] s: system
 * \param[in] n: name of the solver
 * \param[in] parameters: parameters of the linear solver
 */
template <bool is_parallel>
static mfem_mgis::LinearSolverHandler generate(
    LaplaceSystem<is_parallel>& s,
    const std::string& n,
    const mfem_mgis::Parameters& parameters) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto ls = LinearSolverFactory<is_parallel>::getFactory().generate(
      ctx, n, s.fespace, parameters);
  if (isInvalid(ls)) {
    throw std::runtime_error(ctx.getErrorMessage());
  }
  return ls;
}  // end of generate

/*!
 * \return the number of iterations of the resolution of the system of a
 * matrix, with the right hand side of a system, if it converged with a
 * small residual, an empty value otherwise
 * \param[in] s: system
 * \param[in] a: matrix
 * \param[in] ls: linear solver
 * \param[in] name: name of the resolution
 * \param[in] device_rhs: if true, the right hand side uses the memory of
 * the device of MFEM, otherwise the one of the host
 */
template <bool is_parallel>
static std::optional<int> solve(
    LaplaceSystem<is_parallel>& s,
    const typename LaplaceSystem<is_parallel>::Matrix& a,
    mfem::Solver& ls,
    const std::string& name,
    const bool device_rhs = true) {
  using namespace mfem_mgis;
  auto b = mfem::Vector(s.B.Size());
  b.UseDevice(device_rhs);
  b = s.B;
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.SetOperator(a);
  ls.Mult(b, x);
  const auto r = s.getRelativeResidual(a, x);
  const auto converged = hasConverged(ls);
  const auto iterations = getNumberOfIterationsAtConvergence(ls).value_or(-1);
  if (rank == 0) {
    std::cout << name << ": " << iterations << " iterations, residual " << r
              << (converged ? "" : " (not converged)") << '\n';
  }
  if (converged && (r < tolerance)) {
    return iterations;
  }
  return {};
}  // end of solve

//! \return the parameters of a preconditioner
static mfem_mgis::Parameters makePreconditioner(const char* const n) {
  return mfem_mgis::Parameters{{"Name", n}};
}  // end of makePreconditioner

/*!
 * \return the parameters of an iterative solver
 * \param[in] maximum_number_of_iterations: maximum number of iterations
 */
static mfem_mgis::Parameters makeIterativeSolverParameters(
    const int maximum_number_of_iterations = 2000) {
  return mfem_mgis::Parameters{
      {"RelativeTolerance", 1e-10},
      {"AbsoluteTolerance", 0.},
      {"MaximumNumberOfIterations", maximum_number_of_iterations}};
}  // end of makeIterativeSolverParameters

//! \brief names of the Krylov solvers of MAGMA
static constexpr auto magma_solvers = std::array<const char*, 7>{
    "MAGMACG",  "MAGMABiCGSTAB", "MAGMAGMRES", "MAGMAIDR",
    "MAGMACGS", "MAGMAQMR",      "MAGMATFQMR"};

/*!
 * \brief checks that the setup of a preconditioner of MAGMA fails on a
 * matrix with a null diagonal entry (on the first process in parallel)
 * \param[in] s: system
 */
template <bool is_parallel>
static void checkNullDiagonal(LaplaceSystem<is_parallel>& s) {
  using namespace mfem_mgis;
  auto a = typename LaplaceSystem<is_parallel>::Matrix(s.A);
  auto rows = mfem::Array<int>{};
  if (rank == 0) {
    // a row of an unknown which is not eliminated (non null right hand
    // side)
    const auto* const b = s.B.HostRead();
    for (int i = 0; i != s.B.Size(); ++i) {
      if (std::abs(b[i]) > 0) {
        rows.Append(i);
        break;
      }
    }
  }
  if constexpr (is_parallel) {
    a.EliminateRows(rows);
  } else {
    std::ignore = a.HostReadWriteI();
    std::ignore = a.HostReadWriteJ();
    std::ignore = a.HostReadWriteData();
    for (const auto r : rows) {
      a.EliminateRow(r, mfem::Operator::DIAG_ZERO);
    }
  }
  auto reported = false;
  auto collective = false;
  try {
    auto ctx = Context{};
    auto p = PreconditionerFactory<is_parallel>::getFactory().generate(
        ctx, "MAGMAJacobi", s.fespace, Parameters{});
    if (p == nullptr) {
      throw std::runtime_error(ctx.getErrorMessage());
    }
    p->SetOperator(a);
  } catch (CollectiveFailure&) {
    reported = true;
    collective = true;
  } catch (std::exception&) {
    reported = true;
  }
  if (rank == 0) {
    std::cout << "MAGMAJacobi, null diagonal entry: "
              << (reported ? "failure reported" : "not reported") << '\n';
  }
  check(reported && ((!is_parallel) || collective),
        "the null diagonal entry is not reported (on all the processes)");
}  // end of checkNullDiagonal

//! \brief Krylov solvers of MAGMA, sequential
static void runKrylovSolvers(LaplaceSystem<false>& s,
                             LaplaceSystem<false>& s2) {
  using namespace mfem_mgis;
  auto a2 = mfem::SparseMatrix(s.A);
  a2 *= 2.;
  const auto preconditioners =
      std::array<const char*, 4>{nullptr, "MAGMAJacobi", "MAGMAILU", "MAGMAIC"};
  auto cg_iterations = std::array<int, 4>{};
  for (const auto* const n : magma_solvers) {
    for (std::size_t i = 0; i != preconditioners.size(); ++i) {
      const auto* const p = preconditioners[i];
      // QMR of MAGMA, without look-ahead, breaks down on this system
      // without preconditioner (reported as not converged), and is refused
      // with the incomplete Cholesky factorization
      if ((std::string{n} == "MAGMAQMR") &&
          ((p == nullptr) || (std::string{p} == "MAGMAIC"))) {
        continue;
      }
      auto parameters = makeIterativeSolverParameters();
      if (p != nullptr) {
        parameters.insert(throwing,
                          {{"Preconditioner", makePreconditioner(p)}});
      }
      auto ls = generate(s, n, parameters);
      const auto name =
          std::string{n} + (p != nullptr ? " + " + std::string{p} : "");
      const auto it = solve(s, s.A, *(ls.linear_solver), name);
      check(it.has_value(), name + ": resolution");
      check(solve(s, a2, *(ls.linear_solver), name + ", second operator")
                .has_value(),
            name + ": resolution with a second operator");
      check(solve(s2, s2.A, *(ls.linear_solver), name + ", smaller mesh")
                .has_value(),
            name + ": resolution on a smaller mesh");
      if ((std::string{n} == "MAGMACG") && (it.has_value())) {
        cg_iterations[i] = *it;
      }
    }
  }
  check(cg_iterations[3] < cg_iterations[0],
        "the incomplete Cholesky factorization does not reduce the number "
        "of iterations of the conjugate gradient");
  // right hand side in the memory of the host: copied
  {
    auto parameters = makeIterativeSolverParameters();
    parameters.insert(throwing,
                      {{"Preconditioner", makePreconditioner("MAGMAJacobi")}});
    auto ls = generate(s, "MAGMACG", parameters);
    check(solve(s, s.A, *(ls.linear_solver),
                "MAGMACG + MAGMAJacobi, right hand side on the host", false)
              .has_value(),
          "resolution with a right hand side on the host");
  }
  // resolution from the solution: computed with a tighter tolerance, then
  // solved again with the initial guess
  {
    const auto parameters = Parameters{
        {"RelativeTolerance", 1e-12},
        {"AbsoluteTolerance", 0.},
        {"MaximumNumberOfIterations", 2000},
        {"Preconditioner", makePreconditioner("MAGMAJacobi")}};
    auto ls = generate(s, "MAGMACG", parameters);
    auto x = mfem::Vector(s.B.Size());
    x.UseDevice(true);
    x = 0.;
    ls.linear_solver->SetOperator(s.A);
    ls.linear_solver->Mult(s.B, x);
    check(setRelativeTolerance(*(ls.linear_solver), 1e-10),
          "the relative tolerance can't be changed");
    ls.linear_solver->iterative_mode = true;
    ls.linear_solver->Mult(s.B, x);
    check(hasConverged(*(ls.linear_solver)) &&
              (getNumberOfIterationsAtConvergence(*(ls.linear_solver)) == 0),
          "the initial guess is not used");
    ls.linear_solver->iterative_mode = false;
    ls.linear_solver->Mult(s.B, x);
    check(
        getNumberOfIterationsAtConvergence(*(ls.linear_solver)).value_or(0) > 0,
        "the initial guess is used while iterative_mode is false");
    // null right hand side, the initial guess not being null
    auto b = mfem::Vector(s.B.Size());
    b.UseDevice(true);
    b = 0.;
    ls.linear_solver->iterative_mode = true;
    ls.linear_solver->Mult(b, x);
    check(hasConverged(*(ls.linear_solver)) && !(x.Normlinf() > 0),
          "null right hand side: the solution is not null");
  }
  // non convergence of each solver
  for (const auto* const n : magma_solvers) {
    auto ls = generate(s, n, makeIterativeSolverParameters(2));
    check(!solve(s, s.A, *(ls.linear_solver), std::string{n} + ", 2 iterations")
                  .has_value() &&
              !hasConverged(*(ls.linear_solver)),
          std::string{n} +
              ": resolution stopped by the maximum number of iterations");
  }
  // dimensions of the Krylov space of GMRES and of the shadow space of IDR
  {
    auto p5 = makeIterativeSolverParameters();
    p5.insert(throwing, {{"KDim", 5}});
    auto ls50 = generate(s, "MAGMAGMRES", makeIterativeSolverParameters());
    auto ls5 = generate(s, "MAGMAGMRES", p5);
    const auto i50 = solve(s, s.A, *(ls50.linear_solver), "MAGMAGMRES(50)");
    const auto i5 = solve(s, s.A, *(ls5.linear_solver), "MAGMAGMRES(5)");
    check(i50.has_value() && i5.has_value() && (*i5 > *i50),
          "the dimension of the Krylov space of GMRES has no effect");
    auto p4 = makeIterativeSolverParameters();
    p4.insert(throwing, {{"ShadowSpaceDimension", 4}});
    auto ls1 = generate(s, "MAGMAIDR", makeIterativeSolverParameters());
    auto ls4 = generate(s, "MAGMAIDR", p4);
    const auto j1 = solve(s, s.A, *(ls1.linear_solver), "MAGMAIDR(1)");
    const auto j4 = solve(s, s.A, *(ls4.linear_solver), "MAGMAIDR(4)");
    check(j1.has_value() && j4.has_value() && (*j1 != *j4),
          "the dimension of the shadow space of IDR has no effect");
  }
  // refusals
  const auto refused = [&s](const char* const n, const Parameters& p,
                            const std::string& msg) {
    auto ctx = Context{};
    const auto ls =
        LinearSolverFactory<false>::getFactory().generate(ctx, n, s.fespace, p);
    const auto error = std::string{ctx.getErrorMessage()};
    check(isInvalid(ls) && (error.find(msg) != std::string::npos),
          std::string{n} + ": invalid options not refused ('" + error + "')");
  };
  refused("MAGMAIDR", Parameters{{"ShadowSpaceDimension", 50}},
          "shadow space of dimension 50");
  refused("MAGMACG", Parameters{{"RelativeTolerance", -1.}},
          "negative tolerance");
  refused("MAGMACG", Parameters{{"KDim", 10}}, "KDim");
  refused("MAGMACG",
          Parameters{{"Preconditioner", makePreconditioner("DSmoother")}},
          "only accepts the preconditioners of MAGMA");
  refused("MAGMAQMR",
          Parameters{{"Preconditioner", makePreconditioner("MAGMAIC")}},
          "does not accept the preconditioner 'MAGMAIC'");
}  // end of runKrylovSolvers

//! \brief preconditioners of MAGMA applied by the conjugate gradient of MFEM
template <bool is_parallel>
static void runPreconditioners(LaplaceSystem<is_parallel>& s) {
  using namespace mfem_mgis;
  const auto preconditioners = std::array<const char*, 6>{
      nullptr, "MAGMAJacobi", "MAGMAILU", "MAGMAIC", "DSmoother", "BlockILU"};
  auto iterations = std::array<int, 6>{};
  for (std::size_t i = 0; i != preconditioners.size(); ++i) {
    const auto* const p = preconditioners[i];
    auto parameters = makeIterativeSolverParameters();
    if (p != nullptr) {
      // the references: the Jacobi smoother and the incomplete LU
      // factorization without fill-in of MFEM, on the same matrix (the
      // diagonal block of each process in parallel)
      const auto prec =
          std::string{p} == "BlockILU"
              ? Parameters{{"Name", "BlockILU"},
                           {"Options", Parameters{{"BlockSize", 1},
                                                  {"Reordering", "None"}}}}
              : makePreconditioner(p);
      parameters.insert(throwing, {{"Preconditioner", prec}});
    }
    auto ls = generate(s, "CGSolver", parameters);
    const auto name =
        std::string{"CGSolver"} + (p != nullptr ? " + " + std::string{p} : "");
    const auto it = solve(s, s.A, *(ls.linear_solver), name);
    check(it.has_value(), name + ": resolution");
    iterations[i] = it.value_or(-10);
  }
  check(std::abs(iterations[1] - iterations[4]) <= 1,
        "the Jacobi preconditioners of MAGMA and MFEM differ");
  check((std::abs(iterations[2] - iterations[5]) <= 1) &&
            (std::abs(iterations[3] - iterations[5]) <= 1),
        "the incomplete factorizations of MAGMA and MFEM differ");
  if constexpr (!is_parallel) {
    check((iterations[2] < iterations[0]) && (iterations[3] < iterations[0]),
          "the incomplete factorizations of MAGMA do not reduce the number "
          "of iterations of the conjugate gradient of MFEM");
  }
  checkNullDiagonal(s);
}  // end of runPreconditioners

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
#ifdef MFEM_USE_MPI
      auto s = LaplaceSystem<true>(mesh, fec);
      runPreconditioners<true>(s);
#else
      throw std::runtime_error("MFEM is not built with MPI");
#endif
    } else {
      auto s = LaplaceSystem<false>(mesh, fec);
      auto s2 = LaplaceSystem<false>(mesh2, fec);
      runKrylovSolvers(s, s2);
      runPreconditioners<false>(s);
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
