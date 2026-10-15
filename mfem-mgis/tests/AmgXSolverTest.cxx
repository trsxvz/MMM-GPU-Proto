/*!
 * \file   AmgXSolverTest.cxx
 * \brief  Test of the solvers of the AmgX library built by the factory of
 * linear solvers (`AmgX` preconditioner and linear solver, `AmgXPCG` linear
 * solver), on the Laplace equation on a cube, the essential degrees of
 * freedom being eliminated with a diagonal set to one. AmgX computes on the
 * GPU whatever the device used by MFEM. With the `--parallel` option, the
 * computations are parallel. With the `--no-gpu` option, no GPU is expected
 * to be visible: the solvers of AmgX must be refused.
 * \date   10/10/2026
 */

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <functional>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/MPI.hxx"

//! \brief number of failed checks
static int failures = 0;
//! \brief if true, the computations are parallel
static bool parallel = false;
//! \brief rank of the process
static int rank = 0;
//! \brief if true, no GPU is expected
static bool no_gpu = false;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

//! \brief configuration of AmgX: one V-cycle of the aggregation multigrid
static const char* const aggregation =
    "{\"config_version\": 2, \"solver\": {\"solver\": \"AMG\", "
    "\"scope\": \"main\", \"algorithm\": \"AGGREGATION\", "
    "\"selector\": \"SIZE_2\", \"smoother\": {\"scope\": \"jacobi\", "
    "\"solver\": \"BLOCK_JACOBI\", \"relaxation_factor\": 0.8}, "
    "\"presweeps\": 0, \"postsweeps\": 3, \"coarse_solver\": \"NOSOLVER\", "
    "\"max_iters\": 1, \"max_levels\": 50, \"cycle\": \"V\"}}";

//! \brief result of a resolution
struct Resolution {
  //! \brief solution
  mfem::Vector x;
  //! \brief if true, the resolution succeeded
  bool converged = false;
  //! \brief number of iterations, if given by the solver
  std::optional<int> iterations;
  //! \brief norm of the residual relative to the norm of the right hand side
  double residual = 1;
};

/*!
 * \return the scalar product of two vectors
 * \param[in] fespace: finite element space
 * \param[in] a: first vector
 * \param[in] b: second vector
 */
template <bool is_parallel>
static double dot(mfem_mgis::FiniteElementSpace<is_parallel>& fespace,
                  const mfem::Vector& a,
                  const mfem::Vector& b) {
  if constexpr (is_parallel) {
    return mfem::InnerProduct(fespace.GetComm(), a, b);
  } else {
    static_cast<void>(fespace);
    return mfem::InnerProduct(a, b);
  }
}  // end of dot

/*!
 * \return the norm of a vector relative to the norm of a reference vector,
 * computed on copies multiplied by the same power of two (the scalar
 * products of very small or very large vectors underflow or overflow)
 * \param[in] fespace: finite element space
 * \param[in] v: vector
 * \param[in] r: reference vector
 */
template <bool is_parallel>
static double getRelativeNorm(
    mfem_mgis::FiniteElementSpace<is_parallel>& fespace,
    const mfem::Vector& v,
    const mfem::Vector& r) {
  auto m = r.Normlinf();
  if constexpr (is_parallel) {
    MPI_Allreduce(MPI_IN_PLACE, &m, 1, mfem_mgis::mpi_type<mfem::real_t>,
                  MPI_MAX, fespace.GetComm());
  }
  const auto f = (m > 0) ? std::ldexp(1., -std::ilogb(m)) : 1.;
  auto sv = mfem::Vector(v.Size());
  sv.UseDevice(true);
  sv.Set(f, v);
  auto sr = mfem::Vector(r.Size());
  sr.UseDevice(true);
  sr.Set(f, r);
  return std::sqrt(dot<is_parallel>(fespace, sv, sv) /
                   dot<is_parallel>(fespace, sr, sr));
}  // end of getRelativeNorm

/*!
 * \return the resolution of a linear system by a linear solver
 * \param[in] fespace: finite element space
 * \param[in, out] ls: linear solver
 * \param[in] A: operator
 * \param[in] B: right hand side
 * \param[in] x0: initial value of the solution, ignored by the solvers of
 * the factory (zero initial guess)
 */
template <bool is_parallel>
static Resolution solve(mfem_mgis::FiniteElementSpace<is_parallel>& fespace,
                        mfem::Solver& ls,
                        const mfem::Operator& A,
                        const mfem::Vector& B,
                        const double x0 = 0) {
  using namespace mfem_mgis;
  auto r = Resolution{};
  r.x.SetSize(B.Size());
  r.x.UseDevice(true);
  r.x = x0;
  ls.SetOperator(A);
  ls.Mult(B, r.x);
  auto e = mfem::Vector(B.Size());
  e.UseDevice(true);
  A.Mult(r.x, e);
  e -= B;
  r.converged = hasConverged(ls);
  r.iterations = getNumberOfIterationsAtConvergence(ls);
  r.residual = getRelativeNorm<is_parallel>(fespace, e, B);
  return r;
}  // end of solve

/*!
 * \return a linear solver of the factory, throws on failure (the solvers
 * alive are destroyed before the exit)
 * \param[in] fespace: finite element space
 * \param[in] name: name of the linear solver
 * \param[in] parameters: options of the linear solver
 */
template <bool is_parallel>
static mfem_mgis::LinearSolverHandler generate(
    mfem_mgis::FiniteElementSpace<is_parallel>& fespace,
    const std::string& name,
    const mfem_mgis::Parameters& parameters) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto ls = LinearSolverFactory<is_parallel>::getFactory().generate(
      ctx, name, fespace, parameters);
  if (isInvalid(ls)) {
    throw std::runtime_error(ctx.getErrorMessage());
  }
  return ls;
}  // end of generate

/*!
 * \return the conjugate gradient of MFEM with the given preconditioner
 * \param[in] fespace: finite element space
 * \param[in] preconditioner: preconditioner
 */
template <bool is_parallel>
static mfem_mgis::LinearSolverHandler makeConjugateGradient(
    mfem_mgis::FiniteElementSpace<is_parallel>& fespace,
    const mfem_mgis::Parameters& preconditioner) {
  using namespace mfem_mgis;
  return generate<is_parallel>(fespace, "CGSolver",
                               Parameters{{"RelativeTolerance", 1e-10},
                                          {"AbsoluteTolerance", 0.},
                                          {"MaximumNumberOfIterations", 5000},
                                          {"Preconditioner", preconditioner}});
}  // end of makeConjugateGradient

/*!
 * \return if two numbers of iterations are close: the setup of AmgX on the
 * GPU may not be deterministic
 */
static bool close(const std::optional<int> a, const std::optional<int> b) {
  return a.has_value() && b.has_value() && (std::abs(*a - *b) <= *b / 20 + 1);
}  // end of close

/*!
 * \brief checks the solvers of AmgX on a linear system
 * \param[in] fespace: finite element space
 * \param[in] A: assembled operator
 * \param[in] B: right hand side
 * \param[in] release: function destroying the other solvers of AmgX of the
 * process, called before the destruction of all the solvers
 */
template <bool is_parallel>
static void checkAmgX(mfem_mgis::FiniteElementSpace<is_parallel>& fespace,
                      const mfem::Operator& A,
                      const mfem::Vector& B,
                      const std::function<void()>& release = {}) {
  using namespace mfem_mgis;
  auto print = [](const std::string& name, const Resolution& r) {
    if (rank == 0) {
      std::cout << name << ": ";
      if (r.iterations.has_value()) {
        std::cout << *(r.iterations) << " iterations, ";
      }
      std::cout << "residual " << r.residual
                << (r.converged ? "" : " (not converged)") << '\n';
    }
  };
  // two solvers, the first one being destroyed before the second one is
  // used (in sequential, the first solvers of AmgX of the process)
  {
    auto f1 = generate<is_parallel>(fespace, "AmgXPCG",
                                    Parameters{{"RelativeTolerance", 1e-10}});
    auto f2 = generate<is_parallel>(fespace, "AmgXPCG",
                                    Parameters{{"RelativeTolerance", 1e-10}});
    f1 = LinearSolverHandler{};
    const auto r = solve<is_parallel>(fespace, *(f2.linear_solver), A, B);
    check(r.converged && (r.residual < 1e-9),
          "AmgXPCG: resolution after the destruction of the first solver");
  }
  // reference: conjugate gradient preconditioned by a Jacobi smoother
  auto jacobi = makeConjugateGradient<is_parallel>(
      fespace, Parameters{{"Name", "OperatorJacobiSmoother"}});
  const auto rj = solve<is_parallel>(fespace, *(jacobi.linear_solver), A, B);
  print("CG, Jacobi", rj);
  check(rj.converged && rj.iterations.has_value(), "CG, Jacobi: resolution");
  // conjugate gradient of MFEM preconditioned by AmgX, default configuration
  // of MFEM, set up twice
  auto cg =
      makeConjugateGradient<is_parallel>(fespace, Parameters{{"Name", "AmgX"}});
  const auto r1 = solve<is_parallel>(fespace, *(cg.linear_solver), A, B);
  const auto r2 = solve<is_parallel>(fespace, *(cg.linear_solver), A, B);
  print("CG, AmgX", r1);
  print("CG, AmgX, set up again", r2);
  // the conjugate gradient of MFEM checks the norm of the preconditioned
  // residual: the norm of the residual is checked with a margin
  for (const auto& r : {r1, r2}) {
    check(r.converged && (r.residual < 1e-8), "CG, AmgX: resolution");
  }
  check(close(r2.iterations, r1.iterations),
        "CG, AmgX: different numbers of iterations after a new setup");
  check(rj.iterations.has_value() && r1.iterations.has_value() &&
            (2 * *(r1.iterations) < *(rj.iterations)),
        "CG, AmgX: the multigrid does not reduce the number of iterations");
  // configuration given as a string, then by a file
  auto cgs = makeConjugateGradient<is_parallel>(
      fespace, Parameters{{"Name", "AmgX"},
                          {"Options", Parameters{{"Configuration",
                                                  std::string{aggregation}}}}});
  const auto rs = solve<is_parallel>(fespace, *(cgs.linear_solver), A, B);
  print("CG, AmgX, aggregation (string)", rs);
  const auto file = "AmgXSolverTest-" + std::string{is_parallel ? "p" : "s"} +
                    std::to_string(rank) + ".json";
  std::ofstream{file} << aggregation << '\n';
  auto cgf = makeConjugateGradient<is_parallel>(
      fespace,
      Parameters{{"Name", "AmgX"},
                 {"Options", Parameters{{"ConfigurationFile", file}}}});
  const auto rf = solve<is_parallel>(fespace, *(cgf.linear_solver), A, B);
  print("CG, AmgX, aggregation (file)", rf);
  for (const auto& r : {rs, rf}) {
    check(r.converged && (r.residual < 1e-8),
          "CG, AmgX, aggregation: resolution");
  }
  check(close(rf.iterations, rs.iterations),
        "CG, AmgX: different numbers of iterations for the same "
        "configuration given as a string and by a file");
  check(!close(rs.iterations, r1.iterations),
        "CG, AmgX: same numbers of iterations for the aggregation and the "
        "default configurations (configuration ignored ?)");
  // one V-cycle as a linear solver, as `HypreBoomerAMG`: zero initial
  // guess, no number of iterations
  auto amg = generate<is_parallel>(fespace, "AmgX", Parameters{});
  const auto rc0 = solve<is_parallel>(fespace, *(amg.linear_solver), A, B);
  const auto rc1 = solve<is_parallel>(fespace, *(amg.linear_solver), A, B, 1e3);
  print("AmgX, one cycle", rc0);
  check(rc0.converged && (!rc0.iterations.has_value()),
        "AmgX: convergence and number of iterations of one cycle");
  check((rc0.residual < 0.9) &&
            (std::abs(rc1.residual - rc0.residual) <= 1e-6 * rc0.residual),
        "AmgX: one cycle does not reduce the residual or depends on the "
        "initial value of the solution");
  // statistics printed by AmgX
  {
    auto v = generate<is_parallel>(fespace, "AmgX",
                                   Parameters{{"VerbosityLevel", 1}});
    const auto r = solve<is_parallel>(fespace, *(v.linear_solver), A, B);
    check(std::abs(r.residual - rc0.residual) <= 1e-6 * rc0.residual,
          "AmgX: one cycle with the statistics of AmgX");
    // a small right hand side: the residual computed for the statistics
    // must not stop the cycle
    auto Bs = mfem::Vector(B.Size());
    Bs.UseDevice(true);
    Bs.Set(1e-20, B);
    const auto rsmall = solve<is_parallel>(fespace, *(v.linear_solver), A, Bs);
    check(std::abs(rsmall.residual - rc0.residual) <= 1e-6 * rc0.residual,
          "AmgX: one cycle with the statistics of AmgX, small right hand "
          "side");
  }
  // parameters of AmgX completing the default configuration: two cycles
  // reduce the residual more than one
  auto amg2 = generate<is_parallel>(
      fespace, "AmgX",
      Parameters{{"AmgXParameters", Parameters{{"max_iters", 2}}}});
  const auto rc2 = solve<is_parallel>(fespace, *(amg2.linear_solver), A, B);
  print("AmgX, two cycles", rc2);
  check(rc2.residual < 0.9 * rc0.residual,
        "AmgX: the parameters of AmgX are not taken into account (two "
        "cycles)");
  // conjugate gradient of AmgX
  auto pcg =
      generate<is_parallel>(fespace, "AmgXPCG",
                            Parameters{{"RelativeTolerance", 1e-10},
                                       {"MaximumNumberOfIterations", 200}});
  const auto rp = solve<is_parallel>(fespace, *(pcg.linear_solver), A, B, 1e3);
  print("AmgXPCG", rp);
  check(rp.converged && (rp.residual < 1e-9) && rp.iterations.has_value() &&
            (*(rp.iterations) > 2),
        "AmgXPCG: resolution");
  check(rp.iterations.has_value() && rj.iterations.has_value() &&
            (2 * *(rp.iterations) < *(rj.iterations)),
        "AmgXPCG: the multigrid does not reduce the number of iterations");
  // aggregation multigrid of the default configurations: conjugate
  // gradient of AmgX, and AmgX as a preconditioner of the conjugate gradient
  // of MFEM
  {
    auto apcg =
        generate<is_parallel>(fespace, "AmgXPCG",
                              Parameters{{"RelativeTolerance", 1e-10},
                                         {"MaximumNumberOfIterations", 200},
                                         {"Multigrid", "Aggregation"}});
    const auto ra = solve<is_parallel>(fespace, *(apcg.linear_solver), A, B);
    print("AmgXPCG, aggregation", ra);
    check(ra.converged && (ra.residual < 1e-9) && ra.iterations.has_value() &&
              rj.iterations.has_value() &&
              (*(ra.iterations) < *(rj.iterations)),
          "AmgXPCG: resolution with the aggregation multigrid");
    auto acg = generate<is_parallel>(
        fespace, "CGSolver",
        Parameters{{"RelativeTolerance", 1e-10},
                   {"AbsoluteTolerance", 0.},
                   {"MaximumNumberOfIterations", 200},
                   {"Preconditioner",
                    Parameters{{"Name", "AmgX"},
                               {"Options",
                                Parameters{{"Multigrid", "Aggregation"}}}}}});
    const auto rcg = solve<is_parallel>(fespace, *(acg.linear_solver), A, B);
    print("CGSolver, AmgX, aggregation", rcg);
    check(rcg.converged && (rcg.residual < 1e-9),
          "CGSolver: resolution with the aggregation multigrid of AmgX");
    // the option gives the configuration written above, which differs from
    // the default one (the option is not ignored)
    check(close(rcg.iterations, rs.iterations),
          "CGSolver: the aggregation multigrid of AmgX differs from its "
          "configuration given as a string");
  }
  // statistics printed by AmgX: same resolution
  {
    auto v =
        generate<is_parallel>(fespace, "AmgXPCG",
                              Parameters{{"RelativeTolerance", 1e-10},
                                         {"MaximumNumberOfIterations", 200},
                                         {"VerbosityLevel", 1}});
    const auto r = solve<is_parallel>(fespace, *(v.linear_solver), A, B);
    check(r.converged && (r.residual < 1e-9) &&
              close(r.iterations, rp.iterations),
          "AmgXPCG: resolution with the statistics of AmgX");
  }
  // default options (relative tolerance 1e-12, 150 iterations)
  {
    auto d = generate<is_parallel>(fespace, "AmgXPCG", Parameters{});
    const auto r = solve<is_parallel>(fespace, *(d.linear_solver), A, B);
    print("AmgXPCG, default options", r);
    check(r.converged && (r.residual < 1e-11) && r.iterations.has_value() &&
              rp.iterations.has_value() && (*(r.iterations) > *(rp.iterations)),
          "AmgXPCG: resolution with the default options");
  }
  // too few iterations: not converged
  auto pcg2 =
      generate<is_parallel>(fespace, "AmgXPCG",
                            Parameters{{"RelativeTolerance", 1e-12},
                                       {"MaximumNumberOfIterations", 2}});
  const auto rp2 = solve<is_parallel>(fespace, *(pcg2.linear_solver), A, B);
  print("AmgXPCG, two iterations", rp2);
  check((!rp2.converged) && (rp2.iterations == 2) && (rp2.residual > 1e-9),
        "AmgXPCG: resolution not converged in two iterations");
  // two cycles of the multigrid per iteration: fewer iterations
  auto pcgc = generate<is_parallel>(
      fespace, "AmgXPCG",
      Parameters{
          {"RelativeTolerance", 1e-10},
          {"MaximumNumberOfIterations", 200},
          {"AmgXParameters",
           Parameters{{"max_iters", 2}, {"smoother:relaxation_factor", 0.9}}}});
  const auto rpc = solve<is_parallel>(fespace, *(pcgc.linear_solver), A, B);
  print("AmgXPCG, two cycles", rpc);
  check(rpc.converged && (rpc.residual < 1e-9) && rp.iterations.has_value() &&
            rpc.iterations.has_value() &&
            (*(rpc.iterations) < *(rp.iterations)),
        "AmgXPCG: the parameters of AmgX are not taken into account (two "
        "cycles)");
  // absolute tolerance only
  const auto nB = std::sqrt(dot<is_parallel>(fespace, B, B));
  auto pcga =
      generate<is_parallel>(fespace, "AmgXPCG",
                            Parameters{{"RelativeTolerance", 0.},
                                       {"AbsoluteTolerance", 1e-5 * nB},
                                       {"MaximumNumberOfIterations", 200}});
  const auto rpa = solve<is_parallel>(fespace, *(pcga.linear_solver), A, B);
  print("AmgXPCG, absolute tolerance", rpa);
  check(rpa.converged && rpa.iterations.has_value() &&
            rp.iterations.has_value() &&
            (*(rpa.iterations) < *(rp.iterations)) && (rpa.residual <= 1e-5),
        "AmgXPCG: resolution with an absolute tolerance");
  // absolute tolerance only, the right hand side being multiplied by 1000:
  // more iterations
  auto B3 = mfem::Vector(B.Size());
  B3.UseDevice(true);
  B3.Set(1e3, B);
  const auto rpa3 = solve<is_parallel>(fespace, *(pcga.linear_solver), A, B3);
  print("AmgXPCG, absolute tolerance, right hand side multiplied by 1000",
        rpa3);
  check(rpa3.converged && rpa3.iterations.has_value() &&
            rpa.iterations.has_value() &&
            (*(rpa3.iterations) > *(rpa.iterations)) && (rpa3.residual <= 1e-8),
        "AmgXPCG: absolute tolerance, right hand side multiplied by 1000");
  // relative tolerance of 1e-4, alone or with a negligible absolute
  // tolerance: same resolution
  auto pcgr =
      generate<is_parallel>(fespace, "AmgXPCG",
                            Parameters{{"RelativeTolerance", 1e-4},
                                       {"MaximumNumberOfIterations", 200}});
  const auto rpr = solve<is_parallel>(fespace, *(pcgr.linear_solver), A, B);
  print("AmgXPCG, relative tolerance 1e-4", rpr);
  check(rpr.converged && rpr.iterations.has_value() &&
            rp.iterations.has_value() &&
            (*(rpr.iterations) < *(rp.iterations)) && (rpr.residual <= 1e-4) &&
            (rpr.residual > 1e-6),
        "AmgXPCG: resolution with a relative tolerance of 1e-4");
  auto pcgra =
      generate<is_parallel>(fespace, "AmgXPCG",
                            Parameters{{"RelativeTolerance", 1e-4},
                                       {"AbsoluteTolerance", 1e-20},
                                       {"MaximumNumberOfIterations", 200}});
  const auto rpra = solve<is_parallel>(fespace, *(pcgra.linear_solver), A, B);
  print("AmgXPCG, relative tolerance 1e-4, absolute tolerance 1e-20", rpra);
  check(rpra.converged && (rpra.iterations == rpr.iterations),
        "AmgXPCG: a negligible absolute tolerance changes the resolution");
  // single precision: the multigrid in single precision preconditions the
  // conjugate gradient of MFEM as well as in double precision; the
  // conjugate gradient of AmgX in single precision can't reduce the
  // residual below the precision of the single precision
  auto cgm = makeConjugateGradient<is_parallel>(
      fespace, Parameters{{"Name", "AmgX"},
                          {"Options", Parameters{{"Precision", "Single"}}}});
  const auto rm = solve<is_parallel>(fespace, *(cgm.linear_solver), A, B);
  print("CG, AmgX, single precision", rm);
  check(rm.converged && (rm.residual < 1e-8) &&
            close(rm.iterations, r1.iterations),
        "CG, AmgX, single precision: resolution");
  // one cycle in single and in double precision: different, close results
  {
    auto c64 = generate<is_parallel>(fespace, "AmgX", Parameters{});
    auto c32 = generate<is_parallel>(fespace, "AmgX",
                                     Parameters{{"Precision", "Single"}});
    const auto r64 = solve<is_parallel>(fespace, *(c64.linear_solver), A, B);
    const auto r32 = solve<is_parallel>(fespace, *(c32.linear_solver), A, B);
    auto d = mfem::Vector(r64.x.Size());
    d.UseDevice(true);
    subtract(r64.x, r32.x, d);
    const auto e = std::sqrt(dot<is_parallel>(fespace, d, d) /
                             dot<is_parallel>(fespace, r64.x, r64.x));
    if (rank == 0) {
      std::cout << "AmgX, one cycle in single precision: relative "
                << "difference to double precision " << e << '\n';
    }
    check((e > 1e-12) && (e < 1e-4),
          "AmgX, single precision: the computations are not in single "
          "precision");
  }
  auto pcgm =
      generate<is_parallel>(fespace, "AmgXPCG",
                            Parameters{{"RelativeTolerance", 1e-5},
                                       {"MaximumNumberOfIterations", 200},
                                       {"Precision", "Single"}});
  const auto rpm = solve<is_parallel>(fespace, *(pcgm.linear_solver), A, B);
  print("AmgXPCG, single precision", rpm);
  check(rpm.converged && (rpm.residual < 1e-4),
        "AmgXPCG, single precision: resolution");
  // several solvers alive at the same time, the first one destroyed before
  // the last one is used
  {
    auto a1 = generate<is_parallel>(fespace, "AmgXPCG",
                                    Parameters{{"RelativeTolerance", 1e-10}});
    auto a2 = generate<is_parallel>(fespace, "AmgXPCG",
                                    Parameters{{"RelativeTolerance", 1e-10}});
    const auto ra1 = solve<is_parallel>(fespace, *(a1.linear_solver), A, B);
    a1 = LinearSolverHandler{};
    const auto ra2 = solve<is_parallel>(fespace, *(a2.linear_solver), A, B);
    auto a3 = generate<is_parallel>(fespace, "AmgXPCG",
                                    Parameters{{"RelativeTolerance", 1e-10}});
    const auto ra3 = solve<is_parallel>(fespace, *(a3.linear_solver), A, B);
    for (const auto& r : {ra1, ra2, ra3}) {
      check(r.converged && (r.residual < 1e-9),
            "AmgXPCG: resolution with several solvers");
    }
  }
  // failed setup: the solver can't be used with the previous operator;
  // vectors of invalid sizes
  {
    auto a = generate<is_parallel>(fespace, "AmgX", Parameters{});
    auto& ls = *(a.linear_solver);
    ls.SetOperator(A);
    auto x = mfem::Vector(B.Size());
    x.UseDevice(true);
    auto y = mfem::Vector(B.Size() + 1);
    y.UseDevice(true);
    auto raises = [](const std::function<void()>& f) {
      try {
        f();
      } catch (std::exception&) {
        return true;
      }
      return false;
    };
    check(raises([&ls, &B, &y] { ls.Mult(B, y); }),
          "AmgX: solution of invalid size accepted");
    check(raises([&ls, &x, &y] { ls.Mult(y, x); }),
          "AmgX: right hand side of invalid size accepted");
    const auto identity = mfem::IdentityOperator(B.Size());
    check(raises([&ls, &identity] { ls.SetOperator(identity); }),
          "AmgX: unsupported operator accepted");
    check(raises([&ls, &B, &x] { ls.Mult(B, x); }),
          "AmgX: resolution after a failed setup");
    ls.SetOperator(A);
    ls.Mult(B, x);
    check(hasConverged(ls), "AmgX: resolution after a new setup");
  }
  // failed resolution of an iterative solver: not converged, even after a
  // successful resolution
  {
    auto a = generate<is_parallel>(fespace, "AmgXPCG",
                                   Parameters{{"RelativeTolerance", 1e-10}});
    const auto r = solve<is_parallel>(fespace, *(a.linear_solver), A, B);
    auto y = mfem::Vector(B.Size() + 1);
    y.UseDevice(true);
    auto failed = false;
    try {
      a.linear_solver->Mult(y, y);
    } catch (std::exception&) {
      failed = true;
    }
    check(r.converged && failed && (!hasConverged(*(a.linear_solver))) &&
              (!getNumberOfIterationsAtConvergence(*(a.linear_solver))
                    .has_value()),
          "AmgXPCG: convergence reported after a failed resolution");
  }
  // scale of the right hand side: the floor of 1e-20 of the convergence
  // criteria of AmgX, the range of the single precision
  {
    auto pcg10 =
        generate<is_parallel>(fespace, "AmgXPCG",
                              Parameters{{"RelativeTolerance", 1e-10},
                                         {"MaximumNumberOfIterations", 200}});
    auto pcg5 =
        generate<is_parallel>(fespace, "AmgXPCG",
                              Parameters{{"RelativeTolerance", 1e-5},
                                         {"MaximumNumberOfIterations", 200},
                                         {"Precision", "Single"}});
    auto cg32 = makeConjugateGradient<is_parallel>(
        fespace, Parameters{{"Name", "AmgX"},
                            {"Options", Parameters{{"Precision", "Single"}}}});
    // the conjugate gradient of MFEM can't solve the systems whose right
    // hand side is multiplied by 1e-158 or 1e158 (scalar products), where
    // the norm of the right hand side underflows or overflows
    for (const auto f : {1e-15, 1e-25, 1e-40, 1e45, 1e-158, 1e158}) {
      const auto extreme = (f < 1e-100) || (f > 1e100);
      auto Bf = mfem::Vector(B.Size());
      Bf.UseDevice(true);
      Bf.Set(f, B);
      const auto r10 =
          solve<is_parallel>(fespace, *(pcg10.linear_solver), A, Bf);
      const auto r5 = solve<is_parallel>(fespace, *(pcg5.linear_solver), A, Bf);
      const auto r32 =
          extreme ? rm
                  : solve<is_parallel>(fespace, *(cg32.linear_solver), A, Bf);
      // the right hand side is scaled by a power of two by the solvers
      // (the rounding of the multiplication by `f` may change the number
      // of iterations)
      auto os = std::ostringstream{};
      os << f;
      const auto factor = os.str();
      check(r10.converged && (r10.residual < 1e-9) &&
                close(r10.iterations, rp.iterations),
            "AmgXPCG: right hand side multiplied by " + factor);
      check(
          r5.converged && (r5.residual < 1e-4) &&
              close(r5.iterations, rpm.iterations),
          "AmgXPCG, single precision: right hand side multiplied by " + factor);
      check(r32.converged && (r32.residual < 1e-8) &&
                close(r32.iterations, rm.iterations),
            "CG, AmgX, single precision: right hand side multiplied by " +
                factor);
    }
    // components of different scales on the processes: the factor is the
    // same for all the processes
    {
      auto Bp = mfem::Vector(B.Size());
      Bp.UseDevice(true);
      Bp.Set(rank == 0 ? 1. : 1e-3, B);
      const auto r = solve<is_parallel>(fespace, *(pcg10.linear_solver), A, Bp);
      check(r.converged && (r.residual < 1e-9),
            "AmgXPCG: right hand side of different scales on the processes");
    }
    // zero right hand side: zero solution, without iteration
    auto Z = mfem::Vector(B.Size());
    Z.UseDevice(true);
    Z = 0.;
    const auto rz =
        solve<is_parallel>(fespace, *(pcg10.linear_solver), A, Z, 1.);
    check(rz.converged && (rz.iterations == 0) &&
              (dot<is_parallel>(fespace, rz.x, rz.x) <= 0),
          "AmgXPCG: zero right hand side");
  }
  // all the solvers of AmgX destroyed, then a new one created (AmgX is
  // finalised, then initialised again)
  if (release) {
    release();
  }
  cg = LinearSolverHandler{};
  cgs = LinearSolverHandler{};
  cgf = LinearSolverHandler{};
  amg = LinearSolverHandler{};
  pcg = LinearSolverHandler{};
  pcg2 = LinearSolverHandler{};
  pcga = LinearSolverHandler{};
  amg2 = LinearSolverHandler{};
  pcgc = LinearSolverHandler{};
  pcgr = LinearSolverHandler{};
  pcgra = LinearSolverHandler{};
  cgm = LinearSolverHandler{};
  pcgm = LinearSolverHandler{};
  {
    auto a = generate<is_parallel>(fespace, "AmgXPCG",
                                   Parameters{{"RelativeTolerance", 1e-10}});
    const auto r = solve<is_parallel>(fespace, *(a.linear_solver), A, B);
    check(r.converged && (r.residual < 1e-9),
          "AmgXPCG: resolution after the destruction of all the solvers");
  }
  // invalid options
  const auto invalid = std::vector<std::pair<std::string, Parameters>>{
      {"AmgX",
       Parameters{{"Configuration", aggregation}, {"ConfigurationFile", file}}},
      {"AmgX",
       Parameters{{"Configuration", aggregation}, {"VerbosityLevel", 1}}},
      {"AmgX", Parameters{{"Configuration", ""}}},
      {"AmgX", Parameters{{"ConfigurationFile", "missing-file.json"}}},
      {"AmgX", Parameters{{"Tolerance", 1e-8}}},
      {"AmgXPCG", Parameters{{"RelativeTolerance", -1.}}},
      {"AmgXPCG", Parameters{{"RelativeTolerance", 0.}}},
      {"AmgXPCG", Parameters{{"MaximumNumberOfIterations", 0}}},
      {"AmgXPCG", Parameters{{"Configuration", aggregation}}},
      {"AmgX", Parameters{{"Configuration", aggregation},
                          {"AmgXParameters", Parameters{{"max_iters", 2}}}}},
      {"AmgX",
       Parameters{{"AmgXParameters", Parameters{{"max_iters", Parameters{}}}}}},
      {"AmgX", Parameters{{"AmgXParameters", Parameters{{"max=iters", 2}}}}},
      {"AmgX", Parameters{{"AmgXParameters",
                           Parameters{{"interpolator", "D2, amg:cycle=W"}}}}},
      {"AmgXPCG",
       Parameters{{"AmgXParameters", Parameters{{"main:max_iters", 2}}}}},
      {"AmgX", Parameters{{"AmgXParameters",
                           Parameters{{"interpolator", "D2;amg:cycle=W"}}}}},
      {"AmgX", Parameters{{"Precision", "Quadruple"}}},
      {"AmgX", Parameters{{"Multigrid", "Smoothed"}}},
      {"AmgX",
       Parameters{{"Configuration", aggregation}, {"Multigrid", "Aggregation"}}},
      {"AmgXPCG", Parameters{{"Multigrid", "Classic"}}},
      {"AmgXPCG",
       Parameters{{"RelativeTolerance", 1e-8}, {"Precision", "Single"}}},
      {"AmgXPCG", Parameters{{"RelativeTolerance", 0.},
                             {"AbsoluteTolerance", 1e-8},
                             {"Precision", "Single"}}},
      {"AmgXPCG", Parameters{{"RelativeTolerance", 1e-13}}},
      {"AmgXPCG",
       Parameters{{"RelativeTolerance", 0.}, {"AbsoluteTolerance", 1e-21}}},
      {"AmgXPCG", Parameters{{"RelativeTolerance", 1e-5},
                             {"AbsoluteTolerance", 1e-8},
                             {"Precision", "Single"}}},
      {"AmgX", Parameters{{"Configuration",
                           "config_version=2, solver(s)=NOT_A_SOLVER"}}},
      {"CGSolver",
       Parameters{{"Preconditioner",
                   Parameters{{"Name", "AmgX"},
                              {"Options", Parameters{{"Depth", 1}}}}}}}};
  for (const auto& [name, parameters] : invalid) {
    auto ctx = Context{};
    const auto ls = LinearSolverFactory<is_parallel>::getFactory().generate(
        ctx, name, fespace, parameters);
    check(isInvalid(ls), "invalid options of '" + name + "' accepted");
  }
  if constexpr (is_parallel) {
    // AmgX is not a preconditioner of hypre
    auto ctx = Context{};
    const auto ls = LinearSolverFactory<true>::getFactory().generate(
        ctx, "HyprePCG", fespace,
        Parameters{{"Preconditioner", Parameters{{"Name", "AmgX"}}}});
    check(isInvalid(ls), "AmgX accepted as a preconditioner of HyprePCG");
  }
}  // end of checkAmgX

//! \brief Laplace equation on a sequential mesh
struct SequentialSystem {
  /*!
   * \param[in] mesh: mesh
   * \param[in] fec: finite element collection
   */
  SequentialSystem(mfem::Mesh& mesh, mfem::FiniteElementCollection& fec)
      : fespace(&mesh, &fec), a(&fespace) {
    auto ess_bdr = mfem::Array<int>(mesh.bdr_attributes.Max());
    ess_bdr = 1;
    auto ess_tdofs = mfem::Array<int>{};
    this->fespace.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
    this->a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
    this->a.SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
    this->a.Assemble();
    auto one = mfem::ConstantCoefficient(1.);
    auto b = mfem::LinearForm(&(this->fespace));
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    b.Assemble();
    auto x = mfem::GridFunction(&(this->fespace));
    x = 0.;
    // without prolongation, the right hand side given by MFEM refers to the
    // memory of the linear form: it is copied
    auto X = mfem::Vector{};
    auto rhs = mfem::Vector{};
    this->a.FormLinearSystem(ess_tdofs, x, b, this->A, X, rhs);
    this->B = rhs;
  }
  //! \brief finite element space
  mfem::FiniteElementSpace fespace;
  //! \brief bilinear form
  mfem::BilinearForm a;
  //! \brief matrix
  mfem::SparseMatrix A;
  //! \brief right hand side
  mfem::Vector B;
};

#ifdef MFEM_USE_MPI

//! \brief Laplace equation on a parallel mesh
struct ParallelSystem {
  /*!
   * \param[in] comm: communicator
   * \param[in] mesh: sequential mesh, partitioned
   * \param[in] fec: finite element collection
   */
  ParallelSystem(MPI_Comm comm,
                 mfem::Mesh& mesh,
                 mfem::FiniteElementCollection& fec)
      : pmesh(comm, mesh), fespace(&pmesh, &fec), a(&fespace) {
    auto ess_bdr = mfem::Array<int>(this->pmesh.bdr_attributes.Max());
    ess_bdr = 1;
    auto ess_tdofs = mfem::Array<int>{};
    this->fespace.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
    this->a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
    this->a.Assemble();
    auto one = mfem::ConstantCoefficient(1.);
    auto b = mfem::ParLinearForm(&(this->fespace));
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    b.Assemble();
    auto x = mfem::ParGridFunction(&(this->fespace));
    x = 0.;
    auto X = mfem::Vector{};
    this->a.FormLinearSystem(ess_tdofs, x, b, this->A, X, this->B);
  }
  //! \brief mesh
  mfem::ParMesh pmesh;
  //! \brief finite element space
  mfem::ParFiniteElementSpace fespace;
  //! \brief bilinear form
  mfem::ParBilinearForm a;
  //! \brief matrix
  mfem::HypreParMatrix A;
  //! \brief right hand side
  mfem::Vector B;
};

/*!
 * \brief checks the communicators accepted by the parallel solvers of AmgX
 * and the matrices refused
 * \param[in] mesh: sequential mesh
 * \param[in] fec: finite element collection
 * \param[in] world: system on `MPI_COMM_WORLD`
 * \param[in] sequential: sequential system of the process
 */
static void checkCommunicators(mfem::Mesh& mesh,
                               mfem::FiniteElementCollection& fec,
                               ParallelSystem& world,
                               SequentialSystem& sequential) {
  using namespace mfem_mgis;
  auto size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // communicator of one process: each process solves its own system
  {
    auto self = ParallelSystem(MPI_COMM_SELF, mesh, fec);
    auto ls = generate<true>(self.fespace, "AmgXPCG",
                             Parameters{{"RelativeTolerance", 1e-10}});
    const auto r =
        solve<true>(self.fespace, *(ls.linear_solver), self.A, self.B);
    check(r.converged && (r.residual < 1e-9),
          "AmgXPCG: resolution on a communicator of one process");
  }
  // processes in reverse order: refused
  if (size > 1) {
    auto reversed = MPI_Comm{};
    MPI_Comm_split(MPI_COMM_WORLD, 0, -rank, &reversed);
    {
      auto system = ParallelSystem(reversed, mesh, fec);
      auto ctx = Context{};
      const auto ls = LinearSolverFactory<true>::getFactory().generate(
          ctx, "AmgX", system.fespace, Parameters{});
      check(isInvalid(ls),
            "AmgX: communicator with the processes in reverse order "
            "accepted");
    }
    MPI_Comm_free(&reversed);
  }
  // a process owning no row of the matrix (the first process owns all the
  // rows): refused by all the processes
  if (size > 1) {
    const auto n = static_cast<HYPRE_BigInt>(sequential.A.Height());
    auto row_starts = std::array<HYPRE_BigInt, 2>{0, n};
    if (rank != 0) {
      row_starts = {n, n};
    }
    auto empty = mfem::SparseMatrix(0, 0);
    empty.Finalize();
    auto M = mfem::HypreParMatrix(MPI_COMM_WORLD, n, row_starts.data(),
                                  rank == 0 ? &(sequential.A) : &empty);
    auto ls = generate<true>(world.fespace, "AmgX", Parameters{});
    auto refused = false;
    try {
      ls.linear_solver->SetOperator(M);
    } catch (std::exception&) {
      refused = true;
    }
    check(refused, "AmgX: matrix with a process owning no row accepted");
  }
}  // end of checkCommunicators

#endif /* MFEM_USE_MPI */

/*!
 * \brief checks that the solvers of AmgX are refused when no GPU is visible
 * \param[in] fespace: finite element space
 */
static void checkNoGPU(mfem_mgis::FiniteElementSpace<false>& fespace) {
  using namespace mfem_mgis;
  const auto solvers = std::vector<std::pair<std::string, Parameters>>{
      {"AmgX", Parameters{}},
      {"AmgXPCG", Parameters{}},
      {"CGSolver",
       Parameters{{"Preconditioner", Parameters{{"Name", "AmgX"}}}}}};
  for (const auto& [name, parameters] : solvers) {
    auto ctx = Context{};
    const auto ls = LinearSolverFactory<false>::getFactory().generate(
        ctx, name, fespace, parameters);
    check(isInvalid(ls), "'" + name + "' accepted without GPU");
    if (isInvalid(ls)) {
      const auto msg = ctx.getErrorMessage();
      std::cout << name << ": " << msg << '\n';
      check(msg.find("no GPU available") != std::string::npos,
            "'" + name + "': GPU not checked by the factory");
    }
  }
}  // end of checkNoGPU

/*!
 * \brief runs the checks
 * \param[in] mesh: sequential mesh
 * \param[in] fec: finite element collection
 */
static void run(mfem::Mesh& mesh, mfem::FiniteElementCollection& fec) {
  using namespace mfem_mgis;
  if (no_gpu) {
    auto sequential = SequentialSystem(mesh, fec);
    checkNoGPU(sequential.fespace);
  } else if (parallel) {
#ifdef MFEM_USE_MPI
    auto sequential = SequentialSystem(mesh, fec);
    // a sequential solver of AmgX on the first process only, the first
    // solver of AmgX of this process: the creation and the destruction of
    // the resources of AmgX are not collective
    if (rank == 0) {
      auto ls0 = generate<false>(sequential.fespace, "AmgXPCG",
                                 Parameters{{"RelativeTolerance", 1e-10}});
      const auto r = solve<false>(sequential.fespace, *(ls0.linear_solver),
                                  sequential.A, sequential.B);
      check(r.converged && (r.residual < 1e-9),
            "AmgXPCG: sequential resolution on the first process only");
    }
    // the other processes wait for the first one by point to point
    // messages, which can't match a collective operation made by the
    // destruction of the resources of AmgX (a deadlock would follow)
    auto token = 0;
    if (rank == 0) {
      auto size = 0;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      for (auto p = 1; p != size; ++p) {
        MPI_Send(&token, 1, MPI_INT, p, 0, MPI_COMM_WORLD);
      }
    } else {
      MPI_Recv(&token, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // a sequential solver of AmgX, created before the parallel ones, each
    // process solving its own sequential system
    auto ls = generate<false>(sequential.fespace, "AmgXPCG",
                              Parameters{{"RelativeTolerance", 1e-10}});
    const auto rs0 = solve<false>(sequential.fespace, *(ls.linear_solver),
                                  sequential.A, sequential.B);
    check(rs0.converged && (rs0.residual < 1e-9),
          "AmgXPCG: sequential resolution in a parallel computation");
    auto world = ParallelSystem(MPI_COMM_WORLD, mesh, fec);
    // the sequential solver is used again, then destroyed before the
    // destruction of all the solvers
    checkAmgX<true>(world.fespace, world.A, world.B, [&sequential, &ls] {
      const auto rs1 = solve<false>(sequential.fespace, *(ls.linear_solver),
                                    sequential.A, sequential.B);
      check(rs1.converged && (rs1.residual < 1e-9),
            "AmgXPCG: sequential resolution in a parallel computation, "
            "after the parallel ones");
      ls = LinearSolverHandler{};
    });
    checkCommunicators(mesh, fec, world, sequential);
#else  /* MFEM_USE_MPI */
    throw std::runtime_error("parallel computations are not supported");
#endif /* MFEM_USE_MPI */
  } else {
    auto sequential = SequentialSystem(mesh, fec);
    checkAmgX<false>(sequential.fespace, sequential.A, sequential.B);
  }
}  // end of run

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "parallel computations");
  args.AddOption(&no_gpu, "-no-gpu", "--no-gpu", "-gpu", "--gpu",
                 "no GPU expected: the solvers of AmgX must be refused");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
#ifdef MFEM_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif /* MFEM_USE_MPI */
  auto mesh =
      mfem::Mesh::MakeCartesian3D(32, 32, 32, mfem::Element::HEXAHEDRON);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    run(mesh, fec);
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (rank == 0) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
