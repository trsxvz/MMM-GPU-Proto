/*!
 * \file   OperatorSmoothersTest.cxx
 * \brief  Test of the smoothers of the factory of the linear solvers which
 * only use products by the operator and its diagonal
 * (`OperatorJacobiSmoother`, `OperatorChebyshevSmoother`), as
 * preconditioners of the conjugate gradient of MFEM: Laplace equation on a
 * square, the essential degrees of freedom being eliminated with a diagonal
 * set to one. Each smoother gives the solution, with the same number of
 * iterations for the assembled operator and for an operator applied without
 * assembly; the Chebyshev smoother reduces the number of iterations, the
 * more so as its order is higher. With the `--parallel` option, the
 * computations are parallel.
 * \date   10/10/2026
 */

#include <cmath>
#include <memory>
#include <string>
#include <cstdlib>
#include <iostream>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"

//! \brief number of failed checks
static int failures = 0;
//! \brief if true, the computations are parallel
static bool parallel = false;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \brief operator applied by an assembled operator, but which is neither a
 * sparse matrix nor a parallel matrix of hypre: the smoothers compute its
 * diagonal by `AssembleDiagonal`, as for an operator applied without
 * assembly
 */
struct WrappedOperator final : mfem::Operator {
  //! \param[in] a: assembled operator
  explicit WrappedOperator(const mfem::Operator& a)
      : mfem::Operator(a.Height(), a.Width()), op(a) {}
  void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
    this->op.Mult(x, y);
  }
  void AssembleDiagonal(mfem::Vector& d) const override {
    if (const auto* const m = dynamic_cast<const mfem::SparseMatrix*>(&op);
        m != nullptr) {
      m->GetDiag(d);
    } else {
      this->op.AssembleDiagonal(d);
    }
  }

 private:
  //! \brief assembled operator
  const mfem::Operator& op;
};

/*!
 * \brief operator `S A S` applied without assembly, `S` being a diagonal
 * scaling: the Jacobi and Chebyshev smoothers, built on the scaled
 * diagonal, are invariant by this scaling (same numbers of iterations of the
 * conjugate gradient), which fails if the diagonal is not taken into
 * account
 */
struct ScaledOperator final : mfem::Operator {
  /*!
   * \param[in] a: operator
   * \param[in] s: scaling
   */
  ScaledOperator(const mfem::Operator& a, const mfem::Vector& s)
      : mfem::Operator(a.Height(), a.Width()), op(a), scaling(s) {}
  // the values are scaled on the host, where they are made valid (the
  // vectors may be used on a device)
  void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
    auto sx = mfem::Vector(x.Size());
    const auto* const px = x.HostRead();
    const auto* const ps = this->scaling.HostRead();
    auto* const psx = sx.HostWrite();
    for (int i = 0; i != x.Size(); ++i) {
      psx[i] = ps[i] * px[i];
    }
    this->op.Mult(sx, y);
    auto* const py = y.HostReadWrite();
    for (int i = 0; i != y.Size(); ++i) {
      py[i] *= ps[i];
    }
  }
  void AssembleDiagonal(mfem::Vector& d) const override {
    if (const auto* const m = dynamic_cast<const mfem::SparseMatrix*>(&op);
        m != nullptr) {
      m->GetDiag(d);
    } else {
      this->op.AssembleDiagonal(d);
    }
    const auto* const ps = this->scaling.HostRead();
    auto* const pd = d.HostReadWrite();
    for (int i = 0; i != d.Size(); ++i) {
      pd[i] *= ps[i] * ps[i];
    }
  }

 private:
  //! \brief operator
  const mfem::Operator& op;
  //! \brief scaling
  const mfem::Vector& scaling;
};

//! \brief result of a resolution
struct Resolution {
  //! \brief if true, the resolution succeeded
  bool converged = false;
  //! \brief number of iterations
  int iterations = 0;
  //! \brief norm of the residual relative to the norm of the right hand side
  double residual = 1;
};

/*!
 * \return the resolution of a linear system by a linear solver
 * \param[in] fespace: finite element space
 * \param[in, out] ls: linear solver
 * \param[in] A: operator
 * \param[in] B: right hand side
 */
template <bool is_parallel>
static Resolution solve(mfem_mgis::FiniteElementSpace<is_parallel>& fespace,
                        mfem::Solver& ls,
                        const mfem::Operator& A,
                        const mfem::Vector& B) {
  using namespace mfem_mgis;
  auto X = mfem::Vector(B.Size());
  X = 0.;
  ls.SetOperator(A);
  ls.Mult(B, X);
  auto r = mfem::Vector(B.Size());
  A.Mult(X, r);
  r -= B;
  auto dot = [&fespace](const mfem::Vector& a, const mfem::Vector& b) {
    if constexpr (is_parallel) {
      return mfem::InnerProduct(fespace.GetComm(), a, b);
    } else {
      static_cast<void>(fespace);
      return mfem::InnerProduct(a, b);
    }
  };
  return {.converged = hasConverged(ls),
          .iterations = getNumberOfIterationsAtConvergence(ls).value_or(0),
          .residual = std::sqrt(dot(r, r) / dot(B, B))};
}  // end of solve

/*!
 * \return the conjugate gradient of the factory preconditioned by a
 * smoother
 * \param[in] fespace: finite element space
 * \param[in] smoother: preconditioner
 */
template <bool is_parallel>
static mfem_mgis::LinearSolverHandler makeConjugateGradient(
    mfem_mgis::FiniteElementSpace<is_parallel>& fespace,
    const mfem_mgis::Parameters& smoother) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto ls = LinearSolverFactory<is_parallel>::getFactory().generate(
      ctx, "CGSolver", fespace,
      Parameters{{"RelativeTolerance", 1e-10},
                 {"AbsoluteTolerance", 0.},
                 {"MaximumNumberOfIterations", 5000},
                 {"Preconditioner", smoother}});
  if (isInvalid(ls)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    std::exit(EXIT_FAILURE);
  }
  return ls;
}  // end of makeConjugateGradient

/*!
 * \brief checks the smoothers on a linear system
 * \param[in] fespace: finite element space
 * \param[in] A: assembled operator, the essential degrees of freedom being
 * eliminated with a diagonal set to one
 * \param[in] B: right hand side
 */
template <bool is_parallel>
static void checkSmoothers(mfem_mgis::FiniteElementSpace<is_parallel>& fespace,
                           const mfem::Operator& A,
                           const mfem::Vector& B) {
  using namespace mfem_mgis;
  const auto wrapped = WrappedOperator{A};
  // diagonal scaling from 1 to 10, and the scaled right hand side
  auto scaling = mfem::Vector(B.Size());
  auto SB = mfem::Vector(B.Size());
  const auto* const pB = B.HostRead();
  auto* const ps = scaling.HostWrite();
  auto* const pSB = SB.HostWrite();
  for (int i = 0; i != B.Size(); ++i) {
    ps[i] = 1 + 9 * ((7 * i) % 11) / 10.;
    pSB[i] = ps[i] * pB[i];
  }
  const auto scaled = ScaledOperator{wrapped, scaling};
  auto iterations = std::vector<int>{};
  for (const auto& [name, smoother] :
       {std::pair{std::string{"Jacobi"},
                  Parameters{{"Name", "OperatorJacobiSmoother"}}},
        std::pair{std::string{"Chebyshev of order 2"},
                  Parameters{{"Name", "OperatorChebyshevSmoother"},
                             {"Options", Parameters{{"Order", 2}}}}},
        std::pair{std::string{"Chebyshev of order 4"},
                  Parameters{{"Name", "OperatorChebyshevSmoother"},
                             {"Options", Parameters{{"Order", 4}}}}}}) {
    auto ls = makeConjugateGradient<is_parallel>(fespace, smoother);
    const auto ra = solve<is_parallel>(fespace, *(ls.linear_solver), A, B);
    // the smoother is set up again for each operator
    const auto rw =
        solve<is_parallel>(fespace, *(ls.linear_solver), wrapped, B);
    const auto rs =
        solve<is_parallel>(fespace, *(ls.linear_solver), scaled, SB);
    std::cout << name << ": " << ra.iterations << " iterations (assembled), "
              << rw.iterations << " (applied without assembly), "
              << rs.iterations << " (scaled), residuals " << ra.residual << ", "
              << rw.residual << ", " << rs.residual << '\n';
    for (const auto& [r, kind] :
         {std::pair{ra, "assembled operator"},
          std::pair{rw, "operator applied without assembly"},
          std::pair{rs, "scaled operator"}}) {
      check(r.converged && (r.residual < 1e-9),
            name + ", " + kind + ": resolution");
    }
    check(std::abs(ra.iterations - rw.iterations) <= 1,
          name + ": different numbers of iterations");
    check(std::abs(rw.iterations - rs.iterations) <= 1,
          name + ": numbers of iterations changed by a diagonal scaling");
    iterations.push_back(ra.iterations);
  }
  check((iterations.size() == 3) && (iterations[1] < iterations[0]) &&
            (iterations[2] < iterations[1]),
        "the Chebyshev smoother does not reduce the number of iterations");
  // invalid options
  for (const auto order : {0, 6}) {
    auto ctx = Context{};
    const auto ls = LinearSolverFactory<is_parallel>::getFactory().generate(
        ctx, "CGSolver", fespace,
        Parameters{{"Preconditioner",
                    Parameters{{"Name", "OperatorChebyshevSmoother"},
                               {"Options", Parameters{{"Order", order}}}}}});
    check(isInvalid(ls), "an order " + std::to_string(order) +
                             " of the Chebyshev smoother is accepted");
  }
}  // end of checkSmoothers

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "parallel computations");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  auto mesh = mfem::Mesh::MakeCartesian2D(32, 32, mfem::Element::QUADRILATERAL);
  auto fec = mfem::H1_FECollection(1, 2);
  auto ess_bdr = mfem::Array<int>(mesh.bdr_attributes.Max());
  ess_bdr = 1;
  auto one = mfem::ConstantCoefficient(1.);
  if (parallel) {
#ifdef MFEM_USE_MPI
    auto pmesh = mfem::ParMesh(MPI_COMM_WORLD, mesh);
    auto fespace = mfem::ParFiniteElementSpace(&pmesh, &fec);
    auto ess_tdofs = mfem::Array<int>{};
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
    auto a = mfem::ParBilinearForm(&fespace);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
    a.Assemble();
    auto b = mfem::ParLinearForm(&fespace);
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    b.Assemble();
    auto x = mfem::ParGridFunction(&fespace);
    x = 0.;
    auto A = mfem::HypreParMatrix{};
    auto B = mfem::Vector{};
    auto X = mfem::Vector{};
    a.FormLinearSystem(ess_tdofs, x, b, A, X, B);
    checkSmoothers<true>(fespace, A, B);
#else  /* MFEM_USE_MPI */
    std::cerr << "parallel computations are not supported\n";
    return EXIT_FAILURE;
#endif /* MFEM_USE_MPI */
  } else {
    auto fespace = mfem::FiniteElementSpace(&mesh, &fec);
    auto ess_tdofs = mfem::Array<int>{};
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
    auto a = mfem::BilinearForm(&fespace);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
    a.SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
    a.Assemble();
    auto b = mfem::LinearForm(&fespace);
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    b.Assemble();
    auto x = mfem::GridFunction(&fespace);
    x = 0.;
    auto A = mfem::SparseMatrix{};
    auto B = mfem::Vector{};
    auto X = mfem::Vector{};
    a.FormLinearSystem(ess_tdofs, x, b, A, X, B);
    checkSmoothers<false>(fespace, A, B);
  }
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
