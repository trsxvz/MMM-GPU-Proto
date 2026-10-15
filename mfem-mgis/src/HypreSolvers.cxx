/*!
 * \file   src/HypreSolvers.cxx
 * \brief  This file implements the solvers and the preconditioners of hypre
 * which are not wrapped by MFEM (see `HypreSolvers.hxx`)
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#ifdef MFEM_USE_MPI

#include <utility>
#include "_hypre_krylov.h"
#include "MGIS/Raise.hxx"
#include "HypreSolvers.hxx"

namespace mfem_mgis::internals {

  /*!
   * \return the matrix of hypre given to a solver of hypre, raises an
   * exception for another operator
   * \param[in] op: operator
   * \param[in] n: name of the solver
   */
  [[nodiscard]] static const mfem::HypreParMatrix& getHypreParMatrix(
      const mfem::Operator& op, const char* const n) {
    const auto* const a = dynamic_cast<const mfem::HypreParMatrix*>(&op);
    if (a == nullptr) {
      mgis::raise(std::string{n} +
                  "::SetOperator: the operator is not a matrix of hypre");
    }
    return *a;
  }  // end of getHypreParMatrix

  HypreKrylovSolverFunctions getHypreBiCGSTABFunctions() {
    return {.create = HYPRE_ParCSRBiCGSTABCreate,
            .destroy = HYPRE_ParCSRBiCGSTABDestroy,
            .setup = HYPRE_ParCSRBiCGSTABSetup,
            .solve = HYPRE_ParCSRBiCGSTABSolve,
            .set_tolerance = HYPRE_ParCSRBiCGSTABSetTol,
            .set_absolute_tolerance = HYPRE_ParCSRBiCGSTABSetAbsoluteTol,
            .set_maximum_number_of_iterations = HYPRE_ParCSRBiCGSTABSetMaxIter,
            .set_print_level = HYPRE_ParCSRBiCGSTABSetPrintLevel,
            .set_preconditioner = HYPRE_ParCSRBiCGSTABSetPrecond,
            .get_number_of_iterations = HYPRE_ParCSRBiCGSTABGetNumIterations,
            .get_final_relative_residual_norm =
                HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm,
            // only given by the internal interface of hypre
            .get_converged = [](HYPRE_Solver s, HYPRE_Int* c) {
              return hypre_BiCGSTABGetConverged(static_cast<void*>(s), c);
            }};
  }  // end of getHypreBiCGSTABFunctions

  HypreKrylovSolverFunctions getHypreCOGMRESFunctions() {
    return {.create = HYPRE_ParCSRCOGMRESCreate,
            .destroy = HYPRE_ParCSRCOGMRESDestroy,
            .setup = HYPRE_ParCSRCOGMRESSetup,
            .solve = HYPRE_ParCSRCOGMRESSolve,
            .set_tolerance = HYPRE_ParCSRCOGMRESSetTol,
            .set_absolute_tolerance = HYPRE_ParCSRCOGMRESSetAbsoluteTol,
            .set_maximum_number_of_iterations = HYPRE_ParCSRCOGMRESSetMaxIter,
            .set_print_level = HYPRE_ParCSRCOGMRESSetPrintLevel,
            .set_preconditioner = HYPRE_ParCSRCOGMRESSetPrecond,
            .get_number_of_iterations = HYPRE_ParCSRCOGMRESGetNumIterations,
            .get_final_relative_residual_norm =
                HYPRE_ParCSRCOGMRESGetFinalRelativeResidualNorm,
            .get_converged = HYPRE_COGMRESGetConverged,
            .set_krylov_dimension = HYPRE_ParCSRCOGMRESSetKDim};
  }  // end of getHypreCOGMRESFunctions

  HypreKrylovSolverFunctions getHypreLGMRESFunctions() {
    return {.create = HYPRE_ParCSRLGMRESCreate,
            .destroy = HYPRE_ParCSRLGMRESDestroy,
            .setup = HYPRE_ParCSRLGMRESSetup,
            .solve = HYPRE_ParCSRLGMRESSolve,
            .set_tolerance = HYPRE_ParCSRLGMRESSetTol,
            .set_absolute_tolerance = HYPRE_ParCSRLGMRESSetAbsoluteTol,
            .set_maximum_number_of_iterations = HYPRE_ParCSRLGMRESSetMaxIter,
            .set_print_level = HYPRE_ParCSRLGMRESSetPrintLevel,
            .set_preconditioner = HYPRE_ParCSRLGMRESSetPrecond,
            .get_number_of_iterations = HYPRE_ParCSRLGMRESGetNumIterations,
            .get_final_relative_residual_norm =
                HYPRE_ParCSRLGMRESGetFinalRelativeResidualNorm,
            .get_converged = HYPRE_LGMRESGetConverged,
            .set_krylov_dimension = HYPRE_ParCSRLGMRESSetKDim};
  }  // end of getHypreLGMRESFunctions

  HypreKrylovSolver::HypreKrylovSolver(MPI_Comm comm,
                                       const HypreKrylovSolverFunctions& f)
      : functions(f) {
    this->functions.create(comm, &(this->solver));
    // the default values of the Krylov solvers of hypre wrapped by MFEM (the
    // ones of `mfem::HypreGMRES` for the dimension of the Krylov space)
    this->SetTol(this->relative_tolerance);
    this->SetMaxIter(this->maximum_number_of_iterations);
    std::ignore = this->SetKDim(50);
    // as the Krylov solvers of hypre wrapped by MFEM, a resolution which
    // does not converge is not an error (see `hasConverged`)
    this->SetErrorMode(IGNORE_HYPRE_ERRORS);
  }  // end of HypreKrylovSolver

  HypreKrylovSolver::operator HYPRE_Solver() const { return this->solver; }

  HYPRE_PtrToParSolverFcn HypreKrylovSolver::SetupFcn() const {
    return this->functions.setup;
  }  // end of SetupFcn

  HYPRE_PtrToParSolverFcn HypreKrylovSolver::SolveFcn() const {
    return this->functions.solve;
  }  // end of SolveFcn

  void HypreKrylovSolver::SetOperator(const mfem::Operator& op) {
    // as mfem::HypreGMRES::SetOperator
    const auto& a = getHypreParMatrix(op, "HypreKrylovSolver");
    this->height = a.Height();
    this->width = a.Width();
    this->A = &a;
    if (this->preconditioner != nullptr) {
      this->preconditioner->SetOperator(a);
      this->SetPreconditioner(*(this->preconditioner));
    }
    this->setup_called = 0;
    delete this->X;
    delete this->B;
    this->B = nullptr;
    this->X = nullptr;
    this->auxB.Delete();
    this->auxB.Reset();
    this->auxX.Delete();
    this->auxX.Reset();
  }  // end of SetOperator

  void HypreKrylovSolver::SetPreconditioner(mfem::HypreSolver& p) {
    this->preconditioner = &p;
    this->functions.set_preconditioner(this->solver, p.SolveFcn(), p.SetupFcn(),
                                       p);
  }  // end of SetPreconditioner

  void HypreKrylovSolver::SetTol(const real tol) {
    this->relative_tolerance = tol;
    this->functions.set_tolerance(this->solver, tol);
  }  // end of SetTol

  void HypreKrylovSolver::SetAbsTol(const real tol) {
    this->functions.set_absolute_tolerance(this->solver, tol);
  }  // end of SetAbsTol

  void HypreKrylovSolver::SetMaxIter(const int n) {
    this->maximum_number_of_iterations = n;
    this->functions.set_maximum_number_of_iterations(this->solver, n);
  }  // end of SetMaxIter

  void HypreKrylovSolver::SetPrintLevel(const int l) {
    this->functions.set_print_level(this->solver, l);
  }  // end of SetPrintLevel

  bool HypreKrylovSolver::SetKDim(const int n) {
    if (this->functions.set_krylov_dimension == nullptr) {
      return false;
    }
    this->functions.set_krylov_dimension(this->solver, n);
    return true;
  }  // end of SetKDim

  void HypreKrylovSolver::Mult(const mfem::HypreParVector& b,
                               mfem::HypreParVector& x) const {
    if (this->A == nullptr) {
      mgis::raise("HypreKrylovSolver::Mult: the operator is not set");
    }
    if (!this->iterative_mode) {
      x.HypreWrite();
      hypre_ParVectorSetConstantValues(x, 0.0);
    }
    b.HypreRead();
    x.HypreReadWrite();
    this->Setup(b, x);
    this->error_code = this->SolveFcn()(*this, *(this->A), b, x);
    hypre_error_flag = 0;
  }  // end of Mult

  bool HypreKrylovSolver::hasConverged() const noexcept {
    // with a null tolerance (a fixed number of iterations, see
    // `boundRelativeTolerance`), hypre reports no error at the maximum number
    // of iterations: the resolution is then converged
    // the flag of hypre is reset at the beginning of each resolution and
    // only set by the convergence test of the iterations: the solvers of
    // hypre return before it when the initial residual is null, or when it
    // satisfies the absolute tolerance, without error. A breakdown, or the
    // maximum number of iterations, is reported by an error.
    auto converged = HYPRE_Int{};
    this->functions.get_converged(this->solver, &converged);
    return (converged != 0) || (this->error_code == 0);
  }  // end of hasConverged

  std::optional<int> HypreKrylovSolver::getNumberOfIterations() const noexcept {
    auto n = HYPRE_Int{};
    this->functions.get_number_of_iterations(this->solver, &n);
    return static_cast<int>(n);
  }  // end of getNumberOfIterations

  std::optional<real> HypreKrylovSolver::getRelativeTolerance() const noexcept {
    return this->relative_tolerance;
  }  // end of getRelativeTolerance

  bool HypreKrylovSolver::setRelativeTolerance(const real tol) noexcept {
    this->SetTol(tol);
    return true;
  }  // end of setRelativeTolerance

  HypreKrylovSolver::~HypreKrylovSolver() {
    this->functions.destroy(this->solver);
  }  // end of ~HypreKrylovSolver

  HypreHybridSolver::HypreHybridSolver(std::vector<HypreSetting> s,
                                       const real tol,
                                       const int n)
      : settings(std::move(s)),
        relative_tolerance(tol),
        maximum_number_of_iterations(n) {
    this->create();
    // a resolution which does not converge is not an error (see
    // `hasConverged`)
    this->SetErrorMode(IGNORE_HYPRE_ERRORS);
  }  // end of HypreHybridSolver

  void HypreHybridSolver::create() {
    if (this->solver != nullptr) {
      HYPRE_ParCSRHybridDestroy(this->solver);
    }
    HYPRE_ParCSRHybridCreate(&(this->solver));
    // the final relative residual norm is only computed with logging
    HYPRE_ParCSRHybridSetLogging(this->solver, 1);
    // the dimension of the Krylov space of `mfem::HypreGMRES` (hypre: 5)
    HYPRE_ParCSRHybridSetKDim(this->solver, 50);
    HYPRE_ParCSRHybridSetTol(this->solver, this->relative_tolerance);
    HYPRE_ParCSRHybridSetDSCGMaxIter(this->solver,
                                     this->maximum_number_of_iterations);
    HYPRE_ParCSRHybridSetPCGMaxIter(this->solver,
                                    this->maximum_number_of_iterations);
    for (const auto& setting : this->settings) {
      setting(this->solver);
    }
  }  // end of create

  HypreHybridSolver::operator HYPRE_Solver() const { return this->solver; }

  HYPRE_PtrToParSolverFcn HypreHybridSolver::SetupFcn() const {
    return HYPRE_ParCSRHybridSetup;
  }  // end of SetupFcn

  HYPRE_PtrToParSolverFcn HypreHybridSolver::SolveFcn() const {
    return HYPRE_ParCSRHybridSolve;
  }  // end of SolveFcn

  void HypreHybridSolver::SetOperator(const mfem::Operator& op) {
    const auto& a = getHypreParMatrix(op, "HypreHybridSolver");
    this->height = a.Height();
    this->width = a.Width();
    this->A = &a;
    // the inner Krylov solver of the hybrid solver of hypre is created at
    // the first resolution, with the current tolerance: the solver of hypre
    // is created again for each operator and each change of the tolerance
    this->create();
    this->setup_called = 0;
    delete this->X;
    delete this->B;
    this->B = nullptr;
    this->X = nullptr;
    this->auxB.Delete();
    this->auxB.Reset();
    this->auxX.Delete();
    this->auxX.Reset();
  }  // end of SetOperator

  bool HypreHybridSolver::hasConverged() const noexcept {
    auto r = HYPRE_Real{};
    HYPRE_ParCSRHybridGetFinalRelativeResidualNorm(this->solver, &r);
    return r <= this->relative_tolerance;
  }  // end of hasConverged

  std::optional<int> HypreHybridSolver::getNumberOfIterations() const noexcept {
    auto n = HYPRE_Int{};
    HYPRE_ParCSRHybridGetNumIterations(this->solver, &n);
    return static_cast<int>(n);
  }  // end of getNumberOfIterations

  std::optional<real> HypreHybridSolver::getRelativeTolerance() const noexcept {
    return this->relative_tolerance;
  }  // end of getRelativeTolerance

  bool HypreHybridSolver::setRelativeTolerance(const real tol) noexcept {
    this->relative_tolerance = tol;
    this->create();
    this->setup_called = 0;
    return true;
  }  // end of setRelativeTolerance

  HypreHybridSolver::~HypreHybridSolver() {
    HYPRE_ParCSRHybridDestroy(this->solver);
  }  // end of ~HypreHybridSolver

  HypreFSAI::HypreFSAI(std::vector<HypreSetting> s) : settings(std::move(s)) {
    this->create();
  }  // end of HypreFSAI

  void HypreFSAI::create() {
    if (this->solver != nullptr) {
      HYPRE_FSAIDestroy(this->solver);
    }
    HYPRE_FSAICreate(&(this->solver));
    // one application from a null initial guess, as a preconditioner
    HYPRE_FSAISetMaxIterations(this->solver, 1);
    HYPRE_FSAISetTolerance(this->solver, 0);
    HYPRE_FSAISetZeroGuess(this->solver, 1);
    for (const auto& setting : this->settings) {
      setting(this->solver);
    }
  }  // end of create

  HypreFSAI::operator HYPRE_Solver() const { return this->solver; }

  HYPRE_PtrToParSolverFcn HypreFSAI::SetupFcn() const {
    return HYPRE_FSAISetup;
  }  // end of SetupFcn

  HYPRE_PtrToParSolverFcn HypreFSAI::SolveFcn() const {
    return HYPRE_FSAISolve;
  }  // end of SolveFcn

  void HypreFSAI::SetOperator(const mfem::Operator& op) {
    const auto& a = getHypreParMatrix(op, "HypreFSAI");
    this->height = a.Height();
    this->width = a.Width();
    this->A = &a;
    this->create();
    this->setup_called = 0;
    delete this->X;
    delete this->B;
    this->B = nullptr;
    this->X = nullptr;
    this->auxB.Delete();
    this->auxB.Reset();
    this->auxX.Delete();
    this->auxX.Reset();
  }  // end of SetOperator

  HypreFSAI::~HypreFSAI() { HYPRE_FSAIDestroy(this->solver); }

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_USE_MPI */
