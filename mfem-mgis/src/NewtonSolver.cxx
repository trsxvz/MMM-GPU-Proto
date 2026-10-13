/*!
 * \file   src/NewtonSolver.cxx
 * \brief  This file implements the `NewtonSolver` class
 * \author Thomas Helfer
 * \date   29/03/2021
 */

#include <array>
#include <string>
#include <vector>
#include <iomanip>
#include <optional>
#include <tuple>
#include <utility>
#include "MFEMMGIS/Parameters.hxx"
#include "MGIS/Raise.hxx"
#include "MGIS/ErrorBacktrace.hxx"
#include "MGIS/Profiling.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/NonLinearSolvers/NewtonSolver.hxx"

namespace mfem_mgis {

  template <bool parallel>
  static void checkSolverOperator(
      const NonLinearEvolutionProblemImplementation<parallel> &p) {
    MFEM_ASSERT(p.Height() == p.Width(),
                "checkSolverOperator: "
                "a square operator is required.");
    static_cast<void>(p);
  }  // end of checkSolverOperator

#ifdef MFEM_USE_MPI

  NewtonSolver::NewtonSolver(NonLinearEvolutionProblemImplementation<true> &p)
      : NonLinearSolverBase(p) {}  // end of NewtonSolver

#endif /* MFEM_USE_MPI */

  NewtonSolver::NewtonSolver(NonLinearEvolutionProblemImplementation<false> &p)
      : NonLinearSolverBase(p) {}  // end of NewtonSolver

  void NewtonSolver::Mult(const mfem::Vector &, mfem::Vector &x) const {
    auto profiler_mult =
        this->ctx_ptr != nullptr
            ? this->ctx_ptr->startNewProfiling(
                  "NS::Mult", this->ctx_ptr->isProfilingEnabled())
            : mgis::ProfilingSection{};  // fallback RAII
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    MFEM_ASSERT(this->prec != nullptr,
                "the Solver is not set (use setLinearSolver).");
    //
    this->iterations_information.clear();
    // relative tolerance of the last linear solve, if adapted, and norm of
    // its residual, if computed (see `AdaptiveLinearSolverTolerance`)
    auto linear_solver_tolerance = std::optional<real>{};
    auto linear_residual_norm = std::optional<real>{};
    auto add_to_iterations_information = [this, &linear_solver_tolerance,
                                          &linear_residual_norm](
                                             Parameters& iteration,
                                             const bool ls_cv) {
      auto linear_solver_iterations =
          getNumberOfIterationsAtConvergence(*(this->prec));
      if (isValid(linear_solver_iterations)) {
        auto ls = Parameters{};
        ls.replaceOrInsert("NumberOfIterations",
                           linear_solver_iterations.value());
        ls.replaceOrInsert("HasConverged", ls_cv);
        if (linear_solver_tolerance.has_value()) {
          ls.replaceOrInsert("RelativeTolerance", *linear_solver_tolerance);
        }
        if (linear_residual_norm.has_value()) {
          ls.replaceOrInsert("ResidualNorm", *linear_residual_norm);
        }
        iteration.replaceOrInsert("LinearSolver", ls);
      }
      this->iterations_information.push_back(iteration);
    };
    auto &log = this->getLogStream();
    const auto shall_print = this->shallPrint();

    mfem::Vector r;  // residual vector
    mfem::Vector c;  // opposite of the Newton's correction
    r.SetSize(this->oper->Width());
    c.SetSize(this->oper->Width());

    auto updateResidual = [this, &r, &x] {
      if (this->ctx_ptr != nullptr) {
        CatchTimeSection(*(this->ctx_ptr), "NS::computeResidual");
        this->computeResidual(r, x);
      } else {
        this->computeResidual(r, x);
      }
      return this->Norm(r);
    };

    this->final_iter = size_type{};
    this->final_norm = std::numeric_limits<real>::max();

    // relative tolerance of the linear solver adapted to the decrease of the
    // residual, restored at the end of the resolution
    auto adaptive_tolerance =
        std::optional<AdaptiveLinearSolverToleranceUpdater>{};
    if (this->adaptive_linear_solver_tolerance.isEnabled()) {
      const auto created =
          this->callAndHandleFailures([this, &adaptive_tolerance] {
            adaptive_tolerance.emplace(this->adaptive_linear_solver_tolerance,
                                       *(this->prec));
            return true;
          });
      if (!created) {
        this->converged = 0;
        return;
      }
    }

    if (!this->processNewUnknownsEstimate(x)) {
      this->converged = 0;
      return;
    }
    auto norm = updateResidual();
    if (!this->reference_residual_norm.has_value()) {
      this->reference_residual_norm = norm;
    }
    // this data member is not used, but we define it by
    // consistency
    this->initial_norm = *(this->reference_residual_norm);
    auto previous_norms =
        std::array<real, 2u>{this->initial_norm, this->initial_norm};

    const auto norm_goal =
        std::max(rel_tol * (*(this->reference_residual_norm)), abs_tol);
    auto it = size_type{};

    this->converged = 0;
    while (true) {
      auto profiler_while =
          this->ctx_ptr != nullptr
              ? this->ctx_ptr->startNewProfiling(
                    "NS::Mult::WhileLoop", this->ctx_ptr->isProfilingEnabled())
              : mgis::ProfilingSection{};
      MFEM_ASSERT(mfem::IsFinite(norm), "norm = " << norm);
      if (shall_print) {
        log << "Newton iteration " << std::setw(2) << it
            << " : ||r|| = " << norm
            << ", ||r||/||r_0|| = " << norm / (*(this->reference_residual_norm))
            << '\n';
      }
      this->Monitor(it, norm, r, x);
      //
      if (norm <= norm_goal) {
        this->converged = 1;
        break;
      }
      //
      if (it >= this->max_iter) {
        break;
      }
      //
      if (adaptive_tolerance.has_value()) {
        linear_solver_tolerance =
            adaptive_tolerance->update(it, norm, norm_goal);
        if (shall_print) {
          log << "Newton iteration " << std::setw(2) << it
              << " : relative tolerance of the linear solver "
              << *linear_solver_tolerance << '\n';
        }
      }
      auto* const t =
          adaptive_tolerance.has_value() ? &(*adaptive_tolerance) : nullptr;
      const auto corrected = this->computeNewtonCorrection(c, r, x, t);
      if (t != nullptr) {
        linear_residual_norm = t->getLinearResidualNorm();
      }
      if (!corrected) {
        auto iteration = Parameters{};
        add_to_iterations_information(iteration,
                                      false);  // Reaching this line means the
                                               // computation of the
                                               // correction failed
        break;
      }
      //
      // x_{i+1} = x_i - c * [DF(x_i)]^{-1} [F(x_i)-b]
      //      add(x, -1, c, x);
      x -= c;

      if (!this->processNewUnknownsEstimate(x)) {
        ++it;
        // basic line-search in case of integration failure
        //
        // we keep the current direction, but reduce the amplitude by a factor
        // two until we find an estimate of the solution that do not lead
        // to an integration failure or that the number of iterations reaches
        // the maximum value
        while (true) {
          if (it >= this->max_iter) {
            break;
          }
          if (shall_print) {
            log << "Newton iteration " << std::setw(2) << it
                << ": reducing the amplitude of the correction by a "
                   "factor 2\n";
          }
          c *= real{1} / 2;
          x += c;
          if (this->processNewUnknownsEstimate(x)) {
            break;
          }
          ++it;
        }
        if (it >= this->max_iter) {
          break;
        }
      }

      updateResidual();
      previous_norms[0] = previous_norms[1];
      previous_norms[1] = norm;
      norm = this->Norm(r);
      //
      auto iteration = Parameters{};
      iteration.replaceOrInsert("Norm", norm);
      // Called after computeNewtonCorrection which uses prec->Mult
      add_to_iterations_information(
          iteration, true);  // here the linear solver did converge
      //
      ++it;
    }
    this->final_iter = it;
    this->final_norm = norm;
    if (this->converged == 1) {
      // estimation of the convergence order
      if (shall_print) {
        if ((it >= 2) && (norm > 100 * std::numeric_limits<real>::min()) &&
            (previous_norms[0] > 100 * std::numeric_limits<real>::min()) &&
            (previous_norms[1] > 100 * std::numeric_limits<real>::min())) {
          const auto e1 = std::log(norm / previous_norms[1]);
          const auto e2 = std::log(previous_norms[1] / previous_norms[0]);
          if (std::abs(e2) > 100 * std::numeric_limits<real>::min()) {
            log << "Convergence order " << e1 / e2 << "\n\n";
          } else {
            log << "Convergence order undefined\n\n";
          }
        } else {
          log << "Convergence order undefined\n\n";
        }
      }
    }
  }  // end of Mult

  bool NewtonSolver::setSolverParameters(Context& ctx,
                                         const Parameters& params) noexcept {
    // the options are read before being applied
    auto t = this->adaptive_linear_solver_tolerance;
    if (contains(params, "AdaptiveLinearSolverTolerance")) {
      const auto ot = AdaptiveLinearSolverTolerance::read(ctx, params);
      if (isInvalid(ot)) {
        return false;
      }
      t = *ot;
    }
    const auto others = remove(
        params, std::vector<std::string>{"AdaptiveLinearSolverTolerance"});
    if (!NonLinearSolverBase::setSolverParameters(ctx, others)) {
      return false;
    }
    this->adaptive_linear_solver_tolerance = t;
    return true;
  }  // end of setSolverParameters

  bool NewtonSolver::computeNewtonCorrection(
      mfem::Vector& c,
      const mfem::Vector& r,
      const mfem::Vector& u,
      AdaptiveLinearSolverToleranceUpdater* const t) const noexcept {
    auto profiler = this->ctx_ptr != nullptr
                        ? this->ctx_ptr->startNewProfiling(
                              "NS::computeNewtonCorrection",
                              this->ctx_ptr->isProfilingEnabled())
                        : mgis::ProfilingSection{};
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    MFEM_ASSERT(this->prec != nullptr,
                "the Solver is not set (use setLinearSolver).");
    // the assembly of the tangent operator, the setup of the linear solver
    // (and of its preconditioner) and the resolution are timed separately.
    // The failures reported on all the processes, for instance by the
    // computations by kernels, make the computation of the correction fail
    // (see `callAndHandleFailures`).
    const auto success = this->callAndHandleFailures([this, &c, &r, &u, t] {
      auto& J = this->getJacobian(u);
      {
        auto profiler_setup = this->ctx_ptr != nullptr
                                  ? this->ctx_ptr->startNewProfiling(
                                        "NS::setLinearSolverOperator",
                                        this->ctx_ptr->isProfilingEnabled())
                                  : mgis::ProfilingSection{};
        this->prec->SetOperator(J);
      }
      auto profiler_mfem =
          this->ctx_ptr != nullptr
              ? this->ctx_ptr->startNewProfiling(
                    "MFEM::Mult(r,c)", this->ctx_ptr->isProfilingEnabled())
              : mgis::ProfilingSection{};
      this->prec->Mult(r, c);  // c = [DF(x_i)]^{-1} [F(x_i)-b]
      if ((t != nullptr) && (t->requiresLinearResidualNorm())) {
        // norm of the residual of the linear solve, J c - r
        auto lr = mfem::Vector(r.Size());
        lr.UseDevice(true);
        J.Mult(c, lr);
        lr -= r;
        t->setLinearResidualNorm(this->Norm(lr));
      }
      return true;
    });
    if (!success) {
      return false;
    }
    clearHypreErrors();
    if (this->discardLinearSolverFailure) {
      return true;
    }
    return hasConverged(*(this->prec));
  }  // end of computeNewtonCorrection

  NewtonSolver::~NewtonSolver() = default;

}  // end of namespace mfem_mgis
