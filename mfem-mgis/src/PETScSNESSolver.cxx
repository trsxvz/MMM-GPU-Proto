/*!
 * \file   src/PETScSNESSolver.cxx
 * \brief  This file implements the `PETScSNESSolver` class
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC)

#include <cmath>
#include <tuple>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <iomanip>
#include <utility>
#include <string_view>
#include <algorithm>
#include "petscsnes.h"
#include "mfem/general/device.hpp"
#include "mfem/linalg/hypre.hpp"
#include "mfem/linalg/petsc.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "PETScSupport.hxx"
#include "PETScSNESSolver.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief the linear solver of the problem, seen by PETSc as a shell
   * preconditioner, built by `ProblemLinearSolverFactory` at each setup of
   * the preconditioner by PETSc. The operator of the linear solver is set
   * when the tangent operator is computed
   * (`PETScSNESSolverImplementation::GetGradient`). A failure of the linear
   * solver is reported to PETSc as a failure of the preconditioner, which
   * stops the resolution, unless the failures of the linear solver are
   * discarded.
   */
  struct ProblemLinearSolverPreconditioner final : public mfem::Solver {
    //! \param[in] i: implementation of the non linear solver
    explicit ProblemLinearSolverPreconditioner(
        const PETScSNESSolverImplementation& i);
    //! \brief the operator of the linear solver is set elsewhere
    void SetOperator(const mfem::Operator&) override {}
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override;
    //! \brief not available: reported as a failure of the preconditioner
    void MultTranspose(const mfem::Vector& b, mfem::Vector& x) const override;

   private:
    //! \brief implementation of the non linear solver
    const PETScSNESSolverImplementation& implementation;
    //! \brief right hand side
    mutable mfem::Vector rhs;
    //! \brief solution
    mutable mfem::Vector solution;
  };  // end of struct ProblemLinearSolverPreconditioner

  //! \brief factory of the shell preconditioners of PETSc
  struct ProblemLinearSolverFactory final
      : public mfem::PetscPreconditionerFactory {
    //! \param[in] i: implementation of the non linear solver
    explicit ProblemLinearSolverFactory(const PETScSNESSolverImplementation& i)
        : mfem::PetscPreconditionerFactory("linear solver of mfem-mgis"),
          implementation(i) {}
    mfem::Solver* NewPreconditioner(const mfem::OperatorHandle&) override {
      return new ProblemLinearSolverPreconditioner(this->implementation);
    }

   private:
    //! \brief implementation of the non linear solver
    const PETScSNESSolverImplementation& implementation;
  };  // end of struct ProblemLinearSolverFactory

  /*!
   * \brief tangent operator seen by PETSc as a shell matrix: the vectors of
   * PETSc are copied into vectors of mfem-mgis on the host (see
   * `PETScSNESSolverImplementation::copyFromPETSc`)
   */
  struct TangentOperatorForPETSc final : public mfem::Operator {
    /*!
     * \param[in] s: non linear solver (which handles the failures of the
     * products, see `callAndHandleFailures`: the products are called by
     * PETSc, and fill the result with NaN on failure)
     */
    explicit TangentOperatorForPETSc(const PETScSNESSolver& s)
        : mfem::Operator(s.Height()), solver(s) {}
    //! \param[in] op: tangent operator
    void setOperator(const mfem::Operator& op) noexcept { this->tangent = &op; }
    void Mult(const mfem::Vector& px, mfem::Vector& py) const override;
    /*!
     * \brief product by the transposed tangent operator (trust region
     * methods of SNES), only available for an assembled tangent operator:
     * otherwise, the result is filled with NaN, which stops SNES
     */
    void MultTranspose(const mfem::Vector& px, mfem::Vector& py) const override;

   private:
    //! \brief non linear solver
    const PETScSNESSolver& solver;
    //! \brief tangent operator
    const mfem::Operator* tangent = nullptr;
    //! \brief copy of the argument
    mutable mfem::Vector x;
    //! \brief result
    mutable mfem::Vector y;
  };  // end of struct TangentOperatorForPETSc

  /*!
   * \brief the operator of the non linear problem given to the SNES of MFEM
   * (residual and tangent operator) and the SNES itself, built for one
   * resolution. The functions called by PETSc report the failures to SNES
   * (domain errors, failure of the preconditioner) instead of throwing
   * through it (see `callAndHandleFailures`).
   */
  struct PETScSNESSolverImplementation final : public mfem::Operator {
    /*!
     * \param[in] s: non linear solver
     * \param[in] c: communicator
     */
    PETScSNESSolverImplementation(const PETScSNESSolver& s, MPI_Comm c)
        : mfem::Operator(s.Height()),
          solver(s),
          communicator(c),
          options(c, "snes"),
          factory(*this),
          identity(s.Height()),
          tangent_operator(s) {}
    /*!
     * \brief create the SNES of MFEM, with the options of PETSc
     * \param[in] t: tolerance on the norm of the residual
     * \param[in] u0: initial guess
     * \param[in] r0: residual at the initial guess, given to SNES at its
     * first evaluation of the residual, instead of integrating the
     * behaviours again
     * \param[in] norm: norm of the residual at the initial guess
     */
    void initialize(const real t,
                    const mfem::Vector& u0,
                    const mfem::Vector& r0,
                    const real norm) {
      const auto& s = this->solver;
      auto ctx = Context{};
      if (!this->options.setOptions(ctx, s.options)) {
        throw CollectiveFailure("SNES: " + ctx.getErrorMessage());
      }
      this->snes = std::make_unique<mfem::PetscNonlinearSolver>(
          this->communicator, *this, this->options.prefix);
      // the criterion on the length of the step is disabled: newtonls
      // reports a convergence on it when its line search fails
      check(SNESSetTolerances(*(this->snes), PETSC_CURRENT, PETSC_CURRENT, 0,
                              s.max_iter, PETSC_CURRENT),
            "SNESSetTolerances");
      if (s.isLinearSolverFailureDiscarded()) {
        check(SNESSetMaxLinearSolveFailures(
                  *(this->snes), std::numeric_limits<PetscInt>::max()),
              "SNESSetMaxLinearSolveFailures");
      }
      if (s.petsc_linear_solver) {
        this->snes->SetJacobianType(mfem::Operator::PETSC_MATAIJ);
      } else {
        this->snes->SetJacobianType(mfem::Operator::PETSC_MATSHELL);
        this->snes->SetPreconditionerFactory(&(this->factory));
        // the shell preconditioner of the linear solver of the problem is
        // applied once, unless the options give a Krylov solver (set in the
        // code, not in the options database, where the methods of SNES
        // without linear solver would leave it unread)
        auto ksp = KSP{};
        check(SNESGetKSP(*(this->snes), &ksp), "SNESGetKSP");
        check(KSPSetType(ksp, KSPPREONLY), "KSPSetType");
      }
      check(SNESSetConvergenceTest(*(this->snes), checkConvergence, this,
                                   nullptr),
            "SNESSetConvergenceTest");
      // the domain errors of the tangent operator are only checked by
      // default by the debug builds of PETSc
      check(SNESSetCheckJacobianDomainError(*(this->snes), PETSC_TRUE),
            "SNESSetCheckJacobianDomainError");
      this->tolerance = t;
      this->initial_unknowns = u0;
      this->initial_residual = r0;
      this->initial_residual_available = true;
      this->window_norm = norm;
      this->last_norm = norm;
    }  // end of initialize
    /*!
     * \brief solve the non linear problem
     * \return the reason of the end of the resolution given by SNES
     * \param[in, out] x: initial guess, then solution
     * \param[out] iterations: number of iterations
     */
    SNESConvergedReason solve(mfem::Vector& x, int& iterations) {
      this->options.markAsUsed();
      // the vectors of SNES are vectors of PETSc of the host, whatever the
      // device of MFEM: the functions called by PETSc copy them into
      // vectors of mfem-mgis (see `copyFromPETSc`). On a device, the vectors
      // of PETSc of the GPU (memory of MFEM managed by Umpire) were not
      // handled by MFEM 4.10: `mfem::PetscNonlinearSolver::Mult`, which
      // places the memory of the unknowns in the vector of PETSc, left it in
      // a state unknown to the memory manager of MFEM, and so did the
      // destruction of a vector of MFEM wrapping a vector of PETSc whose
      // arrays PETSc changed; the trust region method of SNES and the
      // Krylov solvers of PETSc failed. SNES is thus called on a vector of
      // the host owned by PETSc, written and read through its array
      const auto size = x.Size();
      auto y = PETScVector{};
      check(VecCreate(this->communicator, &(y.v)), "VecCreate");
      check(VecSetSizes(y.v, size, PETSC_DECIDE), "VecSetSizes");
      check(VecSetType(y.v, VECSTANDARD), "VecSetType");
      copyToPETScVector(y.v, x);
      // the options are read here (as `mfem::PetscSolver::Customize`, which
      // stops the program on an invalid option)
      check(SNESSetFromOptions(*(this->snes)),
            "SNESSetFromOptions (invalid option of PETSc ?)");
      // the trust region methods multiply the residual by the transposed
      // tangent operator, which is not available without assembly (see
      // `GetGradient`)
      auto type = SNESType{};
      check(SNESGetType(*(this->snes), &type), "SNESGetType");
      this->trust_region = std::string_view{type}.starts_with(SNESNEWTONTR);
      check(SNESSolve(*(this->snes), nullptr, y.v), "SNESSolve");
      copyFromPETScVector(x, y.v);
      auto reason = SNESConvergedReason{};
      auto n = PetscInt{};
      check(SNESGetConvergedReason(*(this->snes), &reason),
            "SNESGetConvergedReason");
      check(SNESGetIterationNumber(*(this->snes), &n),
            "SNESGetIterationNumber");
      iterations = static_cast<int>(n);
      return reason;
    }  // end of solve
    /*!
     * \brief compute the residual, called by SNES. A failure of the
     * integration of the behaviours is reported to SNES as a domain error,
     * the residual being filled with NaN (a norm which no criterion
     * accepts)
     * \param[in] u: unknowns
     * \param[out] r: residual
     */
    void Mult(const mfem::Vector& pu, mfem::Vector& pr) const override {
      try {
        copyFromPETSc(pu, this->unknowns);
        const auto& u = this->unknowns;
        this->last_unknowns = u;
        this->last_evaluation_succeeded = false;
        // first evaluation, at the initial guess (collective test)
        if (this->initial_residual_available) {
          this->initial_residual_available = false;
          if (!(this->getMaximumDifference(u, this->initial_unknowns) > 0)) {
            copyToPETSc(this->initial_residual, pr);
            this->last_evaluation_succeeded = true;
            return;
          }
        }
        const auto& s = this->solver;
        auto& r = this->residual;
        r.UseDevice(true);
        r.SetSize(u.Size());
        const auto success = s.callAndHandleFailures([&s, &u, &r] {
          if (!s.processNewUnknownsEstimate(u)) {
            return false;
          }
          s.computeResidual(r, u);
          return true;
        });
        ++(this->number_of_evaluations);
        if (success) {
          copyToPETSc(r, pr);
          this->last_evaluation_succeeded = true;
          return;
        }
        if (s.shallPrint()) {
          s.getLogStream() << "SNES, evaluation " << std::setw(2)
                           << this->number_of_evaluations
                           << " : integration failure\n";
        }
      } catch (...) {
      }
      std::ignore = SNESSetFunctionDomainError(*(this->snes));
      auto* const r = pr.HostWrite();
      std::fill(r, r + pr.Size(), std::numeric_limits<real>::quiet_NaN());
    }  // end of Mult
    /*!
     * \brief compute the tangent operator, called by SNES, and give it to
     * the linear solver of the problem, if used. A failure is reported to
     * SNES as a domain error, the identity being returned.
     * \param[in] u: unknowns
     */
    mfem::Operator& GetGradient(const mfem::Vector& pu) const override {
      // the linear solvers of PETSc need the assembled operator, converted
      // by MFEM; otherwise, the products by the tangent operator made by
      // PETSc (shell matrix) are given to `tangent_operator`, which copies
      // the vectors of PETSc
      try {
        copyFromPETSc(pu, this->jacobian_unknowns);
        const auto& u = this->jacobian_unknowns;
        const auto& s = this->solver;
        // the tangent operator is the one of the last integration of the
        // behaviours, which may have been made elsewhere (last trial point
        // of a line search stopped on a small step, for instance)
        const auto integrated =
            this->last_evaluation_succeeded &&
            (!(this->getMaximumDifference(u, this->last_unknowns) > 0));
        auto* j = static_cast<mfem::Operator*>(nullptr);
        const auto tr = this->trust_region;
        const auto success =
            s.callAndHandleFailures([&s, &u, &j, integrated, tr] {
              if ((!integrated) && (!s.processNewUnknownsEstimate(u))) {
                return false;
              }
              j = &(s.getJacobian(u));
              if ((tr) && (dynamic_cast<mfem::HypreParMatrix*>(j) == nullptr)) {
                throw CollectiveFailure(
                    "SNES: the trust region methods require an assembled "
                    "tangent operator");
              }
              if (s.petsc_linear_solver) {
                if (dynamic_cast<mfem::HypreParMatrix*>(j) == nullptr) {
                  throw CollectiveFailure(
                      "SNES: the linear solvers of PETSc require an assembled "
                      "tangent operator");
                }
              } else {
                s.prec->SetOperator(*j);
              }
              return true;
            });
        ++(this->number_of_setups);
        if (!integrated) {
          this->last_unknowns = u;
          this->last_evaluation_succeeded = success;
        }
        if (success) {
          if (s.petsc_linear_solver) {
            if (!mfem::Device::IsEnabled()) {
              return *j;
            }
            // when MFEM computes on a GPU, the tangent operator is converted
            // to the AIJ format of PETSc of the host here (MFEM chooses the
            // format of the GPU), the solvers of PETSc computing on the host
            // (see the description of PETScLinearSolvers.cxx)
            this->petsc_tangent_operator =
                std::make_unique<mfem::PetscParMatrix>(
                    this->communicator, j, mfem::Operator::PETSC_MATAIJ);
            convertToHostMatrix(*(this->petsc_tangent_operator));
            return *(this->petsc_tangent_operator);
          }
          this->tangent_operator.setOperator(*j);
          return this->tangent_operator;
        }
      } catch (...) {
      }
      std::ignore = SNESSetJacobianDomainError(*(this->snes));
      return this->identity;
    }  // end of GetGradient
    /*!
     * \brief set or reset the failure of the shell preconditioner of the
     * linear solver of the problem, as a failure of the preconditioner of
     * PETSc
     * \param[in] f: if true, failure
     */
    void setPreconditionerFailure(const bool f) const noexcept {
      auto ksp = KSP{};
      auto pc = PC{};
      if ((SNESGetKSP(*(this->snes), &ksp) == 0) && (KSPGetPC(ksp, &pc) == 0)) {
        std::ignore = PCSetFailedReason(pc, f ? PC_SUBPC_ERROR : PC_NOERROR);
      }
    }  // end of setPreconditionerFailure
    //! \return if the linear solver of the problem is applied once
    [[nodiscard]] bool isPreconditionerAppliedOnce() const noexcept {
      auto ksp = KSP{};
      auto type = KSPType{};
      return (SNESGetKSP(*(this->snes), &ksp) == 0) &&
             (KSPGetType(ksp, &type) == 0) && (type != nullptr) &&
             (std::string_view{type} == KSPPREONLY);
    }  // end of isPreconditionerAppliedOnce
    /*!
     * \return the maximum absolute value of the difference of two vectors,
     * on all the processes
     */
    [[nodiscard]] real getMaximumDifference(const mfem::Vector& a,
                                            const mfem::Vector& b) const {
      // the maximum is reduced on all the processes, even if the sizes
      // differ on some of them
      if (a.Size() != b.Size()) {
        auto d = mfem::Vector(1);
        d = std::numeric_limits<real>::max();
        return this->getMaximumNorm(d);
      }
      auto d = mfem::Vector(a.Size());
      d.UseDevice(true);
      subtract(a, b, d);
      return this->getMaximumNorm(d);
    }  // end of getMaximumDifference
    //! \return the maximum norm of a vector, on all the processes
    [[nodiscard]] real getMaximumNorm(const mfem::Vector& v) const {
      auto n = static_cast<real>(v.Normlinf());
      auto l = n;
      MPI_Allreduce(&l, &n, 1, mfem::MPITypeMap<real>::mpi_type, MPI_MAX,
                    this->communicator);
      return n;
    }  // end of getMaximumNorm
    /*!
     * \brief copy a vector given by MFEM for PETSc (the vectors given to the
     * functions called by PETSc wrap the arrays of PETSc) into a vector of
     * mfem-mgis, on the host: `__mfem_snes_jacobian` (MFEM 4.10) gives to
     * `GetGradient` the array of the host of PETSc, unknown to the memory
     * manager of MFEM, which stops the program when it is read on a device
     * \param[in] src: vector given by MFEM
     * \param[out] dst: vector of mfem-mgis
     */
    static void copyFromPETSc(const mfem::Vector& src, mfem::Vector& dst) {
      dst.UseDevice(true);
      dst.SetSize(src.Size());
      const auto* const s = src.HostRead();
      std::copy(s, s + src.Size(), dst.HostWrite());
    }  // end of copyFromPETSc
    /*!
     * \brief copy a vector of mfem-mgis into a vector given by MFEM for
     * PETSc, on the host
     * \param[in] src: vector of mfem-mgis
     * \param[out] dst: vector given by MFEM
     */
    static void copyToPETSc(const mfem::Vector& src, mfem::Vector& dst) {
      const auto* const s = src.HostRead();
      std::copy(s, s + src.Size(), dst.HostWrite());
    }  // end of copyToPETSc
    //! \brief non linear solver
    const PETScSNESSolver& solver;
    //! \brief unknowns of the last evaluation of the residual
    mutable mfem::Vector last_unknowns;
    //! \brief if false, the last evaluation of the residual failed
    mutable bool last_evaluation_succeeded = true;
    //! \brief number of evaluations of the residual
    mutable int number_of_evaluations = 0;
    //! \brief number of updates of the tangent operator
    mutable int number_of_setups = 0;
    //! \brief if true, the stagnation of the residual stopped SNES
    bool stagnated = false;
    //! \brief norm of the residual given to the last convergence test
    real last_norm = 0;

   private:
    /*!
     * \brief throw a failure, on all the processes, if a function of PETSc
     * failed
     * \param[in] e: error code
     * \param[in] n: name of the function
     */
    static void check(const PetscErrorCode e, const char* const n) {
      if (e != 0) {
        throw CollectiveFailure("SNES: " + std::string{n} + " failed");
      }
    }  // end of check
    /*!
     * \brief convergence test of SNES (see `PETScSNESSolver`)
     * \return 0
     */
    static PetscErrorCode checkConvergence(SNES,
                                           PetscInt it,
                                           PetscReal,
                                           PetscReal,
                                           PetscReal fnorm,
                                           SNESConvergedReason* reason,
                                           void* const d) noexcept {
      auto& i = *static_cast<PETScSNESSolverImplementation*>(d);
      const auto& s = i.solver;
      const auto norm = static_cast<real>(fnorm);
      i.last_norm = norm;
      *reason = SNES_CONVERGED_ITERATING;
      // the initial residual (iteration 0), already known, and the trial
      // steps rejected by the trust region methods, which test the
      // convergence at the same iteration, are not counted
      const auto new_iteration = (it != 0) && (it > i.last_iteration);
      if (new_iteration) {
        i.last_iteration = it;
      }
      try {
        if (new_iteration) {
          auto iteration = Parameters{};
          iteration.replaceOrInsert("Norm", norm);
          s.iterations_information.push_back(iteration);
        }
        if (s.shallPrint()) {
          s.getLogStream() << "SNES, iteration " << std::setw(2) << it
                           << " : ||r|| = " << norm << ", ||r||/||r_0|| = "
                           << norm / (*(s.reference_residual_norm)) << '\n';
        }
      } catch (...) {
      }
      if (!std::isfinite(norm)) {
        *reason = SNES_DIVERGED_FNORM_NAN;
        return 0;
      }
      if (norm <= i.tolerance) {
        *reason = SNES_CONVERGED_FNORM_ABS;
        return 0;
      }
      // stagnation of the residual at the floor of the rounding errors:
      // three iterations without a decrease of more than 10 %, below the
      // stagnation tolerance
      constexpr auto window_decrease = real{9} / 10;
      if (norm < window_decrease * i.window_norm) {
        i.window_norm = norm;
        i.iterations_without_decrease = 0;
      } else if (new_iteration) {
        ++(i.iterations_without_decrease);
      }
      if ((i.iterations_without_decrease >= 3) &&
          (norm <=
           s.stagnation_relative_tolerance * (*(s.reference_residual_norm)))) {
        i.stagnated = true;
        *reason = SNES_CONVERGED_FNORM_ABS;
        return 0;
      }
      if (it >= s.max_iter) {
        *reason = SNES_DIVERGED_MAX_IT;
      }
      return 0;
    }  // end of checkConvergence
    //! \brief communicator
    const MPI_Comm communicator;
    /*!
     * \brief options of PETSc (declared before the SNES, which is destroyed
     * before the options are removed from the options database of PETSc)
     */
    PETScOptions options;
    //! \brief factory of the shell preconditioners
    ProblemLinearSolverFactory factory;
    //! \brief operator returned when the tangent operator can't be computed
    mutable mfem::IdentityOperator identity;
    //! \brief tangent operator given to PETSc as a shell matrix
    mutable TangentOperatorForPETSc tangent_operator;
    /*!
     * \brief tangent operator converted to the AIJ format of PETSc of the
     * host, when MFEM computes on a GPU and the linear solvers of PETSc are
     * used
     */
    mutable std::unique_ptr<mfem::PetscParMatrix> petsc_tangent_operator;

    //! \brief SNES of MFEM
    std::unique_ptr<mfem::PetscNonlinearSolver> snes;
    //! \brief tolerance on the norm of the residual
    real tolerance = 0;
    //! \brief initial guess
    mfem::Vector initial_unknowns;
    //! \brief residual at the initial guess
    mfem::Vector initial_residual;
    //! \brief if true, the initial residual is given to SNES
    mutable bool initial_residual_available = false;
    //! \brief unknowns given by SNES to the evaluation of the residual
    mutable mfem::Vector unknowns;
    //! \brief residual
    mutable mfem::Vector residual;
    //! \brief unknowns given by SNES to the computation of the tangent
    mutable mfem::Vector jacobian_unknowns;
    /*!
     * \brief norm of the residual at the beginning of the current window of
     * iterations without a decrease of more than 10 % (detection of a
     * stagnation of the residual)
     */
    real window_norm = std::numeric_limits<real>::max();
    //! \brief number of iterations of the current window
    int iterations_without_decrease = 0;
    //! \brief last iteration given to the convergence test
    PetscInt last_iteration = 0;
    //! \brief if true, the method of SNES is a trust region method
    bool trust_region = false;
  };  // end of struct PETScSNESSolverImplementation

  //! \brief fill a vector given by MFEM for PETSc with NaN
  static void fillWithNaN(mfem::Vector& v) noexcept {
    auto* const p = v.HostWrite();
    std::fill(p, p + v.Size(), std::numeric_limits<real>::quiet_NaN());
  }  // end of fillWithNaN

  void TangentOperatorForPETSc::Mult(const mfem::Vector& px,
                                     mfem::Vector& py) const {
    try {
      const auto success = this->solver.callAndHandleFailures([this, &px] {
        PETScSNESSolverImplementation::copyFromPETSc(px, this->x);
        this->y.UseDevice(true);
        this->y.SetSize(this->Height());
        this->tangent->Mult(this->x, this->y);
        return true;
      });
      if (success) {
        PETScSNESSolverImplementation::copyToPETSc(this->y, py);
        return;
      }
    } catch (...) {
    }
    fillWithNaN(py);
  }  // end of Mult

  void TangentOperatorForPETSc::MultTranspose(const mfem::Vector& px,
                                              mfem::Vector& py) const {
    try {
      const auto* const a =
          dynamic_cast<const mfem::HypreParMatrix*>(this->tangent);
      if (a != nullptr) {
        const auto success = this->solver.callAndHandleFailures([this, a, &px] {
          PETScSNESSolverImplementation::copyFromPETSc(px, this->x);
          this->y.UseDevice(true);
          this->y.SetSize(this->Width());
          a->MultTranspose(this->x, this->y);
          return true;
        });
        if (success) {
          PETScSNESSolverImplementation::copyToPETSc(this->y, py);
          return;
        }
      }
    } catch (...) {
    }
    fillWithNaN(py);
  }  // end of MultTranspose

  ProblemLinearSolverPreconditioner::ProblemLinearSolverPreconditioner(
      const PETScSNESSolverImplementation& i)
      : mfem::Solver(i.Height()), implementation(i) {}

  void ProblemLinearSolverPreconditioner::Mult(const mfem::Vector& pb,
                                               mfem::Vector& px) const {
    // the failure of a previous application is forgotten (a shell
    // preconditioner of PETSc keeps it)
    this->implementation.setPreconditionerFailure(false);
    auto failed = true;
    try {
      // copies of the vectors of PETSc (see
      // `PETScSNESSolverImplementation::copyFromPETSc`)
      PETScSNESSolverImplementation::copyFromPETSc(pb, this->rhs);
      this->solution.UseDevice(true);
      this->solution.SetSize(pb.Size());
      const auto& s = this->implementation.solver;
      const auto& b = this->rhs;
      auto& x = this->solution;
      const auto success = s.callAndHandleFailures([&s, &b, &x] {
        s.prec->Mult(b, x);
        return true;
      });
      PETScSNESSolverImplementation::copyToPETSc(x, px);
      clearHypreErrors();
      failed = (!success) || ((!s.isLinearSolverFailureDiscarded()) &&
                              (!hasConverged(*(s.prec))));
    } catch (...) {
    }
    // a failure stops the resolution if the linear solver of the problem is
    // applied once; a Krylov solver of PETSc which it preconditions checks
    // its own convergence
    if (failed && this->implementation.isPreconditionerAppliedOnce()) {
      this->implementation.setPreconditionerFailure(true);
    }
  }  // end of Mult

  void ProblemLinearSolverPreconditioner::MultTranspose(
      const mfem::Vector&, mfem::Vector& px) const {
    // the linear solvers of mfem-mgis are not transposed: the Krylov
    // solvers of PETSc which need it (bicg, cgne, lsqr...) fail
    this->implementation.setPreconditionerFailure(true);
    fillWithNaN(px);
  }  // end of MultTranspose

  PETScSNESSolver::PETScSNESSolver(
      NonLinearEvolutionProblemImplementation<true>& p)
      : NonLinearSolverBase(p) {
    initializePETSc();
  }  // end of PETScSNESSolver

  PETScSNESSolver::PETScSNESSolver(
      NonLinearEvolutionProblemImplementation<false>& p)
      : NonLinearSolverBase(p) {
    mgis::raise("SNES: the non linear solver is only available in parallel");
  }  // end of PETScSNESSolver

  void PETScSNESSolver::Mult(const mfem::Vector&, mfem::Vector& x) const {
    MFEM_ASSERT(this->petsc_linear_solver || (this->prec != nullptr),
                "the Solver is not set (use setLinearSolver).");
    this->iterations_information.clear();
    this->converged = 0;
    this->final_iter = 0;
    this->final_norm = std::numeric_limits<real>::max();
    auto& log = this->getLogStream();
    const auto shall_print = this->shallPrint();
    // initial residual
    auto r = mfem::Vector(this->Height());
    r.UseDevice(true);
    const auto initialized = this->callAndHandleFailures([this, &x, &r] {
      if (!this->processNewUnknownsEstimate(x)) {
        return false;
      }
      this->computeResidual(r, x);
      return true;
    });
    if (!initialized) {
      return;
    }
    const auto norm = this->Norm(r);
    if (!this->reference_residual_norm.has_value()) {
      this->reference_residual_norm = norm;
    }
    this->initial_norm = *(this->reference_residual_norm);
    if (shall_print) {
      log << "SNES, initial residual : ||r|| = " << norm
          << ", ||r||/||r_0|| = " << norm / (*(this->reference_residual_norm))
          << '\n';
    }
    const auto tolerance =
        std::max(this->rel_tol * (*(this->reference_residual_norm)),
                 static_cast<real>(this->abs_tol));
    if (norm <= tolerance) {
      this->converged = 1;
      this->final_norm = norm;
      return;
    }
    // resolution
    auto i = std::unique_ptr<PETScSNESSolverImplementation>{};
    auto reason = SNESConvergedReason{};
    auto iterations = int{};
    const auto solved = this->callAndHandleFailures([&] {
      i = std::make_unique<PETScSNESSolverImplementation>(*this,
                                                          this->GetComm());
      i->initialize(tolerance, x, r, norm);
      reason = i->solve(x, iterations);
      return true;
    });
    if (!solved) {
      return;
    }
    this->final_iter = iterations;
    // the behaviours are integrated again if the last evaluation of the
    // residual failed or was not made at the solution (a trial point of the
    // line search)
    const auto evaluated_at_solution =
        (i->last_evaluation_succeeded) &&
        (!(i->getMaximumDifference(x, i->last_unknowns) > 0));
    if (!evaluated_at_solution) {
      const auto updated = this->callAndHandleFailures([this, &x, &r] {
        if (!this->processNewUnknownsEstimate(x)) {
          return false;
        }
        this->computeResidual(r, x);
        return true;
      });
      if (!updated) {
        return;
      }
      this->final_norm = this->Norm(r);
    } else {
      // the residual of the last evaluation is the one of the solution:
      // its norm is the last one given to the convergence test
      this->final_norm = i->last_norm;
    }
    const auto stagnation =
        this->final_norm <= this->stagnation_relative_tolerance *
                                (*(this->reference_residual_norm));
    // the positive reasons of SNES which the convergence test does not give
    // (for instance the convergence on the length of the step of the
    // method ksponly) are not accepted if the residual is too large
    this->converged = (this->final_norm <= tolerance) || stagnation;
    if (shall_print) {
      log << "SNES : " << SNESConvergedReasons[reason] << " after "
          << iterations << " iterations (" << i->number_of_setups
          << " updates of the tangent operator"
          << (i->stagnated ? ", stagnation of the residual" : "")
          << "), ||r|| = " << this->final_norm << " (tolerance " << tolerance
          << ")" << (this->converged ? "" : ", not converged") << "\n\n";
    }
  }  // end of Mult

  bool PETScSNESSolver::setSolverParameters(Context& ctx,
                                            const Parameters& params) noexcept {
    const auto own = std::vector<std::string>{"LinearSolver", "Options",
                                              "StagnationRelativeTolerance"};
    // the options are read and checked before being applied: a refused call
    // leaves the solver unchanged
    auto new_petsc_linear_solver = this->petsc_linear_solver;
    if (contains(params, "LinearSolver")) {
      const auto ov = get<std::string>(ctx, params, "LinearSolver");
      if (isInvalid(ov)) {
        return false;
      }
      if ((*ov != "Problem") && (*ov != "PETSc")) {
        return ctx.registerErrorMessage("SNES: invalid linear solver '" + *ov +
                                        "' (expected 'Problem' or 'PETSc')");
      }
      new_petsc_linear_solver = *ov == "PETSc";
    }
    auto new_options = this->options;
    if (contains(params, "Options")) {
      const auto ov = get<Parameters>(ctx, params, "Options");
      if (isInvalid(ov)) {
        return false;
      }
      new_options = *ov;
    }
    for (const auto& [n, p] : new_options) {
      if (n.empty() || n.starts_with('-')) {
        return ctx.registerErrorMessage(
            "SNES: invalid option '" + n +
            "' of PETSc (the options are given without the leading dash)");
      }
      if (n.starts_with("snes_convergence_test")) {
        return ctx.registerErrorMessage(
            "SNES: the option '" + n +
            "' is refused (the convergence test is the one of mfem-mgis)");
      }
      if ((!new_petsc_linear_solver) && (n.starts_with("pc_"))) {
        return ctx.registerErrorMessage(
            "SNES: the option '" + n +
            "' of the preconditioner of PETSc is refused with the linear "
            "solver of the problem (option LinearSolver set to 'Problem')");
      }
    }
    const auto ostagnation = get_if<double>(
        ctx, params, "StagnationRelativeTolerance",
        static_cast<double>(this->stagnation_relative_tolerance));
    if (isInvalid(ostagnation)) {
      return false;
    }
    if ((*ostagnation < 0) || (*ostagnation >= 1)) {
      return ctx.registerErrorMessage(
          "SNES: the option StagnationRelativeTolerance must be in [0, 1[");
    }
    const auto oiterations =
        get_if<int>(ctx, params, "MaximumNumberOfIterations", this->max_iter);
    if (isInvalid(oiterations)) {
      return false;
    }
    if (*oiterations < 1) {
      return ctx.registerErrorMessage(
          "SNES: the maximum number of iterations must be strictly positive");
    }
    if (!NonLinearSolverBase::setSolverParameters(ctx, remove(params, own))) {
      return false;
    }
    this->petsc_linear_solver = new_petsc_linear_solver;
    this->options = std::move(new_options);
    this->stagnation_relative_tolerance = static_cast<real>(*ostagnation);
    return true;
  }  // end of setSolverParameters

  PETScSNESSolver::~PETScSNESSolver() = default;

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC) */
