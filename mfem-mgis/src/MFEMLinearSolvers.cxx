/*!
 * \file   src/MFEMLinearSolvers.cxx
 * \brief  This file implements the module of the linear solvers and of the
 * preconditioners of MFEM: Krylov solvers, UMFPACK and KLU (SuiteSparse),
 * Jacobi and Chebyshev smoothers, smoothers and incomplete factorization of
 * sparse matrices (Gauss-Seidel, Jacobi, block ILU), the orthogonalizing
 * wrapper of a solver, and the `AssembledOperator` preconditioner of
 * mfem-mgis (see `LinearSolverModules.hxx`)
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <tuple>
#include <memory>
#include <string>
#include <algorithm>
#include <vector>
#include <optional>
#include <utility>
#include <functional>
#include <string_view>
#include "mfem/config/config.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/linalg/sparsesmoothers.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/MatrixFreeTangentOperator.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis {

  namespace {

    /*!
     * \brief Jacobi smoother of MFEM (`mfem::OperatorJacobiSmoother`), whose
     * diagonal is computed by the operator (`AssembleDiagonal`), or by
     * `GetDiag` for a sparse matrix of MFEM, which does not implement
     * `AssembleDiagonal`. The diagonal entries of the essential degrees of
     * freedom are the ones of the operator. The smoother of MFEM is kept
     * from an operator to the next one of the same size: only its inverse
     * diagonal is computed again.
     */
    struct JacobiSmoother final : mfem::Solver {
      //! \param[in] d: damping factor
      explicit JacobiSmoother(const real d) : damping(d) {}
      void SetOperator(const mfem::Operator& op) override {
        this->height = op.Height();
        this->width = op.Width();
        this->diagonal.UseDevice(true);
        if (const auto* const m = dynamic_cast<const mfem::SparseMatrix*>(&op);
            m != nullptr) {
          m->GetDiag(this->diagonal);
        } else {
          this->diagonal.SetSize(op.Height());
          op.AssembleDiagonal(this->diagonal);
        }
        if ((this->smoother == nullptr) ||
            (this->smoother->Height() != op.Height())) {
          this->smoother = std::make_unique<mfem::OperatorJacobiSmoother>(
              this->diagonal, this->essential_dofs, this->damping);
        } else {
          this->smoother->Setup(this->diagonal);
        }
      }
      void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
        if (this->smoother == nullptr) {
          raise("JacobiSmoother::Mult: the operator is not set");
        }
        this->smoother->Mult(x, y);
      }

     private:
      //! \brief damping factor
      const real damping;
      //! \brief diagonal of the operator
      mfem::Vector diagonal;
      /*!
       * \brief essential degrees of freedom, empty (the operator is
       * constrained), referenced by the smoother of MFEM
       */
      const mfem::Array<int> essential_dofs;
      //! \brief smoother of MFEM
      std::unique_ptr<mfem::OperatorJacobiSmoother> smoother;
    };

    /*!
     * \brief Chebyshev smoother of MFEM (`mfem::OperatorChebyshevSmoother`):
     * polynomial of the given order in the operator scaled by the inverse of
     * its diagonal, computed as for the Jacobi smoother (`JacobiSmoother`).
     * The largest eigenvalue of the scaled operator is estimated by the power
     * method at each call to `SetOperator`. Only products by the operator
     * are used: usable with an operator applied without assembly (see
     * `MatrixFreeTangentOperator`), whose essential degrees of freedom are
     * eliminated with a diagonal set to one.
     */
    struct ChebyshevSmoother final : mfem::Solver {
#ifdef MFEM_USE_MPI
      /*!
       * \param[in] o: order of the polynomial
       * \param[in] n: number of iterations of the power method
       * \param[in] c: communicator (`MPI_COMM_NULL` in sequential)
       */
      ChebyshevSmoother(const int o, const int n, MPI_Comm c)
          : order(o), power_iterations(n), comm(c) {}
#else  /* MFEM_USE_MPI */
      /*!
       * \param[in] o: order of the polynomial
       * \param[in] n: number of iterations of the power method
       */
      ChebyshevSmoother(const int o, const int n)
          : order(o), power_iterations(n) {}
#endif /* MFEM_USE_MPI */
      void SetOperator(const mfem::Operator& op) override {
        this->height = op.Height();
        this->width = op.Width();
        this->diagonal.UseDevice(true);
        if (const auto* const m = dynamic_cast<const mfem::SparseMatrix*>(&op);
            m != nullptr) {
          m->GetDiag(this->diagonal);
        } else {
          this->diagonal.SetSize(op.Height());
          op.AssembleDiagonal(this->diagonal);
        }
        // the smoother of MFEM estimates the largest eigenvalue at its
        // construction and references the diagonal
        this->smoother.reset();
#ifdef MFEM_USE_MPI
        this->smoother = std::make_unique<mfem::OperatorChebyshevSmoother>(
            op, this->diagonal, this->essential_dofs, this->order, this->comm,
            this->power_iterations);
#else  /* MFEM_USE_MPI */
        this->smoother = std::make_unique<mfem::OperatorChebyshevSmoother>(
            op, this->diagonal, this->essential_dofs, this->order,
            this->power_iterations);
#endif /* MFEM_USE_MPI */
      }
      void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
        if (this->smoother == nullptr) {
          raise("ChebyshevSmoother::Mult: the operator is not set");
        }
        this->smoother->Mult(x, y);
      }

     private:
      //! \brief order of the polynomial
      const int order;
      //! \brief number of iterations of the power method
      const int power_iterations;
#ifdef MFEM_USE_MPI
      //! \brief communicator of the power method
      MPI_Comm comm;
#endif /* MFEM_USE_MPI */
      //! \brief diagonal of the operator
      mfem::Vector diagonal;
      /*!
       * \brief essential degrees of freedom, empty (the operator is
       * constrained), referenced by the smoother of MFEM
       */
      const mfem::Array<int> essential_dofs;
      //! \brief smoother of MFEM
      std::unique_ptr<mfem::OperatorChebyshevSmoother> smoother;
    };

    /*!
     * \brief preconditioner built on the assembled version of an operator
     * applied without assembly (see `AssembledOperatorProvider`), set up
     * every `period` calls to `SetOperator`, and kept between them.
     *
     * The preconditioner is also set up again when the assembled operator
     * on which it was set up is no longer valid, i.e. when the operator has
     * been assembled again in between (see
     * `AssembledOperatorProvider::getAssembledOperatorRevision`).
     *
     * An operator which can't give an assembled version (for instance the
     * matrix of the prediction) keeps the preconditioner set up on a valid
     * assembled operator, if any: it is still a preconditioner of a close
     * operator. Otherwise, the preconditioner is set up on this operator,
     * and set up again at the next call (this operator may be destroyed
     * after its resolution).
     */
    struct AssembledOperatorPreconditioner final : mfem::Solver {
      /*!
       * \param[in] p: preconditioner
       * \param[in] n: period of the updates
       */
      AssembledOperatorPreconditioner(std::unique_ptr<mfem::Solver> p,
                                      const int n)
          : preconditioner(std::move(p)), period(n) {}
      void SetOperator(const mfem::Operator& op) override {
        this->height = op.Height();
        this->width = op.Width();
        const auto* const provider =
            dynamic_cast<const AssembledOperatorProvider*>(&op);
        if (provider == nullptr) {
          if (this->isSetUpOnAValidAssembledOperator()) {
            return;
          }
          this->preconditioner->SetOperator(op);
          this->revision.reset();
          this->setup_revision.reset();
          return;
        }
        const auto r = provider->getAssembledOperatorRevision();
        const auto valid =
            (r == this->revision) && this->isSetUpOnAValidAssembledOperator();
        if ((!valid) || (this->count % this->period == 0)) {
          this->preconditioner->SetOperator(provider->getAssembledOperator());
          this->revision = r;
          this->setup_revision = *r;
          this->count = 0;
        }
        ++(this->count);
      }
      void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
        this->preconditioner->Mult(x, y);
      }

     private:
      /*!
       * \return if the preconditioner is set up on an assembled operator which
       * has not been assembled again since
       */
      [[nodiscard]] bool isSetUpOnAValidAssembledOperator() const noexcept {
        return (this->revision != nullptr) &&
               (this->setup_revision.has_value()) &&
               (*(this->revision) == *(this->setup_revision));
      }
      //! \brief preconditioner
      std::unique_ptr<mfem::Solver> preconditioner;
      //! \brief period of the updates
      const int period;
      //! \brief number of calls to `SetOperator` since the last setup
      int count = 0;
      //! \brief number of assemblies of the assembled operator
      std::shared_ptr<const size_type> revision;
      /*!
       * \brief number of assemblies of the assembled operator when the
       * preconditioner was set up on it, if it was
       */
      std::optional<size_type> setup_revision;
    };

    /*!
     * \brief solver or preconditioner applying a solver of MFEM defined on a
     * sparse matrix: to the sparse matrix of MFEM in sequential, and to the
     * diagonal block of the matrix of hypre owned by each process in
     * parallel (block Jacobi method, one block per process).
     *
     * These solvers compute on the host, and some of them read the memory
     * of the host without synchronizing it with the device:
     * `mfem::KLUSolver` sorts the columns of its matrix (on the device when
     * MFEM computes on a GPU) then reads the arrays of the host,
     * `mfem::BlockILU` and `mfem::KLUSolver` read and write the vectors on
     * the host. They are given a copy of the matrix which only lives on the
     * host (which also leaves the matrix of mfem-mgis or of hypre
     * unchanged, the direct solvers of SuiteSparse sorting the columns of
     * their matrix), and copies of the vectors on the host.
     */
    struct HostSparseMatrixSolver final : mfem::Solver {
      /*!
       * \param[in] p: solver
       * \param[in] n: name of the solver, for the error messages
       */
      HostSparseMatrixSolver(std::unique_ptr<mfem::Solver> p, std::string n)
          : solver(std::move(p)), name(std::move(n)) {}
      void SetOperator(const mfem::Operator& op) override {
        const auto* m = dynamic_cast<const mfem::SparseMatrix*>(&op);
#ifdef MFEM_USE_MPI
        const auto* const a = dynamic_cast<const mfem::HypreParMatrix*>(&op);
        auto block = mfem::SparseMatrix{};
        if (a != nullptr) {
          a->GetDiag(block);
          m = &block;
        }
#endif /* MFEM_USE_MPI */
        if (m == nullptr) {
          raise("the solver '" + this->name +
                "' requires an assembled operator");
        }
        // a copy on the host only, whose columns are sorted on the host
        this->matrix = std::make_unique<mfem::SparseMatrix>(
            *m, true, mfem::MemoryType::HOST);
        this->matrix->UseGPUSparse(false);
        std::ignore = this->matrix->HostReadWriteI();
        std::ignore = this->matrix->HostReadWriteJ();
        std::ignore = this->matrix->HostReadWriteData();
        this->height = op.Height();
        this->width = op.Width();
        this->solver->SetOperator(*(this->matrix));
        // the factorization of KLU does not report its failures (a singular
        // matrix): all the processes fail together
        auto failed = 0;
#ifdef MFEM_MGIS_HAVE_SUITESPARSE
        if (const auto* const klu =
                dynamic_cast<const mfem::KLUSolver*>(this->solver.get());
            klu != nullptr) {
          failed = klu->Common.status != KLU_OK ? 1 : 0;
        }
#endif /* MFEM_MGIS_HAVE_SUITESPARSE */
#ifdef MFEM_USE_MPI
        if (a != nullptr) {
          MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_MAX,
                        a->GetComm());
        }
#endif /* MFEM_USE_MPI */
        if (failed != 0) {
          throw CollectiveFailure(
              "the factorization of the solver '" + this->name +
              "' failed (singular matrix or lack of memory)");
        }
      }
      void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
        this->apply(x, y, false);
      }
      void MultTranspose(const mfem::Vector& x,
                         mfem::Vector& y) const override {
        this->apply(x, y, true);
      }

     private:
      /*!
       * \brief apply the solver, or its transpose, to vectors copied on the
       * host
       * \param[in] x: vector
       * \param[out] y: result
       * \param[in] transpose: if true, the transpose of the solver is applied
       */
      void apply(const mfem::Vector& x,
                 mfem::Vector& y,
                 const bool transpose) const {
        this->host_x.SetSize(x.Size());
        this->host_y.SetSize(y.Size());
        std::copy_n(x.HostRead(), x.Size(), this->host_x.HostWrite());
        this->host_y = 0.;
        if (transpose) {
          this->solver->MultTranspose(this->host_x, this->host_y);
        } else {
          this->solver->Mult(this->host_x, this->host_y);
        }
        std::copy_n(this->host_y.HostRead(), y.Size(), y.HostWrite());
      }
      /*!
       * \brief copy of the matrix on the host, declared before the solver,
       * which references it
       */
      std::unique_ptr<mfem::SparseMatrix> matrix;
      //! \brief solver
      std::unique_ptr<mfem::Solver> solver;
      //! \brief name of the solver
      const std::string name;
      //! \brief copies of the vectors on the host
      mutable mfem::Vector host_x, host_y;
    };

  }  // end of namespace

  /*!
   * \brief set the parameters of an iterative solver of MFEM: tolerances,
   * maximum number of iterations, verbosity, dimension of the Krylov space
   * of the GMRES solvers (`KDim`) and preconditioner
   * \return the preconditioner (a null pointer if none is requested), or an
   * empty value on error
   * \param[in, out] ctx: execution context
   * \param[in, out] s: solver
   * \param[in] fespace: finite element space
   * \param[in] n: name of the solver
   * \param[in] d: description of the solver
   * \param[in] params: parameters of the solver
   */
  template <bool parallel>
  static std::optional<std::unique_ptr<LinearSolverPreconditioner>>
  setLinearSolverParameters(Context& ctx,
                            mfem::IterativeSolver& s,
                            FiniteElementSpace<parallel>& fespace,
                            const std::string_view n,
                            const LinearSolverDescription& d,
                            const Parameters& params) {
    const char* const Preconditioner = "Preconditioner";
    const char* const Krylov_Dimension = "KDim";
    const auto is_gmres_solver =
        (dynamic_cast<mfem::GMRESSolver*>(&s) != nullptr) ||
        (dynamic_cast<mfem::FGMRESSolver*>(&s) != nullptr);
    s.iterative_mode = false;
    auto allowed_parameters = getIterativeSolverParametersList();
    allowed_parameters.push_back(Preconditioner);
    if (is_gmres_solver) {
      allowed_parameters.push_back(Krylov_Dimension);
    }
    if (!checkParameters(ctx, params, allowed_parameters)) {
      return {};
    }
    const auto oparams =
        extract(ctx, params, getIterativeSolverParametersList());
    if (isInvalid(oparams)) {
      return {};
    }
    if (!setSolverParameters(ctx, s, *oparams)) {
      return {};
    }
    if ((is_gmres_solver) && (contains(params, Krylov_Dimension))) {
      const auto odim = get<int>(ctx, params, Krylov_Dimension);
      if (isInvalid(odim)) {
        return {};
      }
      if (*odim < 1) {
        std::ignore = ctx.registerErrorMessage(
            "invalid Krylov dimension (must be greater than one)");
        return {};
      }
      if (auto* const ptr = dynamic_cast<mfem::GMRESSolver*>(&s);
          ptr != nullptr) {
        ptr->SetKDim(*odim);
      }
      if (auto* const ptr = dynamic_cast<mfem::FGMRESSolver*>(&s);
          ptr != nullptr) {
        ptr->SetKDim(*odim);
      }
    }
    auto oprec =
        internals::makePreconditioner<parallel>(ctx, fespace, n, d, params);
    if (isInvalid(oprec)) {
      return {};
    }
    if (*oprec != nullptr) {
      s.SetPreconditioner(**oprec);
    }
    return oprec;
  }  // end of setLinearSolverParameters

  /*!
   * \return the generator of an iterative solver of MFEM
   * \param[in] n: name of the solver
   * \param[in] d: description of the solver
   */
  template <bool parallel, typename LinearSolverType>
  static typename LinearSolverFactory<parallel>::Generator
  buildIterativeSolverGenerator(const std::string n,
                                const LinearSolverDescription d) {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      return [n, d](Context& ctx, FiniteElementSpace<true>& fespace,
                    const Parameters& params) {
        auto s = std::make_unique<LinearSolverType>(fespace.GetComm());
        auto oprec =
            setLinearSolverParameters<true>(ctx, *s, fespace, n, d, params);
        if (isInvalid(oprec)) {
          return LinearSolverHandler{};
        }
        return LinearSolverHandler{std::move(s), std::move(*oprec)};
      };
#else  /* MFEM_USE_MPI */
      return {};
#endif /* MFEM_USE_MPI */
    } else {
      return [n, d](Context& ctx, FiniteElementSpace<false>& fespace,
                    const Parameters& params) {
        auto s = std::make_unique<LinearSolverType>();
        auto oprec =
            setLinearSolverParameters<false>(ctx, *s, fespace, n, d, params);
        if (isInvalid(oprec)) {
          return LinearSolverHandler{};
        }
        return LinearSolverHandler{std::move(s), std::move(*oprec)};
      };
    }
  }  // end of buildIterativeSolverGenerator

  /*!
   * \return the Jacobi smoother of MFEM, whose diagonal is computed by the
   * operator (`AssembleDiagonal`): usable with an operator applied without
   * assembly (see `MatrixFreeTangentOperator`)
   * \param[in, out] ctx: execution context
   * \param[in] opts: options (`Damping`)
   */
  template <bool parallel>
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  makeOperatorJacobiSmoother(Context& ctx,
                             FiniteElementSpace<parallel>&,
                             const Parameters& opts) {
    if (!checkParameters(ctx, opts, {"Damping"})) {
      return {};
    }
    const auto odamping = get_if<double>(ctx, opts, "Damping", 1.);
    if (isInvalid(odamping)) {
      return {};
    }
    return std::make_unique<JacobiSmoother>(*odamping);
  }  // end of makeOperatorJacobiSmoother

  /*!
   * \return the Chebyshev smoother of MFEM: usable with an operator applied
   * without assembly
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options (`Order`, `NumberOfPowerIterations`)
   */
  template <bool parallel>
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  makeOperatorChebyshevSmoother(
      Context& ctx,
      [[maybe_unused]] FiniteElementSpace<parallel>& fespace,
      const Parameters& opts) {
    if (!checkParameters(
            ctx, opts,
            std::vector<std::string>{"Order", "NumberOfPowerIterations"})) {
      return {};
    }
    const auto oorder = get_if<int>(ctx, opts, "Order", 2);
    const auto oiterations =
        get_if<int>(ctx, opts, "NumberOfPowerIterations", 10);
    if (isInvalid(oorder) || isInvalid(oiterations)) {
      return {};
    }
    // orders implemented by MFEM
    if ((*oorder < 1) || (*oorder > 5)) {
      return ctx.registerErrorMessage(
          "setLinearSolverPreconditioner: the order of the Chebyshev "
          "smoother must be between 1 and 5");
    }
    if (*oiterations < 1) {
      return ctx.registerErrorMessage(
          "setLinearSolverPreconditioner: the number of iterations of the "
          "power method of the Chebyshev smoother must be strictly "
          "positive");
    }
#ifdef MFEM_USE_MPI
    auto comm = MPI_Comm{MPI_COMM_NULL};
    if constexpr (parallel) {
      comm = fespace.GetComm();
    }
    return std::make_unique<ChebyshevSmoother>(*oorder, *oiterations, comm);
#else  /* MFEM_USE_MPI */
    return std::make_unique<ChebyshevSmoother>(*oorder, *oiterations);
#endif /* MFEM_USE_MPI */
  }  // end of makeOperatorChebyshevSmoother

  /*!
   * \return the preconditioner built on the assembled version of an operator
   * applied without assembly (see `MatrixFreeTangentOperator`)
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options (`Preconditioner`, the wrapped preconditioner,
   * and `UpdatePeriod`)
   */
  template <bool parallel>
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  makeAssembledOperatorPreconditioner(Context& ctx,
                                      FiniteElementSpace<parallel>& fespace,
                                      const Parameters& opts) {
    if (!checkParameters(
            ctx, opts,
            std::vector<std::string>{"Preconditioner", "UpdatePeriod"})) {
      return {};
    }
    const auto oprec = get<Parameters>(ctx, opts, "Preconditioner");
    const auto operiod = get_if<int>(ctx, opts, "UpdatePeriod", 1);
    if ((isInvalid(oprec)) || (isInvalid(operiod))) {
      return {};
    }
    if (*operiod < 1) {
      return ctx.registerErrorMessage(
          "setLinearSolverPreconditioner: invalid update period of the "
          "'AssembledOperator' preconditioner");
    }
    auto p = internals::generatePreconditioner<parallel>(ctx, fespace, *oprec);
    if (p == nullptr) {
      return ctx.registerErrorMessage(
          "setLinearSolverPreconditioner: no preconditioner given to the "
          "'AssembledOperator' preconditioner");
    }
    return std::make_unique<AssembledOperatorPreconditioner>(std::move(p),
                                                             *operiod);
  }  // end of makeAssembledOperatorPreconditioner

  /*!
   * \return the Gauss-Seidel smoother of MFEM
   * \param[in, out] ctx: execution context
   * \param[in] opts: options (`Type`: `Symmetric`, by default, `Forward` or
   * `Backward`; `NumberOfIterations`, one by default)
   */
  template <bool parallel>
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  makeGSSmoother(Context& ctx,
                 FiniteElementSpace<parallel>&,
                 const Parameters& opts) {
    using GSSmoother = mfem::GSSmoother;
    if (!checkParameters(
            ctx, opts,
            std::vector<std::string>{"Type", "NumberOfIterations"})) {
      return {};
    }
    const auto otype =
        get_if<std::string>(ctx, opts, "Type", std::string{"Symmetric"});
    const auto oiterations = get_if<int>(ctx, opts, "NumberOfIterations", 1);
    if (isInvalid(otype) || isInvalid(oiterations)) {
      return {};
    }
    auto type = GSSmoother::SYMMETRIC;
    if (*otype == "Forward") {
      type = GSSmoother::FORWARD;
    } else if (*otype == "Backward") {
      type = GSSmoother::BACKWARD;
    } else if (*otype != "Symmetric") {
      return ctx.registerErrorMessage("invalid type '" + *otype +
                                      "' of the Gauss-Seidel smoother");
    }
    if (*oiterations < 1) {
      return ctx.registerErrorMessage(
          "invalid number of iterations of the Gauss-Seidel smoother");
    }
    return std::make_unique<HostSparseMatrixSolver>(
        std::make_unique<GSSmoother>(type, *oiterations), "GSSmoother");
  }  // end of makeGSSmoother

  /*!
   * \return the Jacobi smoother of the sparse matrices of MFEM
   * \param[in, out] ctx: execution context
   * \param[in] opts: options (`Type`: `Jacobi`, by default, `l1Jacobi` or
   * `LumpedJacobi`; `Damping`, one by default; `NumberOfIterations`, one by
   * default)
   */
  template <bool parallel>
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  makeDSmoother(Context& ctx,
                FiniteElementSpace<parallel>&,
                const Parameters& opts) {
    using DSmoother = mfem::DSmoother;
    if (!checkParameters(ctx, opts,
                         std::vector<std::string>{"Type", "Damping",
                                                  "NumberOfIterations"})) {
      return {};
    }
    const auto otype =
        get_if<std::string>(ctx, opts, "Type", std::string{"Jacobi"});
    const auto odamping = get_if<double>(ctx, opts, "Damping", 1.);
    const auto oiterations = get_if<int>(ctx, opts, "NumberOfIterations", 1);
    if (isInvalid(otype) || isInvalid(odamping) || isInvalid(oiterations)) {
      return {};
    }
    auto type = DSmoother::JACOBI;
    if (*otype == "l1Jacobi") {
      type = DSmoother::L1_JACOBI;
    } else if (*otype == "LumpedJacobi") {
      type = DSmoother::LUMPED_JACOBI;
    } else if (*otype != "Jacobi") {
      return ctx.registerErrorMessage("invalid type '" + *otype +
                                      "' of the Jacobi smoother");
    }
    if (*oiterations < 1) {
      return ctx.registerErrorMessage(
          "invalid number of iterations of the Jacobi smoother");
    }
    if (*odamping <= 0) {
      return ctx.registerErrorMessage(
          "invalid damping factor of the Jacobi smoother");
    }
    return std::make_unique<HostSparseMatrixSolver>(
        std::make_unique<DSmoother>(type, *odamping, *oiterations),
        "DSmoother");
  }  // end of makeDSmoother

  /*!
   * \return the block incomplete LU factorization of MFEM, without fill
   * (the only level of fill supported by MFEM)
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options (`BlockSize`, the number of components of the
   * unknowns, which must be ordered by vector dimension, by default;
   * `Reordering`: `MinimumDiscardedFill`, by default, or `None`)
   */
  template <bool parallel>
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner> makeBlockILU(
      Context& ctx,
      FiniteElementSpace<parallel>& fespace,
      const Parameters& opts) {
    using BlockILU = mfem::BlockILU;
    if (!checkParameters(ctx, opts,
                         std::vector<std::string>{"BlockSize", "Reordering"})) {
      return {};
    }
    // the blocks of the unknowns ordered by vector dimension
    const auto default_block_size =
        fespace.GetOrdering() == mfem::Ordering::byVDIM ? fespace.GetVDim() : 1;
    const auto osize = get_if<int>(ctx, opts, "BlockSize", default_block_size);
    const auto oreordering = get_if<std::string>(
        ctx, opts, "Reordering", std::string{"MinimumDiscardedFill"});
    if (isInvalid(osize) || isInvalid(oreordering)) {
      return {};
    }
    auto reordering = BlockILU::Reordering::MINIMUM_DISCARDED_FILL;
    if (*oreordering == "None") {
      reordering = BlockILU::Reordering::NONE;
    } else if (*oreordering != "MinimumDiscardedFill") {
      return ctx.registerErrorMessage("invalid reordering '" + *oreordering +
                                      "' of the block ILU factorization");
    }
    if (*osize < 1) {
      return ctx.registerErrorMessage(
          "invalid block size of the block ILU factorization");
    }
    return std::make_unique<HostSparseMatrixSolver>(
        std::make_unique<BlockILU>(*osize, reordering), "BlockILU");
  }  // end of makeBlockILU

  /*!
   * \return the wrapper of a preconditioner which orthogonalizes the
   * vectors to the constant vector (`mfem::OrthoSolver`, which computes on
   * the host), for operators whose kernel is the constant vector (scalar
   * problems with pure Neumann boundary conditions)
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options (`Preconditioner`, the wrapped preconditioner)
   */
  template <bool parallel>
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  makeOrthoSolver(Context& ctx,
                  FiniteElementSpace<parallel>& fespace,
                  const Parameters& opts) {
    if (!checkParameters(ctx, opts,
                         std::vector<std::string>{"Preconditioner"})) {
      return {};
    }
    const auto oprec = get<Parameters>(ctx, opts, "Preconditioner");
    if (isInvalid(oprec)) {
      return {};
    }
    const auto oname = get<std::string>(ctx, *oprec, "Name");
    if (isInvalid(oname)) {
      return {};
    }
    if (*oname == "None") {
      return ctx.registerErrorMessage(
          "setLinearSolverPreconditioner: the 'OrthoSolver' preconditioner "
          "requires a preconditioner");
    }
    // the error of the generation, if any, is already registered
    auto p = internals::generatePreconditioner<parallel>(ctx, fespace, *oprec);
    if (p == nullptr) {
      return {};
    }
    // the orthogonalizing solver of MFEM references the wrapped solver
    struct Wrapper final : mfem::Solver {
#ifdef MFEM_USE_MPI
      Wrapper(std::unique_ptr<mfem::Solver> s, MPI_Comm c)
          : solver(std::move(s)), ortho(c) {
        this->ortho.SetSolver(*(this->solver));
      }
#endif /* MFEM_USE_MPI */
      explicit Wrapper(std::unique_ptr<mfem::Solver> s) : solver(std::move(s)) {
        this->ortho.SetSolver(*(this->solver));
      }
      void SetOperator(const mfem::Operator& op) override {
        this->ortho.SetOperator(op);
        this->height = op.Height();
        this->width = op.Width();
      }
      // as a preconditioner, from a null initial guess (`iterative_mode` of
      // the orthogonalizing solver of MFEM, false by default)
      void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
        this->ortho.Mult(x, y);
      }
      std::unique_ptr<mfem::Solver> solver;
      mfem::OrthoSolver ortho;
    };
#ifdef MFEM_USE_MPI
    if constexpr (parallel) {
      return std::make_unique<Wrapper>(std::move(p), fespace.GetComm());
    }
#endif /* MFEM_USE_MPI */
    return std::make_unique<Wrapper>(std::move(p));
  }  // end of makeOrthoSolver

}  // end of namespace mfem_mgis

namespace mfem_mgis::internals {

  template <bool parallel>
  void declareMFEMLinearSolvers(LinearSolverFactory<parallel>& f) {
    auto options = getIterativeSolverParametersList();
    options.push_back("Preconditioner");
    auto gmres_options = options;
    gmres_options.push_back("KDim");
    const auto krylov = [](std::vector<std::string> o, std::string d) {
      return LinearSolverDescription{
          .library = "MFEM",
          .kind = LinearSolverKind::ITERATIVE,
          .placement = LinearSolverPlacement::MFEM,
          .preconditioners = AcceptedPreconditioners::ANY,
          .options = std::move(o),
          .description = std::move(d)};
    };
    const auto add = [&f]<typename SolverType>(const std::string& n,
                                               LinearSolverDescription d) {
      declare(f, n, d,
              buildIterativeSolverGenerator<parallel, SolverType>(n, d));
    };
    add.template operator()<mfem::CGSolver>(
        "CGSolver", krylov(options, "conjugate gradient of MFEM"));
    add.template operator()<mfem::GMRESSolver>(
        "GMRESSolver", krylov(gmres_options, "GMRES of MFEM"));
    add.template operator()<mfem::FGMRESSolver>(
        "FGMRESSolver", krylov(gmres_options, "flexible GMRES of MFEM"));
    add.template operator()<mfem::BiCGSTABSolver>(
        "BiCGSTABSolver", krylov(options, "BiCGSTAB of MFEM"));
    add.template operator()<mfem::MINRESSolver>(
        "MINRESSolver", krylov(options, "MINRES of MFEM"));
    add.template operator()<mfem::SLISolver>(
        "SLISolver", krylov(options, "stationary linear iterations of MFEM"));
    if constexpr (!parallel) {
#ifdef MFEM_MGIS_HAVE_SUITESPARSE
      declare(f, "UMFPackSolver",
              LinearSolverDescription{
                  .library = "SuiteSparse",
                  .kind = LinearSolverKind::DIRECT,
                  .placement = LinearSolverPlacement::HOST,
                  .preconditioners = AcceptedPreconditioners::NONE,
                  .options = {},
                  .description = "UMFPACK direct solver of SuiteSparse"},
              [](Context& ctx, FiniteElementSpace<false>&,
                 const Parameters& params) -> LinearSolverHandler {
                if (!checkParameters(ctx, params, std::vector<std::string>{})) {
                  return {};
                }
                return LinearSolverHandler{
                    std::make_unique<mfem::UMFPackSolver>(),
                    std::unique_ptr<LinearSolverPreconditioner>{}};
              });
      declare(f, "KLUSolver",
              LinearSolverDescription{
                  .library = "SuiteSparse",
                  .kind = LinearSolverKind::DIRECT,
                  .placement = LinearSolverPlacement::HOST,
                  .preconditioners = AcceptedPreconditioners::NONE,
                  .options = {},
                  .description = "KLU direct solver of SuiteSparse (for "
                                 "circuit-like sparse matrices)"},
              [](Context& ctx, FiniteElementSpace<false>&,
                 const Parameters& params) -> LinearSolverHandler {
                if (!checkParameters(ctx, params, std::vector<std::string>{})) {
                  return {};
                }
                return LinearSolverHandler{
                    std::make_unique<HostSparseMatrixSolver>(
                        std::make_unique<mfem::KLUSolver>(), "KLUSolver"),
                    std::unique_ptr<LinearSolverPreconditioner>{}};
              });
#endif /* MFEM_MGIS_HAVE_SUITESPARSE */
    }
  }  // end of declareMFEMLinearSolvers

  template <bool parallel>
  void declareMFEMPreconditioners(PreconditionerFactory<parallel>& f) {
    const auto smoother = [](std::vector<std::string> o, std::string d) {
      return LinearSolverDescription{
          .library = "MFEM",
          .kind = LinearSolverKind::PRECONDITIONER,
          .placement = LinearSolverPlacement::MFEM,
          .preconditioners = AcceptedPreconditioners::NONE,
          .options = std::move(o),
          .description = std::move(d)};
    };
    declare(f, "OperatorJacobiSmoother",
            smoother({"Damping"}, "Jacobi smoother of MFEM"),
            makeOperatorJacobiSmoother<parallel>);
    declare(f, "OperatorChebyshevSmoother",
            smoother({"Order", "NumberOfPowerIterations"},
                     "Chebyshev smoother of MFEM"),
            makeOperatorChebyshevSmoother<parallel>);
    // a wrapper: the location of the wrapped preconditioner is checked
    declare(f, "AssembledOperator",
            LinearSolverDescription{
                .library = "mfem-mgis",
                .kind = LinearSolverKind::PRECONDITIONER,
                .placement = std::nullopt,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Preconditioner", "UpdatePeriod"},
                .description =
                    "preconditioner set up on the assembled version of an "
                    "operator applied without assembly"},
            makeAssembledOperatorPreconditioner<parallel>);
    // preconditioners of the sparse matrices of MFEM, applied to the diagonal
    // block of each process in parallel
    const auto sparse = [](const char* const library,
                           std::vector<std::string> o, std::string d) {
      return LinearSolverDescription{
          .library = library,
          .kind = LinearSolverKind::PRECONDITIONER,
          .placement = LinearSolverPlacement::HOST,
          .preconditioners = AcceptedPreconditioners::NONE,
          .options = std::move(o),
          .description = std::move(d)};
    };
    declare(f, "GSSmoother",
            sparse("MFEM", {"Type", "NumberOfIterations"},
                   "Gauss-Seidel smoother of MFEM, applied to the diagonal "
                   "block of each process in parallel"),
            makeGSSmoother<parallel>);
    declare(f, "DSmoother",
            sparse("MFEM", {"Type", "Damping", "NumberOfIterations"},
                   "Jacobi smoother of the sparse matrices of MFEM, applied "
                   "to the diagonal block of each process in parallel"),
            makeDSmoother<parallel>);
    declare(f, "BlockILU",
            sparse("MFEM", {"BlockSize", "Reordering"},
                   "block incomplete LU factorization of MFEM, applied to the "
                   "diagonal block of each process in parallel"),
            makeBlockILU<parallel>);
#ifdef MFEM_MGIS_HAVE_SUITESPARSE
    const auto direct = [](const char* const n) {
      return [n](Context& ctx, FiniteElementSpace<parallel>&,
                 const Parameters& opts)
                 -> std::unique_ptr<LinearSolverPreconditioner> {
        if (!checkParameters(ctx, opts, std::vector<std::string>{})) {
          return {};
        }
        auto s = std::unique_ptr<mfem::Solver>{};
        if (std::string_view{n} == "UMFPack") {
          s = std::make_unique<mfem::UMFPackSolver>();
        } else {
          s = std::make_unique<mfem::KLUSolver>();
        }
        return std::make_unique<HostSparseMatrixSolver>(std::move(s), n);
      };
    };
    declare(f, "UMFPack",
            sparse("SuiteSparse", {},
                   "UMFPACK direct solver of SuiteSparse, applied to the "
                   "diagonal block of each process in parallel"),
            direct("UMFPack"));
    declare(f, "KLU",
            sparse("SuiteSparse", {},
                   "KLU direct solver of SuiteSparse, applied to the diagonal "
                   "block of each process in parallel"),
            direct("KLU"));
#endif /* MFEM_MGIS_HAVE_SUITESPARSE */
    // the orthogonalization computes on the host (the location of the
    // wrapped preconditioner is checked at its generation)
    declare(f, "OrthoSolver",
            LinearSolverDescription{
                .library = "MFEM",
                .kind = LinearSolverKind::PRECONDITIONER,
                .placement = LinearSolverPlacement::HOST,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Preconditioner"},
                .description =
                    "wrapper of a preconditioner orthogonalizing the vectors "
                    "to the constant vector, on the host (operators whose "
                    "kernel is the constant vector)"},
            makeOrthoSolver<parallel>);
  }  // end of declareMFEMPreconditioners

#ifdef MFEM_USE_MPI
  template void declareMFEMLinearSolvers<true>(LinearSolverFactory<true>&);
  template void declareMFEMPreconditioners<true>(PreconditionerFactory<true>&);
#endif /* MFEM_USE_MPI */
  template void declareMFEMLinearSolvers<false>(LinearSolverFactory<false>&);
  template void declareMFEMPreconditioners<false>(
      PreconditionerFactory<false>&);

}  // end of namespace mfem_mgis::internals
