/*!
 * \file   src/MFEMLinearSolvers.cxx
 * \brief  This file implements the module of the linear solvers and of the
 * preconditioners of MFEM: Krylov solvers, UMFPACK (SuiteSparse), Jacobi
 * and Chebyshev smoothers, and the `AssembledOperator` preconditioner of
 * mfem-mgis (see `LinearSolverModules.hxx`)
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <memory>
#include <string>
#include <vector>
#include <optional>
#include <utility>
#include <functional>
#include <string_view>
#include "mfem/config/config.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Raise.hxx"
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
#ifdef MFEM_USE_SUITESPARSE
      declare(f, "UMFPackSolver",
              LinearSolverDescription{
                  .library = "SuiteSparse",
                  .kind = LinearSolverKind::DIRECT,
                  .placement = LinearSolverPlacement::HOST,
                  .preconditioners = AcceptedPreconditioners::NONE,
                  .options = {},
                  .description = "UMFPACK direct solver of SuiteSparse"},
              [](Context&, FiniteElementSpace<false>&, const Parameters&) {
                return LinearSolverHandler{
                    std::make_unique<mfem::UMFPackSolver>(),
                    std::unique_ptr<LinearSolverPreconditioner>{}};
              });
#endif /* MFEM_USE_SUITESPARSE */
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
  }  // end of declareMFEMPreconditioners

#ifdef MFEM_USE_MPI
  template void declareMFEMLinearSolvers<true>(LinearSolverFactory<true>&);
  template void declareMFEMPreconditioners<true>(PreconditionerFactory<true>&);
#endif /* MFEM_USE_MPI */
  template void declareMFEMLinearSolvers<false>(LinearSolverFactory<false>&);
  template void declareMFEMPreconditioners<false>(
      PreconditionerFactory<false>&);

}  // end of namespace mfem_mgis::internals
