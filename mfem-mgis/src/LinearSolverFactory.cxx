/*!
 * \file   src/LinearSolverFactory.cxx
 * \brief  This file implements the `LinearSolverFactory` class
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <array>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <optional>
#include <utility>
#include <algorithm>
#include "mfem/linalg/solvers.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/linalg/petsc.hpp"
#include "mfem/config/config.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#ifdef MFEM_USE_MUMPS
#include "mfem/linalg/mumps.hpp"
#endif
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/MatrixFreeTangentOperator.hxx"
#include "MGIS/Contract.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "AmgXSolvers.hxx"
#include "MUMPSSolvers.hxx"

#ifdef MFEM_USE_MPI

namespace mfem_mgis::internals {

  /*!
   * \brief BoomerAMG preconditioner whose settings are applied again for each
   * new operator.
   *
   * `HypreBoomerAMG::SetOperator` creates a new solver of hypre, which only
   * keeps some options (`HypreBoomerAMG::ResetAMGPrecond`): the description
   * of the systems (`SetSystemsOptions`, which depends on the size of the
   * operator for unknowns ordered by nodes), the number of cycles per
   * application, the type of cycle and the number of sweeps per cycle are
   * lost. The settings are applied after each `SetOperator`, before the
   * setup of hypre, made at the first application.
   */
  struct ConfiguredHypreBoomerAMG final : public mfem::HypreBoomerAMG {
    //! \brief a setting of the preconditioner
    using Setting = std::function<void(mfem::HypreBoomerAMG&)>;
    //! \param[in] s: settings
    explicit ConfiguredHypreBoomerAMG(std::vector<Setting> s)
        : settings(std::move(s)) {}
    //! \brief set the operator, then the settings
    void SetOperator(const mfem::Operator& op) override {
      mfem::HypreBoomerAMG::SetOperator(op);
      for (const auto& setting : this->settings) {
        setting(*this);
      }
    }

   private:
    //! \brief settings
    const std::vector<Setting> settings;
  };  // end of struct ConfiguredHypreBoomerAMG

#if MFEM_HYPRE_VERSION >= 21900

  /*!
   * \brief incomplete factorisation of hypre whose settings are applied
   * again after each new operator: `mfem::HypreILU::SetOperator` creates a
   * new preconditioner of hypre with the default options of MFEM (a level
   * of fill of one), which lost the settings given by the user
   */
  struct ConfiguredHypreILU final : public mfem::HypreILU {
    //! \brief a setting of the preconditioner
    using Setting = std::function<void(mfem::HypreILU&)>;
    //! \param[in] s: settings
    explicit ConfiguredHypreILU(std::vector<Setting> s)
        : settings(std::move(s)) {
      this->applySettings();
    }
    //! \brief set the operator, then the settings
    void SetOperator(const mfem::Operator& op) override {
      mfem::HypreILU::SetOperator(op);
      this->applySettings();
    }

   private:
    //! \brief apply the settings
    void applySettings() {
      for (const auto& setting : this->settings) {
        setting(*this);
      }
    }
    //! \brief settings
    const std::vector<Setting> settings;
  };  // end of struct ConfiguredHypreILU

#endif /* MFEM_HYPRE_VERSION >= 21900 */

  /*!
   * \brief ParaSails preconditioner of hypre which can be given a new
   * operator after the destruction of the previous one, as the non linear
   * solvers do: `mfem::HypreParaSails::SetOperator` (MFEM 4.10) creates the
   * new preconditioner of hypre with the communicator of the previous
   * operator, whose memory may then be freed (invalid communicator, see
   * reproducers/hypre-parasails-operator in the prototype). The operator is
   * replaced by the new one before calling this method, which reports an
   * operator which is not a matrix of hypre as before.
   */
  struct ReusableHypreParaSails final : public mfem::HypreParaSails {
    using mfem::HypreParaSails::HypreParaSails;
    //! \brief set the operator
    void SetOperator(const mfem::Operator& op) override {
      const auto* const new_A = dynamic_cast<const mfem::HypreParMatrix*>(&op);
      if ((this->A != nullptr) && (new_A != nullptr)) {
        this->A = new_A;
      }
      mfem::HypreParaSails::SetOperator(op);
    }
  };  // end of struct ReusableHypreParaSails

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_USE_MPI */

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  /*!
   * \return the BoomerAMG preconditioner of hypre
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options:
   * - `Strategy`: `System` (systems options, by default for vector unknowns),
   *   `Elasticity` (elasticity options, requiring unknowns ordered by vector
   *   dimension, `System` being used otherwise) or `None`;
   * - `VerbosityLevel`: print level of hypre;
   * - options of hypre, whose values are the ones of hypre (see its
   *   documentation): `StrengthThreshold`, `CoarseningType`,
   *   `AggressiveCoarseningLevels`, `InterpolationType`, `RelaxationType`,
   *   `RelaxationSweeps` (sweeps before and after the coarse levels),
   *   `CycleType` (1: V-cycle, 2: W-cycle), `NumberOfCycles` (cycles per
   *   application, 1 by default) and `MaximumNumberOfLevels`.
   */
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  setHypreBoomerAMGPreconditioner(Context& ctx,
                                  FiniteElementSpace<true>& fespace,
                                  const Parameters& opts) noexcept {
    using Problem = AbstractNonLinearEvolutionProblem;
    using Setting = internals::ConfiguredHypreBoomerAMG::Setting;
    const auto integer_options = std::vector<std::string>{
        "CoarseningType", "AggressiveCoarseningLevels", "InterpolationType",
        "RelaxationType", "RelaxationSweeps",           "CycleType",
        "NumberOfCycles", "MaximumNumberOfLevels"};
    auto allowed = integer_options;
    allowed.insert(allowed.end(), {"Strategy", "StrengthThreshold",
                                   Problem::SolverVerbosityLevel});
    if (!checkParameters(ctx, opts, allowed)) {
      return {};
    }
    auto settings = std::vector<Setting>{};
    // strategy
    const auto ostrategy = get_if<std::string>(
        ctx, opts, "Strategy",
        std::string{fespace.GetVDim() > 1 ? "System" : "None"});
    if (isInvalid(ostrategy)) {
      return {};
    }
    auto elasticity = false;
    if (*ostrategy == "Elasticity") {
      // the elasticity version of BoomerAMG requires Ordering::byVDIM, and
      // is not supported by hypre when it computes on a GPU
      const auto by_vdim = fespace.GetOrdering() == mfem::Ordering::byVDIM;
      elasticity = by_vdim && (!isHypreUsingGPU());
      if ((!elasticity) && (fespace.GetMyRank() == 0)) {
        warning(ctx.log(),
                by_vdim ? "setLinearSolverParameters: strategy 'Elasticity' "
                          "of preconditioner HypreBoomerAMG is not supported "
                          "by hypre on a GPU, using strategy 'System' instead"
                        : "setLinearSolverParameters: strategy 'Elasticity' "
                          "of preconditioner HypreBoomerAMG requires unknowns "
                          "ordered by vector dimension (i.e. the parameter "
                          "'FiniteElementSpaceOrdering' must be set to "
                          "'byVDIM'), using strategy 'System' instead");
      }
    } else if ((*ostrategy != "System") && (*ostrategy != "None")) {
      return ctx.registerErrorMessage(
          "setLinearSolverParameters: "
          "invalid strategy '" +
          *ostrategy + "' for preconditioner HypreBoomerAMG");
    }
    if ((!elasticity) && (*ostrategy != "None")) {
      const auto n = fespace.GetVDim();
      const auto by_nodes = fespace.GetOrdering() == mfem::Ordering::byNODES;
      settings.push_back([n, by_nodes](mfem::HypreBoomerAMG& amg) {
        amg.SetSystemsOptions(n, by_nodes);
      });
    }
    // options of hypre
    if (contains(opts, "StrengthThreshold")) {
      const auto ov = get<double>(ctx, opts, "StrengthThreshold");
      if (isInvalid(ov)) {
        return {};
      }
      settings.push_back(
          [v = *ov](mfem::HypreBoomerAMG& amg) { amg.SetStrengthThresh(v); });
    }
    for (const auto& o : integer_options) {
      if (!contains(opts, o)) {
        continue;
      }
      const auto ov = get<int>(ctx, opts, o);
      if (isInvalid(ov)) {
        return {};
      }
      const auto v = *ov;
      if ((v < 0) || ((v == 0) && (o != "CoarseningType") &&
                      (o != "AggressiveCoarseningLevels") &&
                      (o != "InterpolationType") && (o != "RelaxationType"))) {
        return ctx.registerErrorMessage(
            "setLinearSolverParameters: invalid "
            "value for option '" +
            o + "' of HypreBoomerAMG");
      }
      if (o == "CoarseningType") {
        settings.push_back(
            [v](mfem::HypreBoomerAMG& amg) { amg.SetCoarsening(v); });
      } else if (o == "AggressiveCoarseningLevels") {
        settings.push_back(
            [v](mfem::HypreBoomerAMG& amg) { amg.SetAggressiveCoarsening(v); });
      } else if (o == "InterpolationType") {
        settings.push_back(
            [v](mfem::HypreBoomerAMG& amg) { amg.SetInterpolation(v); });
      } else if (o == "RelaxationType") {
        settings.push_back(
            [v](mfem::HypreBoomerAMG& amg) { amg.SetRelaxType(v); });
      } else if (o == "RelaxationSweeps") {
        settings.push_back(
            [v](mfem::HypreBoomerAMG& amg) { amg.SetCycleNumSweeps(v, v); });
      } else if (o == "CycleType") {
        settings.push_back(
            [v](mfem::HypreBoomerAMG& amg) { amg.SetCycleType(v); });
      } else if (o == "NumberOfCycles") {
        settings.push_back(
            [v](mfem::HypreBoomerAMG& amg) { amg.SetMaxIter(v); });
      } else {
        settings.push_back(
            [v](mfem::HypreBoomerAMG& amg) { amg.SetMaxLevels(v); });
      }
    }
    auto amg = std::make_unique<internals::ConfiguredHypreBoomerAMG>(settings);
    if (elasticity) {
      amg->SetElasticityOptions(&fespace);
    }
    if (contains(opts, Problem::SolverVerbosityLevel)) {
      const auto olvl = get<int>(ctx, opts, Problem::SolverVerbosityLevel);
      if (isInvalid(olvl)) {
        return {};
      }
      amg->SetPrintLevel(*olvl);
    }
    return amg;
  }  // end of setHypreBoomerAMGPreconditioner

  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  setHypreEuclidPreconditioner(Context& ctx,
                               FiniteElementSpace<true>& p,
                               const Parameters& opts) {
    using Problem = AbstractNonLinearEvolutionProblem;
    if (!checkParameters(ctx, opts, {Problem::SolverVerbosityLevel})) {
      return {};
    }
    auto euclid = std::make_unique<mfem::HypreEuclid>(p.GetComm());
    return euclid;
  }  // end of setHypreEuclidPreconditioner

  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  setHypreILUPreconditioner(Context& ctx,
                            FiniteElementSpace<true>&,
                            const Parameters& opts) {
#if MFEM_HYPRE_VERSION >= 21900
    using Problem = AbstractNonLinearEvolutionProblem;
    if (!checkParameters(
            ctx, opts,
            std::vector<std::string>{Problem::SolverVerbosityLevel,
                                     "HypreILUType", "HypreILULevelOfFill"})) {
      return {};
    }
    using Setting = internals::ConfiguredHypreILU::Setting;
    auto settings = std::vector<Setting>{};
    if (contains(opts, Problem::SolverVerbosityLevel)) {
      const auto olvl = get<int>(ctx, opts, Problem::SolverVerbosityLevel);
      if (isInvalid(olvl)) {
        return {};
      }
      settings.push_back(
          [v = *olvl](mfem::HypreILU& ilu) { ilu.SetPrintLevel(v); });
    }
    // type of the factorisation of hypre (`HYPRE_ILUSetType`): ILU(k)
    // (0, the default, 10, 20, 30, 40), ILUT (1, 11, 21, 31, 41), modified
    // ILU(0) (50), with a block Jacobi (0, 1), GMRES on the Schur complement
    // (10, 11, 40, 41), Newton-Schulz-Hotelling (20, 21), restricted additive
    // Schwarz (30, 31) or RAP (50) global method
    const auto otype = get_if<int>(ctx, opts, "HypreILUType", 0);
    if (isInvalid(otype)) {
      return {};
    }
    const auto types = std::array{0, 1, 10, 11, 20, 21, 30, 31, 40, 41, 50};
    if (std::find(types.begin(), types.end(), *otype) == types.end()) {
      return ctx.registerErrorMessage(
          "setHypreILUPreconditioner: invalid type of factorisation " +
          std::to_string(*otype));
    }
    const auto ilut = (*otype % 10) == 1;
    if (contains(opts, "HypreILUType")) {
      settings.push_back([v = *otype](mfem::HypreILU& ilu) { ilu.SetType(v); });
    }
    // the ILU(k) factorisations with a level of fill greater than zero make
    // the triangular solves of hypre fail when hypre computes on a GPU
    // (segmentation fault in hypre_ILUSolveLUDevice, reproduced with hypre
    // 3.2.0 and CUDA 13, see reproducers/hypre-ilu-gpu in the prototype),
    // the ILUT factorisations do not: the level of fill is then zero by
    // default, one being the default of MFEM
    auto level_of_fill = std::optional<int>{};
    if (contains(opts, "HypreILULevelOfFill")) {
      const auto olvl = get<int>(ctx, opts, "HypreILULevelOfFill");
      if (isInvalid(olvl)) {
        return {};
      }
      if (ilut) {
        return ctx.registerErrorMessage(
            "setHypreILUPreconditioner: no level of fill for an ILUT "
            "factorisation");
      }
      if ((*olvl > 0) && (isHypreUsingGPU())) {
        return ctx.registerErrorMessage(
            "setHypreILUPreconditioner: a level of fill greater than zero "
            "is not supported when hypre computes on a GPU (the ILUT "
            "factorisation, HypreILUType 1, is)");
      }
      level_of_fill = *olvl;
    } else if ((!ilut) && (isHypreUsingGPU())) {
      level_of_fill = 0;
    }
    if (level_of_fill.has_value()) {
      settings.push_back(
          [v = *level_of_fill](mfem::HypreILU& ilu) { ilu.SetLevelOfFill(v); });
    }
    auto ilu = std::make_unique<internals::ConfiguredHypreILU>(settings);
    return ilu;
#else  /*  HYPRE_OLD_VERSION */
    return ctx.registerErrorMessage(
        "Support for HypreILU is not available with this version of MFEM");
#endif /* HYPRE_OLD_VERSION */
  }    // end of setHypreILUPreconditioner

  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  setHypreParaSailsPreconditioner(Context& ctx,
                                  FiniteElementSpace<true>& p,
                                  const Parameters& opts) {
    using Problem = AbstractNonLinearEvolutionProblem;
    if (!checkParameters(ctx, opts, {Problem::SolverVerbosityLevel})) {
      return {};
    }
    // ParaSails only computes on the host: its setup reads the matrix stored
    // on the device when hypre computes on a GPU (segmentation fault in
    // MatrixSetRow with hypre 3.2.0 and CUDA 13)
    if (isHypreUsingGPU()) {
      return ctx.registerErrorMessage(
          "setHypreParaSailsPreconditioner: the ParaSails preconditioner is "
          "not supported when hypre computes on a GPU");
    }
    auto parasails =
        std::make_unique<internals::ReusableHypreParaSails>(p.GetComm());
    // statistics of the setups printed by hypre (printed by default by
    // MFEM), kept for the new operators
    if (contains(opts, Problem::SolverVerbosityLevel)) {
      const auto olvl = get<int>(ctx, opts, Problem::SolverVerbosityLevel);
      if (isInvalid(olvl)) {
        return {};
      }
      parasails->SetLogging(*olvl);
    }
    return parasails;
  }  // end of setHypreParaSailsPreconditioner

  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  setHypreDiagScalePreconditioner(Context& ctx,
                                  FiniteElementSpace<true>&,
                                  const Parameters& opts) {
    using Problem = AbstractNonLinearEvolutionProblem;
    if (!checkParameters(ctx, opts, {Problem::SolverVerbosityLevel})) {
      return {};
    }
    return std::make_unique<mfem::HypreDiagScale>();
  }   // end of setHypreDiagScalePreconditioner
#else /* MFEM_USE_MPI */

  [[noreturn]] std::unique_ptr<LinearSolverPreconditioner>
  setHypreBoomerAMGPreconditioner(Context&,
                                  FiniteElementSpace<true>&,
                                  const Parameters&) {
    reportUnsupportedParallelComputations();
  }  // end of setHypreBoomerAMGPreconditioner

  [[noreturn]] std::unique_ptr<LinearSolverPreconditioner>
  setHypreEuclidPreconditioner(Context&,
                               FiniteElementSpace<true>&,
                               const Parameters&) {
    reportUnsupportedParallelComputations();
  }  // end of setHypreEuclidPreconditioner

  [[noreturn]] std::unique_ptr<LinearSolverPreconditioner>
  setHypreILUPreconditioner(Context&,
                            FiniteElementSpace<true>&,
                            const Parameters&) {
    reportUnsupportedParallelComputations();
  }  // end of setHypreILUPreconditioner

  [[noreturn]] std::unique_ptr<LinearSolverPreconditioner>
  setHypreParaSailsPreconditioner(Context&,
                                  FiniteElementSpace<true>&,
                                  const Parameters&) {
    reportUnsupportedParallelComputations();
  }  // end of setHypreParaSailsPreconditioner

#endif /* MFEM_USE_MPI */

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

  template <bool parallel>
  std::unique_ptr<LinearSolverPreconditioner> getLinearSolverPreconditioner(
      Context& ctx,
      FiniteElementSpace<parallel>& fespace,
      const Parameters& pr) {
    if (!checkParameters(ctx, pr,
                         std::vector<std::string>{"Name", "Options"})) {
      return {};
    }
    const auto oname = get<std::string>(ctx, pr, "Name");
    if (isInvalid(oname)) {
      return {};
    }
    if (*oname == "None") {
      if (contains(pr, "Options")) {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "no options expected for preconditioner '" +
            *oname + "'");
      }
      return {};
    } else if (*oname == "AssembledOperator") {
      // preconditioner built on the assembled version of an operator applied
      // without assembly (see `MatrixFreeTangentOperator`)
      const auto oparams = get_if<Parameters>(ctx, pr, "Options", Parameters{});
      if (isInvalid(oparams)) {
        return {};
      }
      if (!checkParameters(ctx, *oparams,
                           std::vector<std::string>{"Preconditioner",
                                                    "UpdatePeriod"})) {
        return {};
      }
      const auto oprec = get<Parameters>(ctx, *oparams, "Preconditioner");
      const auto operiod = get_if<int>(ctx, *oparams, "UpdatePeriod", 1);
      if ((isInvalid(oprec)) || (isInvalid(operiod))) {
        return {};
      }
      if (*operiod < 1) {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: invalid update period of the "
            "'AssembledOperator' preconditioner");
      }
      auto p = getLinearSolverPreconditioner<parallel>(ctx, fespace, *oprec);
      if (p == nullptr) {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: no preconditioner given to the "
            "'AssembledOperator' preconditioner");
      }
      return std::make_unique<AssembledOperatorPreconditioner>(std::move(p),
                                                               *operiod);
    } else if (*oname == "OperatorJacobiSmoother") {
      // Jacobi smoother of MFEM, whose diagonal is computed by the operator
      // (`AssembleDiagonal`): usable with an operator applied without
      // assembly (see `MatrixFreeTangentOperator`)
      const auto oparams = get_if<Parameters>(ctx, pr, "Options", Parameters{});
      if (isInvalid(oparams)) {
        return {};
      }
      if (!checkParameters(ctx, *oparams, {"Damping"})) {
        return {};
      }
      const auto odamping = get_if<double>(ctx, *oparams, "Damping", 1.);
      if (isInvalid(odamping)) {
        return {};
      }
      return std::make_unique<JacobiSmoother>(*odamping);
    } else if (*oname == "OperatorChebyshevSmoother") {
      // Chebyshev smoother of MFEM: usable with an operator applied without
      // assembly
      const auto oparams = get_if<Parameters>(ctx, pr, "Options", Parameters{});
      if (isInvalid(oparams)) {
        return {};
      }
      if (!checkParameters(
              ctx, *oparams,
              std::vector<std::string>{"Order", "NumberOfPowerIterations"})) {
        return {};
      }
      const auto oorder = get_if<int>(ctx, *oparams, "Order", 2);
      const auto oiterations =
          get_if<int>(ctx, *oparams, "NumberOfPowerIterations", 10);
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
    } else if (*oname == "AmgX") {
#ifdef MFEM_MGIS_HAVE_AMGX
      const auto oparams = get_if<Parameters>(ctx, pr, "Options", Parameters{});
      if (isInvalid(oparams)) {
        return {};
      }
      return internals::makeAmgXPreconditioner<parallel>(ctx, fespace,
                                                         *oparams);
#else  /* MFEM_MGIS_HAVE_AMGX */
      return ctx.registerErrorMessage(
          "setLinearSolverPreconditioner: mfem-mgis is not compiled with "
          "AmgX (option enable-amgx)");
#endif /* MFEM_MGIS_HAVE_AMGX */
    } else if (*oname == "HypreBoomerAMG") {
      if constexpr (parallel) {
        const auto oparams =
            get_if<Parameters>(ctx, pr, "Options", Parameters{});
        if (isInvalid(oparams)) {
          return {};
        }
        return setHypreBoomerAMGPreconditioner(ctx, fespace, *oparams);
      } else {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "the 'HypreBoomerAMG' is only available in parallel");
      }
    } else if (*oname == "HypreEuclid") {
      if constexpr (parallel) {
        const auto oparams =
            get_if<Parameters>(ctx, pr, "Options", Parameters{});
        if (isInvalid(oparams)) {
          return {};
        }
        return setHypreEuclidPreconditioner(ctx, fespace, *oparams);
      } else {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "the 'HypreEuclid' is only available in parallel");
      }
    } else if (*oname == "HypreILU") {
      if constexpr (parallel) {
        const auto oparams =
            get_if<Parameters>(ctx, pr, "Options", Parameters{});
        if (isInvalid(oparams)) {
          return {};
        }
        return setHypreILUPreconditioner(ctx, fespace, *oparams);
      } else {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "the 'HypreILU' is only available in parallel");
      }
    } else if (*oname == "HypreParaSails") {
      if constexpr (parallel) {
        const auto oparams =
            get_if<Parameters>(ctx, pr, "Options", Parameters{});
        if (isInvalid(oparams)) {
          return {};
        }
        return setHypreParaSailsPreconditioner(ctx, fespace, *oparams);
      } else {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "the 'HypreParaSails' is only available in parallel");
      }
    } else if (*oname == "HypreDiagScale") {
      if constexpr (parallel) {
        const auto oparams =
            get_if<Parameters>(ctx, pr, "Options", Parameters{});
        if (isInvalid(oparams)) {
          return {};
        }
        return setHypreDiagScalePreconditioner(ctx, fespace, *oparams);
      } else {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "the 'HypreDiagScale' is only available in parallel");
      }
    } else {
      return ctx.registerErrorMessage(
          "setLinearSolverPreconditioner: "
          "unsupported preconditioner '" +
          *oname + "'");
    }
    return {};
  }  // end of getLinearSolverPreconditioner

  template <bool parallel>
  static std::optional<std::unique_ptr<LinearSolverPreconditioner>>
  setLinearSolverParameters(Context& ctx,
                            mfem::IterativeSolver& s,
                            FiniteElementSpace<parallel>& fespace,
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
    if (contains(params, Preconditioner)) {
      const auto opr = get<Parameters>(ctx, params, Preconditioner);
      if (isInvalid(opr)) {
        return {};
      }
      auto prec = getLinearSolverPreconditioner<parallel>(ctx, fespace, *opr);
      if (isInvalid(prec)) {
        return {};
      }
      s.SetPreconditioner(*prec);
      return prec;
    }
    return std::unique_ptr<LinearSolverPreconditioner>{};
  }  // end of setLinearSolverParameters

#ifdef MFEM_USE_MPI

  /*!
   * \return a generator of a linear solver applying one cycle of BoomerAMG,
   * an approximate solver used for instance by the fixed point iterations on
   * the elastic operator (see `FixedPointSolver`). The options are the ones
   * of the `HypreBoomerAMG` preconditioner.
   */
  [[nodiscard]] static LinearSolverFactory<true>::Generator
  buildHypreBoomerAMGSolverGenerator() {
    return [](Context& ctx, FiniteElementSpace<true>& fespace,
              const Parameters& params) noexcept -> LinearSolverHandler {
      auto amg = setHypreBoomerAMGPreconditioner(ctx, fespace, params);
      if (amg == nullptr) {
        return {};
      }
      return LinearSolverHandler{std::move(amg),
                                 std::unique_ptr<LinearSolverPreconditioner>{}};
    };
  }  // end of buildHypreBoomerAMGSolverGenerator

  [[nodiscard]] static std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHyprePCGSolverGenerator() {
    return [](Context& ctx, FiniteElementSpace<true>& fespace,
              const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      const char* const SolverTolerance = "Tolerance";
      const char* const Preconditioner = "Preconditioner";
      const auto allowed_parameters = std::vector<std::string>{
          Problem::SolverVerbosityLevel, SolverTolerance,
          Problem::SolverAbsoluteTolerance,
          Problem::SolverMaximumNumberOfIterations, Preconditioner};
      auto s = std::make_unique<mfem::HyprePCG>(fespace.GetComm());
      s->iterative_mode = false;
      if (!checkParameters(ctx, params, allowed_parameters)) {
        return {};
      }
      if (contains(params, Problem::SolverVerbosityLevel)) {
        const auto olvl = get<int>(ctx, params, Problem::SolverVerbosityLevel);
        if (isInvalid(olvl)) {
          return {};
        }
        s->SetPrintLevel(*olvl);
      }
      if (contains(params, Problem::SolverAbsoluteTolerance)) {
        const auto otol =
            get<double>(ctx, params, Problem::SolverAbsoluteTolerance);
        if (isInvalid(otol)) {
          return {};
        }
        s->SetAbsTol(*otol);
      }
      if (contains(params, SolverTolerance)) {
        const auto otol = get<double>(ctx, params, SolverTolerance);
        if (isInvalid(otol)) {
          return {};
        }
        s->SetTol(*otol);
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        const auto oiterMax =
            get<int>(ctx, params, Problem::SolverMaximumNumberOfIterations);
        if (isInvalid(oiterMax)) {
          return {};
        }
        s->SetMaxIter(*oiterMax);
      }
      auto prec = std::unique_ptr<mfem::Solver>{};
      if (contains(params, Preconditioner)) {
        const auto oparams = get<Parameters>(ctx, params, Preconditioner);
        if (isInvalid(oparams)) {
          return {};
        }
        prec = getLinearSolverPreconditioner<true>(ctx, fespace, *oparams);
        if (isInvalid(prec)) {
          return {};
        }
        auto* const hypre_prec =
            dynamic_cast<mfem::HypreSolver* const>(prec.get());
        if (hypre_prec == nullptr) {
          return ctx.registerErrorMessage(
              "buildHyprePCGSolverGenerator: only Hypre preconditioners "
              "allowed with "
              "the HyprePCG linear solver");
        }
        s->SetPreconditioner(*hypre_prec);
      }
      return LinearSolverHandler{std::move(s), std::move(prec)};
    };
  }  // end of buildHyprePCGSolverGenerator

  [[nodiscard]] static std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHypreGMRESSolverGenerator() {
    return [](Context& ctx, FiniteElementSpace<true>& fespace,
              const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      const char* const SolverTolerance = "Tolerance";
      const char* const Preconditioner = "Preconditioner";
      const char* const Krylov_Dimension = "KDim";
      const auto allowed_parameters = std::vector<std::string>{
          Problem::SolverVerbosityLevel, SolverTolerance,
          Problem::SolverAbsoluteTolerance,
          Problem::SolverMaximumNumberOfIterations, Preconditioner,
          Krylov_Dimension};
      auto s = std::make_unique<mfem::HypreGMRES>(fespace.GetComm());
      s->iterative_mode = false;
      if (!checkParameters(ctx, params, allowed_parameters)) {
        return {};
      }
      if (contains(params, Problem::SolverVerbosityLevel)) {
        const auto olvl = get<int>(ctx, params, Problem::SolverVerbosityLevel);
        if (isInvalid(olvl)) {
          return {};
        }
        s->SetPrintLevel(*olvl);
      }
      if (contains(params, Problem::SolverAbsoluteTolerance)) {
        const auto otol =
            get<double>(ctx, params, Problem::SolverAbsoluteTolerance);
        if (isInvalid(otol)) {
          return {};
        }
        s->SetAbsTol(*otol);
      }
      if (contains(params, SolverTolerance)) {
        const auto otol = get<double>(ctx, params, SolverTolerance);
        if (isInvalid(otol)) {
          return {};
        }
        s->SetTol(*otol);
      }
      if (contains(params, Krylov_Dimension)) {
        const auto odim = get<int>(ctx, params, Krylov_Dimension);
        if (isInvalid(odim)) {
          return {};
        }
        if (*odim < 1) {
          return ctx.registerErrorMessage(
              "invalid Krylov dimension (must be greater than one)");
        }
        s->SetKDim(*odim);
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        const auto oiterMax =
            get<int>(ctx, params, Problem::SolverMaximumNumberOfIterations);
        if (isInvalid(oiterMax)) {
          return {};
        }
        s->SetMaxIter(*oiterMax);
      }
      auto prec = std::unique_ptr<mfem::Solver>{};
      if (contains(params, Preconditioner)) {
        const auto opr = get<Parameters>(ctx, params, Preconditioner);
        if (isInvalid(opr)) {
          return {};
        }
        prec = getLinearSolverPreconditioner<true>(ctx, fespace, *opr);
        if (isInvalid(prec)) {
          return {};
        }
        auto* const hypre_prec =
            dynamic_cast<mfem::HypreSolver* const>(prec.get());
        if (hypre_prec == nullptr) {
          return ctx.registerErrorMessage(
              "buildHypreGMRESSolverGenerator: only Hypre preconditioners "
              "allowed with "
              "the HypreGMRES linear solver");
        }
        s->SetPreconditioner(*hypre_prec);
      }
      return LinearSolverHandler{std::move(s), std::move(prec)};
    };
  }  // end of buildHypreGMRESSolverGenerator

  [[nodiscard]] static std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>& p, const Parameters&)>
  buildHypreFGMRESSolverGenerator() {
    return [](Context& ctx, FiniteElementSpace<true>& fespace,
              const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      const char* const SolverTolerance = "Tolerance";
      const char* const Preconditioner = "Preconditioner";
      const char* const Krylov_Dimension = "KDim";
      const auto allowed_parameters = std::vector<std::string>{
          Problem::SolverVerbosityLevel, SolverTolerance,
          Problem::SolverAbsoluteTolerance,
          Problem::SolverMaximumNumberOfIterations, Preconditioner,
          Krylov_Dimension};
      auto s = std::make_unique<mfem::HypreFGMRES>(fespace.GetComm());
      s->iterative_mode = false;
      if (!checkParameters(ctx, params, allowed_parameters)) {
        return {};
      }
      if (contains(params, Problem::SolverVerbosityLevel)) {
        const auto olvl = get<int>(ctx, params, Problem::SolverVerbosityLevel);
        if (isInvalid(olvl)) {
          return {};
        }
        s->SetPrintLevel(*olvl);
      }
      if (contains(params, Problem::SolverAbsoluteTolerance)) {
        const auto otol =
            get<double>(ctx, params, Problem::SolverAbsoluteTolerance);
        if (isInvalid(otol)) {
          return {};
        }
        // mfem::HypreFGMRES has no `SetAbsTol` method
        HYPRE_ParCSRFlexGMRESSetAbsoluteTol(*s, *otol);
      }
      if (contains(params, SolverTolerance)) {
        const auto otol = get<double>(ctx, params, SolverTolerance);
        if (isInvalid(otol)) {
          return {};
        }
        s->SetTol(*otol);
      }
      if (contains(params, Krylov_Dimension)) {
        const auto odim = get<int>(ctx, params, Krylov_Dimension);
        if (isInvalid(odim)) {
          return {};
        }
        s->SetKDim(*odim);
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        const auto oiterMax =
            get<int>(ctx, params, Problem::SolverMaximumNumberOfIterations);
        if (isInvalid(oiterMax)) {
          return {};
        }
        s->SetMaxIter(*oiterMax);
      }
      auto prec = std::unique_ptr<mfem::Solver>{};
      if (contains(params, Preconditioner)) {
        const auto opr = get<Parameters>(ctx, params, Preconditioner);
        if (isInvalid(opr)) {
          return {};
        }
        prec = getLinearSolverPreconditioner<true>(ctx, fespace, *opr);
        if (isInvalid(prec)) {
          return {};
        }
        auto* const hypre_prec =
            dynamic_cast<mfem::HypreSolver* const>(prec.get());
        if (hypre_prec == nullptr) {
          return ctx.registerErrorMessage(
              "buildHypreFGMRESSolverGenerator: only Hypre preconditioners "
              "allowed with "
              "the HypreFGMRES linear solver");
        }
        s->SetPreconditioner(*hypre_prec);
      }
      return LinearSolverHandler{std::move(s), std::move(prec)};
    };
  }  // end of buildHypreFGMRESSolverGenerator

#else

  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHyprePCGSolverGenerator();

  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHypreGMRESSolverGenerator();

  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHypreFGMRESSolverGenerator();

#endif

  template <bool parallel, typename LinearSolverType>
  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<parallel>&, const Parameters&)>
  buildIterativeSolverGenerator() {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      return [](Context& ctx, FiniteElementSpace<true>& fespace,
                const Parameters& params) {
        auto s = std::make_unique<LinearSolverType>(fespace.GetComm());
        auto oprec = setLinearSolverParameters<true>(ctx, *s, fespace, params);
        if (isInvalid(oprec)) {
          return LinearSolverHandler{};
        }
        return LinearSolverHandler{std::move(s), std::move(*oprec)};
      };
#else  /* MFEM_USE_MPI */
      return {};
#endif /* MFEM_USE_MPI */
    } else {
      return [](Context& ctx, FiniteElementSpace<false>& fespace,
                const Parameters& params) {
        auto s = std::make_unique<LinearSolverType>();
        auto oprec = setLinearSolverParameters<false>(ctx, *s, fespace, params);
        if (isInvalid(oprec)) {
          return LinearSolverHandler{};
        }
        return LinearSolverHandler{std::move(s), std::move(*oprec)};
      };
    }
  }  // end of buildIterativeSolverGenerator

#ifdef MFEM_USE_MUMPS

  [[nodiscard]] static std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildMUMPSSolverGenerator() {
    return [](Context& ctx, FiniteElementSpace<true>& fespace,
              const Parameters& params) -> LinearSolverHandler {
      if (!checkParameters(
              ctx, params,
              std::vector<std::string>{"Symmetric", "PositiveDefinite",
                                       "ReuseAnalysis"})) {
        return {};
      }
      const auto oreuse_analysis =
          get_if<bool>(ctx, params, "ReuseAnalysis", false);
      if (isInvalid(oreuse_analysis)) {
        return {};
      }
      auto s = std::unique_ptr<mfem::MUMPSSolver>{};
      if (*oreuse_analysis) {
        s = std::make_unique<internals::MUMPSSolverReusingAnalysis>(
            fespace.GetComm());
      } else {
        s = std::make_unique<mfem::MUMPSSolver>(fespace.GetComm());
      }
      const auto osymmetric = get_if<bool>(ctx, params, "Symmetric", false);
      if (isInvalid(osymmetric)) {
        return {};
      }
      const auto opositive_definite =
          get_if<bool>(ctx, params, "PositiveDefinite", false);
      if (isInvalid(opositive_definite)) {
        return {};
      }
      s->SetPrintLevel(1);
      if (*osymmetric) {
        if (*opositive_definite) {
          s->SetMatrixSymType(
              mfem::MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
        } else {
          s->SetMatrixSymType(mfem::MUMPSSolver::MatType::SYMMETRIC_INDEFINITE);
        }
      } else {
        s->SetMatrixSymType(mfem::MUMPSSolver::MatType::UNSYMMETRIC);
      }
      return LinearSolverHandler{std::move(s),
                                 std::unique_ptr<LinearSolverPreconditioner>{}};
    };
  }  // end of builMUMPSGenerator

#endif /* MFEM_USE_MUMPS */

  template <bool parallel>
  static void declareDefaultSolvers(LinearSolverFactory<parallel>& f) {
    auto ctx = Context{};
    auto check = [&ctx](const bool b) {
      if (!b) {
        auto eh = mgis::ContractViolationHandler{};
        eh.registerErrorMessage(ctx.getErrorMessage().c_str());
      }
    };
    check(f.add(ctx, "CGSolver",
                buildIterativeSolverGenerator<parallel, mfem::CGSolver>()));
    check(f.add(ctx, "GMRESSolver",
                buildIterativeSolverGenerator<parallel, mfem::GMRESSolver>()));
    check(f.add(ctx, "FGMRESSolver",
                buildIterativeSolverGenerator<parallel, mfem::FGMRESSolver>()));
    check(
        f.add(ctx, "BiCGSTABSolver",
              buildIterativeSolverGenerator<parallel, mfem::BiCGSTABSolver>()));
    check(f.add(ctx, "MINRESSolver",
                buildIterativeSolverGenerator<parallel, mfem::MINRESSolver>()));
    check(f.add(ctx, "SLISolver",
                buildIterativeSolverGenerator<parallel, mfem::SLISolver>()));
#ifdef MFEM_MGIS_HAVE_AMGX
    // cycles of the algebraic multigrid of AmgX (approximate solver, as
    // `HypreBoomerAMG`) and conjugate gradient of AmgX
    check(f.add(ctx, "AmgX",
                [](Context& c, FiniteElementSpace<parallel>& fespace,
                   const Parameters& params) noexcept -> LinearSolverHandler {
                  auto s = internals::makeAmgXPreconditioner<parallel>(
                      c, fespace, params);
                  if (s == nullptr) {
                    return {};
                  }
                  return LinearSolverHandler{
                      std::move(s),
                      std::unique_ptr<LinearSolverPreconditioner>{}};
                }));
    check(f.add(ctx, "AmgXPCG",
                [](Context& c, FiniteElementSpace<parallel>& fespace,
                   const Parameters& params) noexcept -> LinearSolverHandler {
                  auto s = internals::makeAmgXConjugateGradientSolver<parallel>(
                      c, fespace, params);
                  if (s == nullptr) {
                    return {};
                  }
                  return LinearSolverHandler{
                      std::move(s),
                      std::unique_ptr<LinearSolverPreconditioner>{}};
                }));
#endif /* MFEM_MGIS_HAVE_AMGX */
    if constexpr (parallel) {
#ifdef MFEM_USE_MUMPS
      check(f.add(ctx, "MUMPSSolver", buildMUMPSSolverGenerator()));
#endif
      check(f.add(ctx, "HyprePCG", buildHyprePCGSolverGenerator()));
      check(f.add(ctx, "HypreBoomerAMG", buildHypreBoomerAMGSolverGenerator()));
      check(f.add(ctx, "HypreGMRES", buildHypreGMRESSolverGenerator()));
      check(f.add(ctx, "HypreFGMRES", buildHypreFGMRESSolverGenerator()));
    }
    if constexpr (!parallel) {
#ifdef MFEM_USE_SUITESPARSE
      check(f.add(ctx, "UMFPackSolver",
                  [](Context&, FiniteElementSpace<false>&, const Parameters&) {
                    return LinearSolverHandler{
                        std::make_unique<mfem::UMFPackSolver>(),
                        std::unique_ptr<LinearSolverPreconditioner>{}};
                  }));
#endif
    }
  }  // end of declareDefaultSolvers

#ifdef MFEM_USE_MPI

  LinearSolverFactory<true>& LinearSolverFactory<true>::getFactory() {
    static LinearSolverFactory<true> factory;
    return factory;
  }  // end of getFactory

  bool LinearSolverFactory<true>::add(Context& ctx,
                                      std::string_view n,
                                      Generator g) noexcept {
    const auto pg = this->generators.find(n);
    if (pg != this->generators.end()) {
      auto msg = std::string{"LinearSolverFactory<true>::add: "};
      msg += "a linear solver called '";
      msg += n;
      msg += "' has already been declared";
      return ctx.registerErrorMessage(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
    return true;
  }  // end of add

  LinearSolverHandler LinearSolverFactory<true>::generate(
      Context& ctx,
      std::string_view n,
      FiniteElementSpace<true>& fespace,
      const Parameters& params) const {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      std::string msg("LinearSolverFactory<true>::generate: ");
      msg += "no linear solver called '";
      msg += n;
      msg += "' declared. Here is the list of available linear solvers:";
      for (const auto& g : this->generators) {
        msg += " " + g.first;
      }
      return ctx.registerErrorMessage(msg);
    }
    const auto& g = pg->second;
    LinearSolverHandler s;
    try {
      s = g(ctx, fespace, params);
    } catch (std::exception& e) {
      std::string msg("LinearSolverFactory<true>::generate: ");
      msg += "error while generating linear solver '";
      msg += n;
      msg += "'\n";
      msg += e.what();
      return ctx.registerErrorMessage(msg);
    }
    return s;
  }  // end of generate

  LinearSolverFactory<true>::LinearSolverFactory() {
    declareDefaultSolvers(*this);
  }  // end of LinearSolverFactory

  LinearSolverFactory<true>::~LinearSolverFactory() = default;

#endif /* MFEM_USE_MPI */

  LinearSolverFactory<false>& LinearSolverFactory<false>::getFactory() {
    static LinearSolverFactory<false> factory;
    return factory;
  }  // end of getFactory

  bool LinearSolverFactory<false>::add(Context& ctx,
                                       std::string_view n,
                                       Generator g) noexcept {
    const auto pg = this->generators.find(n);
    if (pg != this->generators.end()) {
      auto msg = std::string{"LinearSolverFactory<false>::add: "};
      msg += "a linear solver called '";
      msg += n;
      msg += "' has already been declared";
      return ctx.registerErrorMessage(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
    return true;
  }  // end of add

  LinearSolverHandler LinearSolverFactory<false>::generate(
      Context& ctx,
      std::string_view n,
      FiniteElementSpace<false>& fespace,
      const Parameters& params) const {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      std::string msg("LinearSolverFactory<false>::generate: ");
      msg += "no linear solver called '";
      msg += n;
      msg += "' declared. Here is the list of available linear solvers:";
      for (const auto& g : this->generators) {
        msg += " " + g.first;
      }
      return ctx.registerErrorMessage(msg);
    }
    const auto& g = pg->second;
    LinearSolverHandler s;
    try {
      s = g(ctx, fespace, params);
    } catch (std::exception& e) {
      std::string msg("LinearSolverFactory<false>::generate: ");
      msg += "error while generating linear solver '";
      msg += n;
      msg += "'\n";
      msg += e.what();
      return ctx.registerErrorMessage(msg);
    }
    return s;
  }  // end of generate

  LinearSolverFactory<false>::LinearSolverFactory() {
    declareDefaultSolvers(*this);
  }  // end of LinearSolverFactory

  LinearSolverFactory<false>::~LinearSolverFactory() = default;

}  // end of namespace mfem_mgis
