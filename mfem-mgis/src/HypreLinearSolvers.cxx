/*!
 * \file   src/HypreLinearSolvers.cxx
 * \brief  This file implements the module of the linear solvers and of the
 * preconditioners of hypre (see `LinearSolverModules.hxx`)
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <array>
#include <tuple>
#include <memory>
#include <string>
#include <vector>
#include <optional>
#include <utility>
#include <algorithm>
#include <functional>
#include <string_view>
#include "mfem/config/config.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "ExecutionSettings.hxx"
#include "LinearSolverModules.hxx"

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
   * operator, whose memory may then be freed (invalid communicator). The
   * operator is replaced by the new one before calling this method, which
   * reports an operator which is not a matrix of hypre as before.
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

namespace mfem_mgis {

  /*!
   * \return the options of the BoomerAMG preconditioner of hypre whose
   * values are integers (see `setHypreBoomerAMGPreconditioner`)
   */
  [[nodiscard]] static std::vector<std::string>
  getHypreBoomerAMGIntegerOptions() {
    return {"CoarseningType", "AggressiveCoarseningLevels", "InterpolationType",
            "RelaxationType", "RelaxationSweeps",           "CycleType",
            "NumberOfCycles", "MaximumNumberOfLevels"};
  }  // end of getHypreBoomerAMGIntegerOptions

  //! \return the options of the BoomerAMG preconditioner of hypre
  [[nodiscard]] static std::vector<std::string> getHypreBoomerAMGOptions() {
    using Problem = AbstractNonLinearEvolutionProblem;
    auto options = getHypreBoomerAMGIntegerOptions();
    options.insert(options.end(), {"Strategy", "StrengthThreshold",
                                   Problem::SolverVerbosityLevel});
    return options;
  }  // end of getHypreBoomerAMGOptions

  /*!
   * \return the BoomerAMG preconditioner of hypre
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options:
   * - `Strategy`: `System` (systems options, by default for vector unknowns),
   *   `Elasticity` (elasticity options, requiring unknowns ordered by vector
   *   dimension and hypre computing on the host, `System` being used
   *   otherwise) or `None` (by default for scalar unknowns);
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
    const auto integer_options = getHypreBoomerAMGIntegerOptions();
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
    // 3.2.0 and CUDA 13), the ILUT factorisations do not: the level of fill
    // is then zero by default, one being the default of MFEM
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
  }  // end of setHypreILUPreconditioner

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
  }  // end of setHypreDiagScalePreconditioner

  /*!
   * \brief set the preconditioner of a Krylov solver of hypre, which only
   * accepts the preconditioners of hypre (see
   * `internals::makePreconditioner`)
   * \return the preconditioner (a null pointer if none is requested), or an
   * empty value on error
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] n: name of the linear solver
   * \param[in] d: description of the linear solver
   * \param[in, out] s: linear solver
   * \param[in] params: parameters of the linear solver
   */
  template <typename SolverType>
  [[nodiscard]] static std::optional<
      std::unique_ptr<LinearSolverPreconditioner>>
  setHypreKrylovSolverPreconditioner(Context& ctx,
                                     FiniteElementSpace<true>& fespace,
                                     const std::string_view n,
                                     const LinearSolverDescription& d,
                                     SolverType& s,
                                     const Parameters& params) {
    auto oprec =
        internals::makePreconditioner<true>(ctx, fespace, n, d, params);
    if (isInvalid(oprec)) {
      return {};
    }
    if (*oprec != nullptr) {
      auto* const hypre_prec =
          dynamic_cast<mfem::HypreSolver* const>(oprec->get());
      if (hypre_prec == nullptr) {
        std::ignore = ctx.registerErrorMessage(
            "the linear solver '" + std::string{n} +
            "' only accepts the preconditioners of hypre");
        return {};
      }
      s.SetPreconditioner(*hypre_prec);
    }
    return std::move(*oprec);
  }  // end of setHypreKrylovSolverPreconditioner

  /*!
   * \return a generator of a linear solver applying cycles of BoomerAMG
   * (one by default, see the option `NumberOfCycles`), an approximate
   * solver used for instance by the fixed point iterations on the elastic
   * operator (see `FixedPointSolver`). The options are the ones of the
   * `HypreBoomerAMG` preconditioner.
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
  buildHyprePCGSolverGenerator(const LinearSolverDescription d) {
    return [d](Context& ctx, FiniteElementSpace<true>& fespace,
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
        s->SetTol(boundRelativeTolerance(*otol));
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        const auto oiterMax =
            get<int>(ctx, params, Problem::SolverMaximumNumberOfIterations);
        if (isInvalid(oiterMax)) {
          return {};
        }
        s->SetMaxIter(*oiterMax);
      }
      auto oprec = setHypreKrylovSolverPreconditioner(ctx, fespace, "HyprePCG",
                                                      d, *s, params);
      if (isInvalid(oprec)) {
        return {};
      }
      return LinearSolverHandler{std::move(s), std::move(*oprec)};
    };
  }  // end of buildHyprePCGSolverGenerator

  [[nodiscard]] static std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHypreGMRESSolverGenerator(const LinearSolverDescription d) {
    return [d](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      const char* const SolverTolerance = "Tolerance";
      const char* const Preconditioner = "Preconditioner";
      const char* const Krylov_Dimension = "KDim";
      const auto allowed_parameters =
          std::vector<std::string>{Problem::SolverVerbosityLevel,
                                   SolverTolerance,
                                   Problem::SolverAbsoluteTolerance,
                                   Problem::SolverMaximumNumberOfIterations,
                                   Preconditioner,
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
        s->SetTol(boundRelativeTolerance(*otol));
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
      auto oprec = setHypreKrylovSolverPreconditioner(
          ctx, fespace, "HypreGMRES", d, *s, params);
      if (isInvalid(oprec)) {
        return {};
      }
      return LinearSolverHandler{std::move(s), std::move(*oprec)};
    };
  }  // end of buildHypreGMRESSolverGenerator

  [[nodiscard]] static std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>& p, const Parameters&)>
  buildHypreFGMRESSolverGenerator(const LinearSolverDescription d) {
    return [d](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      const char* const SolverTolerance = "Tolerance";
      const char* const Preconditioner = "Preconditioner";
      const char* const Krylov_Dimension = "KDim";
      const auto allowed_parameters =
          std::vector<std::string>{Problem::SolverVerbosityLevel,
                                   SolverTolerance,
                                   Problem::SolverAbsoluteTolerance,
                                   Problem::SolverMaximumNumberOfIterations,
                                   Preconditioner,
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
        s->SetTol(boundRelativeTolerance(*otol));
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
      auto oprec = setHypreKrylovSolverPreconditioner(
          ctx, fespace, "HypreFGMRES", d, *s, params);
      if (isInvalid(oprec)) {
        return {};
      }
      return LinearSolverHandler{std::move(s), std::move(*oprec)};
    };
  }  // end of buildHypreFGMRESSolverGenerator
}  // end of namespace mfem_mgis

namespace mfem_mgis::internals {

  /*!
   * \return the options of a Krylov solver of hypre
   * \param[in] gmres: if true, the solver is a GMRES solver, whose
   * dimension of the Krylov space can be given
   */
  [[nodiscard]] static std::vector<std::string> getHypreKrylovSolverOptions(
      const bool gmres) {
    using Problem = AbstractNonLinearEvolutionProblem;
    auto options = std::vector<std::string>{
        Problem::SolverVerbosityLevel, "Tolerance",
        Problem::SolverAbsoluteTolerance,
        Problem::SolverMaximumNumberOfIterations, "Preconditioner"};
    if (gmres) {
      options.push_back("KDim");
    }
    return options;
  }  // end of getHypreKrylovSolverOptions

  void declareHypreLinearSolvers(LinearSolverFactory<true>& f) {
    const auto krylov = [](const bool gmres, std::string description) {
      return LinearSolverDescription{
          .library = "hypre",
          .kind = LinearSolverKind::ITERATIVE,
          .placement = LinearSolverPlacement::HYPRE,
          .preconditioners = AcceptedPreconditioners::SAME_LIBRARY,
          .options = getHypreKrylovSolverOptions(gmres),
          .description = std::move(description)};
    };
    const auto pcg = krylov(false, "conjugate gradient of hypre");
    declare(f, "HyprePCG", pcg, buildHyprePCGSolverGenerator(pcg));
    const auto gmres = krylov(true, "GMRES of hypre");
    declare(f, "HypreGMRES", gmres, buildHypreGMRESSolverGenerator(gmres));
    const auto fgmres = krylov(true, "flexible GMRES of hypre");
    declare(f, "HypreFGMRES", fgmres, buildHypreFGMRESSolverGenerator(fgmres));
    declare(f, "HypreBoomerAMG",
            LinearSolverDescription{
                .library = "hypre",
                .kind = LinearSolverKind::APPROXIMATE,
                .placement = LinearSolverPlacement::HYPRE,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = getHypreBoomerAMGOptions(),
                .description = "cycles of the BoomerAMG algebraic multigrid "
                               "of hypre (one by default)"},
            buildHypreBoomerAMGSolverGenerator());
  }  // end of declareHypreLinearSolvers

  void declareHyprePreconditioners(PreconditionerFactory<true>& f) {
    using Problem = AbstractNonLinearEvolutionProblem;
    const auto preconditioner = [](const LinearSolverPlacement p,
                                   std::vector<std::string> options,
                                   std::string description) {
      return LinearSolverDescription{
          .library = "hypre",
          .kind = LinearSolverKind::PRECONDITIONER,
          .placement = p,
          .preconditioners = AcceptedPreconditioners::NONE,
          .options = std::move(options),
          .description = std::move(description)};
    };
    declare(
        f, "HypreBoomerAMG",
        preconditioner(LinearSolverPlacement::HYPRE, getHypreBoomerAMGOptions(),
                       "BoomerAMG algebraic multigrid of hypre"),
        setHypreBoomerAMGPreconditioner);
    // Euclid and ParaSails only compute on the host
    declare(f, "HypreEuclid",
            preconditioner(LinearSolverPlacement::HOST,
                           {Problem::SolverVerbosityLevel},
                           "Euclid incomplete factorization of hypre"),
            setHypreEuclidPreconditioner);
    declare(f, "HypreILU",
            preconditioner(LinearSolverPlacement::HYPRE,
                           {Problem::SolverVerbosityLevel, "HypreILUType",
                            "HypreILULevelOfFill"},
                           "incomplete factorizations of hypre"),
            setHypreILUPreconditioner);
    declare(f, "HypreParaSails",
            preconditioner(LinearSolverPlacement::HOST,
                           {Problem::SolverVerbosityLevel},
                           "ParaSails sparse approximate inverse of hypre"),
            setHypreParaSailsPreconditioner);
    declare(f, "HypreDiagScale",
            preconditioner(LinearSolverPlacement::HYPRE,
                           {Problem::SolverVerbosityLevel},
                           "diagonal scaling of hypre"),
            setHypreDiagScalePreconditioner);
  }  // end of declareHyprePreconditioners

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_USE_MPI */
