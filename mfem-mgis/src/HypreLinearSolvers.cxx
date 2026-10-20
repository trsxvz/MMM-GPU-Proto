/*!
 * \file   src/HypreLinearSolvers.cxx
 * \brief  This file implements the module of the linear solvers and of the
 * preconditioners of hypre (see `LinearSolverModules.hxx`): the ones wrapped
 * by MFEM and the ones of `HypreSolvers.hxx`
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <cmath>
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
#include <type_traits>
#include "mfem/config/config.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "ExecutionSettings.hxx"
#include "HypreSolvers.hxx"
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

  //! \brief an option of BoomerAMG given to hypre
  struct HypreBoomerAMGOption {
    //! \brief name of the option
    const char* name;
    //! \brief minimum value
    double minimum;
    //! \brief maximum value
    double maximum;
    //! \brief if true, the value is an integer
    bool integer;
    //! \brief setting of the option
    std::function<void(mfem::HypreBoomerAMG&, double)> set;
  };

  /*!
   * \return the options of the BoomerAMG preconditioner of hypre given to
   * hypre (see `setHypreBoomerAMGPreconditioner`), whose values are numbers
   */
  [[nodiscard]] static const std::vector<HypreBoomerAMGOption>&
  getHypreBoomerAMGNumericalOptions() {
    using AMG = mfem::HypreBoomerAMG;
    constexpr auto large = 1e30;
    static const auto options = std::vector<HypreBoomerAMGOption>{
        {"StrongThreshold", 0, 1, false,
         [](AMG& amg, const double v) { amg.SetStrengthThresh(v); }},
        {"CoarseningType", 0, large, true,
         [](AMG& amg, const double v) {
           amg.SetCoarsening(static_cast<int>(v));
         }},
        {"AggressiveCoarseningLevels", 0, large, true,
         [](AMG& amg, const double v) {
           amg.SetAggressiveCoarsening(static_cast<int>(v));
         }},
        // types 0 to 25 of hypre (the type 100 is not supported)
        {"InterpolationType", 0, 25, true,
         [](AMG& amg, const double v) {
           amg.SetInterpolation(static_cast<int>(v));
         }},
        {"RelaxationType", 0, large, true,
         [](AMG& amg, const double v) {
           amg.SetRelaxType(static_cast<int>(v));
         }},
        {"RelaxationSweeps", 1, large, true,
         [](AMG& amg, const double v) {
           amg.SetCycleNumSweeps(static_cast<int>(v), static_cast<int>(v));
         }},
        {"CycleType", 1, 2, true,
         [](AMG& amg, const double v) {
           amg.SetCycleType(static_cast<int>(v));
         }},
        {"NumberOfCycles", 1, large, true,
         [](AMG& amg, const double v) { amg.SetMaxIter(static_cast<int>(v)); }},
        {"MaximumNumberOfLevels", 1, large, true,
         [](AMG& amg, const double v) {
           amg.SetMaxLevels(static_cast<int>(v));
         }},
        // options which MFEM does not wrap
        // strictly lower than one for hypre
        {"TruncationFactor", 0, std::nextafter(1., 0.), false,
         [](AMG& amg, const double v) {
           HYPRE_BoomerAMGSetTruncFactor(amg, v);
         }},
        {"MaximumElementsPerRow", 0, large, true,
         [](AMG& amg, const double v) {
           HYPRE_BoomerAMGSetPMaxElmts(amg, static_cast<int>(v));
         }},
        {"AggressiveInterpolationType", 0, 9, true,
         [](AMG& amg, const double v) {
           HYPRE_BoomerAMGSetAggInterpType(amg, static_cast<int>(v));
         }},
        {"NumberOfPaths", 1, large, true,
         [](AMG& amg, const double v) {
           HYPRE_BoomerAMGSetNumPaths(amg, static_cast<int>(v));
         }},
        {"MaximumCoarseSize", 1, large, true,
         [](AMG& amg, const double v) {
           HYPRE_BoomerAMGSetMaxCoarseSize(amg, static_cast<int>(v));
         }},
        {"MinimumCoarseSize", 0, large, true,
         [](AMG& amg, const double v) {
           HYPRE_BoomerAMGSetMinCoarseSize(amg, static_cast<int>(v));
         }},
        {"RelaxationWeight", 0, large, false,
         [](AMG& amg, const double v) { HYPRE_BoomerAMGSetRelaxWt(amg, v); }},
        {"ChebyshevOrder", 1, 4, true, [](AMG& amg, const double v) {
           HYPRE_BoomerAMGSetChebyOrder(amg, static_cast<int>(v));
         }}};
    return options;
  }  // end of getHypreBoomerAMGNumericalOptions

  //! \return the options of the BoomerAMG preconditioner of hypre
  [[nodiscard]] static std::vector<std::string> getHypreBoomerAMGOptions() {
    using Problem = AbstractNonLinearEvolutionProblem;
    auto options = std::vector<std::string>{"Strategy", "KeepTranspose",
                                            Problem::SolverVerbosityLevel};
    for (const auto& o : getHypreBoomerAMGNumericalOptions()) {
      options.push_back(o.name);
    }
    return options;
  }  // end of getHypreBoomerAMGOptions

  /*!
   * \return the BoomerAMG preconditioner of hypre
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options:
   * - `Strategy`: `System` (systems options, by default for vector unknowns),
   *   `Elasticity` (elasticity options, for unknowns with as many components
   *   as the dimension of the mesh, 2 or 3, which is the space dimension,
   *   requiring unknowns ordered by vector dimension and hypre computing on
   *   the host, `System` being used otherwise) or `None` (by default for
   *   scalar unknowns);
   * - `VerbosityLevel`: print level of hypre;
   * - options of hypre, whose values are the ones of hypre (see its
   *   documentation): `StrongThreshold`, `CoarseningType`,
   *   `AggressiveCoarseningLevels`, `InterpolationType`, `RelaxationType`,
   *   `RelaxationSweeps` (sweeps before and after the coarse levels),
   *   `CycleType` (1: V-cycle, 2: W-cycle), `NumberOfCycles` (cycles per
   *   application, 1 by default), `MaximumNumberOfLevels`,
   *   `TruncationFactor` and `MaximumElementsPerRow` (truncation of the
   *   interpolation), `AggressiveInterpolationType` and `NumberOfPaths`
   *   (aggressive coarsening), `MaximumCoarseSize` and `MinimumCoarseSize`
   *   (size of the coarsest level), `RelaxationWeight`, `ChebyshevOrder`
   *   (Chebyshev smoothers, `RelaxationType` 16) and `KeepTranspose`
   *   (transposed interpolations kept, faster on a GPU).
   */
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  setHypreBoomerAMGPreconditioner(Context& ctx,
                                  FiniteElementSpace<true>& fespace,
                                  const Parameters& opts) noexcept {
    using Problem = AbstractNonLinearEvolutionProblem;
    using Setting = internals::ConfiguredHypreBoomerAMG::Setting;
    if (!checkParameters(ctx, opts, getHypreBoomerAMGOptions())) {
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
      // MFEM computes the rigid body modes with as many components as the
      // dimension of the mesh, from the first coordinates of the nodes: they
      // are only defined for unknowns with as many components as this
      // dimension, which shall be the space dimension (MFEM aborts or
      // computes wrong modes otherwise)
      const auto* const mesh = fespace.GetParMesh();
      const auto d = mesh->Dimension();
      if ((fespace.GetVDim() != d) || (mesh->SpaceDimension() != d) ||
          ((d != 2) && (d != 3))) {
        return ctx.registerErrorMessage(
            "setLinearSolverParameters: strategy 'Elasticity' of "
            "preconditioner HypreBoomerAMG requires unknowns with as many "
            "components as the dimension of the mesh (2 or 3), which "
            "shall be the space dimension");
      }
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
    for (const auto& o : getHypreBoomerAMGNumericalOptions()) {
      if (!contains(opts, o.name)) {
        continue;
      }
      auto v = 0.;
      if (o.integer) {
        const auto ov = get<int>(ctx, opts, o.name);
        if (isInvalid(ov)) {
          return {};
        }
        v = *ov;
      } else {
        const auto ov = get<double>(ctx, opts, o.name);
        if (isInvalid(ov)) {
          return {};
        }
        v = *ov;
      }
      if ((v < o.minimum) || (v > o.maximum)) {
        return ctx.registerErrorMessage(
            "setLinearSolverParameters: invalid "
            "value for option '" +
            std::string{o.name} + "' of HypreBoomerAMG");
      }
      settings.push_back(
          [v, set = o.set](mfem::HypreBoomerAMG& amg) { set(amg, v); });
    }
    if (contains(opts, "KeepTranspose")) {
      const auto ov = get<bool>(ctx, opts, "KeepTranspose");
      if (isInvalid(ov)) {
        return {};
      }
      settings.push_back([v = *ov](mfem::HypreBoomerAMG& amg) {
        HYPRE_BoomerAMGSetKeepTranspose(amg, v ? 1 : 0);
      });
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
   * \return the smoother of hypre (`mfem::HypreSmoother`)
   * \param[in, out] ctx: execution context
   * \param[in] opts: options:
   * - `Type`: type of the smoother (`Jacobi`, `l1Jacobi`, `l1GS`, `l1GStr`,
   *   `lumpedJacobi`, `GS`, `OPFS`, `Chebyshev`, `Taubin` or `FIR`), `l1GS`
   *   when hypre computes on the host, `l1Jacobi` on a GPU, by default;
   * - `RelaxationSweeps`: number of sweeps (one by default);
   * - `RelaxationWeight` and `Omega`: damping and SOR parameters (one by
   *   default);
   * - `PolynomialOrder`, `PolynomialFraction` and
   *   `EigenvalueEstimationIterations`: options of the polynomial smoothers
   *   (2, 0.3 and 10 by default; with no iteration, the largest eigenvalue
   *   is estimated by hypre);
   * - `TaubinLambda`, `TaubinMu` and `TaubinIterations`: options of the
   *   Taubin smoother (0.5, -0.5 and 40 by default);
   * - `Window`: window of the FIR smoother (`Rectangular`, by default,
   *   `Hanning`, `Hamming` or `Blackman`).
   */
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  setHypreSmootherPreconditioner(Context& ctx,
                                 FiniteElementSpace<true>&,
                                 const Parameters& opts) {
    using Smoother = mfem::HypreSmoother;
    if (!checkParameters(ctx, opts,
                         std::vector<std::string>{
                             "Type", "RelaxationSweeps", "RelaxationWeight",
                             "Omega", "PolynomialOrder", "PolynomialFraction",
                             "EigenvalueEstimationIterations", "TaubinLambda",
                             "TaubinMu", "TaubinIterations", "Window"})) {
      return {};
    }
    const auto types = std::array<std::pair<const char*, Smoother::Type>, 10>{
        {{"Jacobi", Smoother::Jacobi},
         {"l1Jacobi", Smoother::l1Jacobi},
         {"l1GS", Smoother::l1GS},
         {"l1GStr", Smoother::l1GStr},
         {"lumpedJacobi", Smoother::lumpedJacobi},
         {"GS", Smoother::GS},
         {"OPFS", Smoother::OPFS},
         {"Chebyshev", Smoother::Chebyshev},
         {"Taubin", Smoother::Taubin},
         {"FIR", Smoother::FIR}}};
    auto type = Smoother::DefaultType();
    if (contains(opts, "Type")) {
      const auto otype = get<std::string>(ctx, opts, "Type");
      if (isInvalid(otype)) {
        return {};
      }
      const auto p =
          std::find_if(types.begin(), types.end(),
                       [&otype](const auto& t) { return *otype == t.first; });
      if (p == types.end()) {
        return ctx.registerErrorMessage(
            "setHypreSmootherPreconditioner: invalid type '" + *otype +
            "' of the smoother of hypre");
      }
      type = p->second;
    }
    // the on-processor forward solve sorts the matrix on the host
    // (segmentation fault in hypre_topo_sort with hypre 3.2.0 on a GPU), MFEM
    // applies the Taubin and FIR smoothers by loops on the host and computes
    // the sums of the rows of the lumped Jacobi smoother in a vector of the
    // host wrapping the memory of hypre
    if (((type == Smoother::OPFS) || (type == Smoother::Taubin) ||
         (type == Smoother::FIR) || (type == Smoother::lumpedJacobi)) &&
        (isHypreUsingGPU())) {
      return ctx.registerErrorMessage(
          "setHypreSmootherPreconditioner: the smoothers 'OPFS', 'Taubin', "
          "'FIR' and 'lumpedJacobi' only compute on the host, which hypre "
          "does not");
    }
    const auto osweeps = get_if<int>(ctx, opts, "RelaxationSweeps", 1);
    const auto oweight = get_if<double>(ctx, opts, "RelaxationWeight", 1.);
    const auto oomega = get_if<double>(ctx, opts, "Omega", 1.);
    const auto oorder = get_if<int>(ctx, opts, "PolynomialOrder", 2);
    const auto ofraction = get_if<double>(ctx, opts, "PolynomialFraction", 0.3);
    const auto oeig =
        get_if<int>(ctx, opts, "EigenvalueEstimationIterations", 10);
    const auto olambda = get_if<double>(ctx, opts, "TaubinLambda", 0.5);
    const auto omu = get_if<double>(ctx, opts, "TaubinMu", -0.5);
    const auto otaubin = get_if<int>(ctx, opts, "TaubinIterations", 40);
    if (isInvalid(osweeps) || isInvalid(oweight) || isInvalid(oomega) ||
        isInvalid(oorder) || isInvalid(ofraction) || isInvalid(oeig) ||
        isInvalid(olambda) || isInvalid(omu) || isInvalid(otaubin)) {
      return {};
    }
    if ((*osweeps < 1) || (*oorder < 1) || (*oeig < 0) || (*otaubin < 1)) {
      return ctx.registerErrorMessage(
          "setHypreSmootherPreconditioner: invalid value of the number of "
          "sweeps, of the order of the polynomial or of the number of "
          "iterations");
    }
    auto smoother = std::make_unique<Smoother>();
    smoother->SetType(type, *osweeps);
    smoother->SetSOROptions(*oweight, *oomega);
    smoother->SetPolyOptions(*oorder, *ofraction, *oeig);
    smoother->SetTaubinOptions(*olambda, *omu, *otaubin);
    // the parameters of the window are not initialized by the default
    // constructor of `mfem::HypreSmoother`
    const auto owindow =
        get_if<std::string>(ctx, opts, "Window", std::string{"Rectangular"});
    if (isInvalid(owindow)) {
      return {};
    }
    if ((*owindow != "Rectangular") && (*owindow != "Hanning") &&
        (*owindow != "Hamming") && (*owindow != "Blackman")) {
      return ctx.registerErrorMessage(
          "setHypreSmootherPreconditioner: invalid window '" + *owindow +
          "' of the FIR smoother");
    }
    smoother->SetWindowByName(owindow->c_str());
    return smoother;
  }  // end of setHypreSmootherPreconditioner

  /*!
   * \brief an option of a solver of hypre whose value is given to hypre
   * \tparam ValueType: type of the value (`HYPRE_Int` or `HYPRE_Real`)
   */
  template <typename ValueType>
  struct HypreOption {
    //! \brief name of the option
    const char* name;
    //! \brief minimum value
    ValueType minimum;
    //! \brief maximum value
    ValueType maximum;
    //! \brief function of hypre setting the option
    HYPRE_Int (*set)(HYPRE_Solver, ValueType);
  };

  //! \brief largest value of an option of hypre
  static constexpr auto hypre_option_maximum = HYPRE_Int{1 << 30};

  /*!
   * \brief add the settings of the options of a solver of hypre given by the
   * parameters, refusing the values out of the ranges accepted by hypre
   * (which ignores them)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] settings: settings
   * \param[in] params: parameters
   * \param[in] options: options
   * \param[in] n: name of the solver, for the error messages
   */
  template <typename ValueType, std::size_t N>
  [[nodiscard]] static bool addHypreSettings(
      Context& ctx,
      std::vector<internals::HypreSetting>& settings,
      const Parameters& params,
      const std::array<HypreOption<ValueType>, N>& options,
      const std::string_view n) {
    for (const auto& o : options) {
      if (!contains(params, o.name)) {
        continue;
      }
      auto v = ValueType{};
      if constexpr (std::is_same_v<ValueType, HYPRE_Real>) {
        const auto ov = get<double>(ctx, params, o.name);
        if (isInvalid(ov)) {
          return false;
        }
        v = static_cast<ValueType>(*ov);
      } else {
        const auto ov = get<int>(ctx, params, o.name);
        if (isInvalid(ov)) {
          return false;
        }
        v = static_cast<ValueType>(*ov);
      }
      if ((v < o.minimum) || (v > o.maximum)) {
        return ctx.registerErrorMessage("invalid value for option '" +
                                        std::string{o.name} + "' of " +
                                        std::string{n});
      }
      settings.push_back([f = o.set, v](HYPRE_Solver s) { f(s, v); });
    }
    return true;
  }  // end of addHypreSettings

  /*!
   * \return the factorized sparse approximate inverse of hypre (FSAI)
   * \param[in, out] ctx: execution context
   * \param[in] opts: options, whose values are the ones of hypre (see its
   * documentation): `AlgorithmType`, `LocalSolveType`, `MaximumSteps`,
   * `MaximumStepSize`, `MaximumNonZerosPerRow`, `NumberOfLevels`,
   * `Threshold`, `KapTolerance`, `EigenvalueEstimationIterations` and
   * `VerbosityLevel`
   */
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  setHypreFSAIPreconditioner(Context& ctx,
                             FiniteElementSpace<true>&,
                             const Parameters& opts) {
    using Problem = AbstractNonLinearEvolutionProblem;
    constexpr auto m = hypre_option_maximum;
    const auto integers = std::array<HypreOption<HYPRE_Int>, 8>{
        {{"AlgorithmType", 1, 3, HYPRE_FSAISetAlgoType},
         {"LocalSolveType", 0, 2, HYPRE_FSAISetLocalSolveType},
         {"MaximumSteps", 0, m, HYPRE_FSAISetMaxSteps},
         {"MaximumStepSize", 0, m, HYPRE_FSAISetMaxStepSize},
         {"MaximumNonZerosPerRow", 0, m, HYPRE_FSAISetMaxNnzRow},
         {"NumberOfLevels", 1, m, HYPRE_FSAISetNumLevels},
         {"EigenvalueEstimationIterations", 0, m, HYPRE_FSAISetEigMaxIters},
         {Problem::SolverVerbosityLevel, 0, m, HYPRE_FSAISetPrintLevel}}};
    const auto reals = std::array<HypreOption<HYPRE_Real>, 2>{
        {{"Threshold", 0, HYPRE_Real{m}, HYPRE_FSAISetThreshold},
         {"KapTolerance", 0, HYPRE_Real{m}, HYPRE_FSAISetKapTolerance}}};
    auto allowed = std::vector<std::string>{};
    for (const auto& o : integers) {
      allowed.push_back(o.name);
    }
    for (const auto& o : reals) {
      allowed.push_back(o.name);
    }
    if (!checkParameters(ctx, opts, allowed)) {
      return {};
    }
    auto settings = std::vector<internals::HypreSetting>{};
    const auto n = "the factorized sparse approximate inverse of hypre";
    if ((!addHypreSettings(ctx, settings, opts, integers, n)) ||
        (!addHypreSettings(ctx, settings, opts, reals, n))) {
      return {};
    }
    return std::make_unique<internals::HypreFSAI>(std::move(settings));
  }  // end of setHypreFSAIPreconditioner

  /*!
   * \return the components kept at each level of the multigrid reduction of
   * hypre, given by the option `CoarseComponents` (a list of components, for
   * one level, or a list of lists of components, one per level, each level
   * keeping a subset of the components of the previous one), or an empty
   * value on error. By default, one level keeps all the components but the
   * last one.
   * \param[in, out] ctx: execution context
   * \param[in] opts: options
   * \param[in] vdim: number of components of the unknowns
   */
  [[nodiscard]] static std::optional<std::vector<std::vector<HYPRE_Int>>>
  getMGRCoarseComponents(Context& ctx, const Parameters& opts, const int vdim) {
    const auto n = "the multigrid reduction of hypre";
    auto levels = std::vector<std::vector<HYPRE_Int>>{};
    if (!contains(opts, "CoarseComponents")) {
      if (vdim < 2) {
        return ctx.registerErrorMessage(
            std::string{n} +
            " requires at least two components of the "
            "unknowns");
      }
      auto l = std::vector<HYPRE_Int>{};
      for (int c = 0; c + 1 < vdim; ++c) {
        l.push_back(c);
      }
      levels.push_back(std::move(l));
      return levels;
    }
    const auto ovalues =
        get<std::vector<Parameter>>(ctx, opts, "CoarseComponents");
    if (isInvalid(ovalues)) {
      return {};
    }
    // a list of integers gives one level, a list of lists one level per list
    const auto read = [&ctx, vdim, n](const std::vector<Parameter>& values)
        -> std::optional<std::vector<HYPRE_Int>> {
      auto l = std::vector<HYPRE_Int>{};
      for (const auto& v : values) {
        const auto oc = get<int>(ctx, v);
        if (isInvalid(oc)) {
          return {};
        }
        const auto c = *oc;
        if ((c < 0) || (c >= vdim) ||
            (std::find(l.begin(), l.end(), c) != l.end())) {
          return ctx.registerErrorMessage(
              std::string{"invalid component in the option "
                          "'CoarseComponents' of "} +
              n);
        }
        l.push_back(static_cast<HYPRE_Int>(c));
      }
      return l;
    };
    if ((!ovalues->empty()) && is<std::vector<Parameter>>(ovalues->front())) {
      for (const auto& v : *ovalues) {
        if (!is<std::vector<Parameter>>(v)) {
          return ctx.registerErrorMessage(
              std::string{"invalid option 'CoarseComponents' of "} + n);
        }
        auto ol = read(get<std::vector<Parameter>>(throwing, v));
        if (isInvalid(ol)) {
          return {};
        }
        levels.push_back(std::move(*ol));
      }
    } else {
      auto ol = read(*ovalues);
      if (isInvalid(ol)) {
        return {};
      }
      levels.push_back(std::move(*ol));
    }
    // each level keeps a non empty subset of the components of the previous
    // one, and does not keep all of them
    for (std::size_t i = 0; i != levels.size(); ++i) {
      const auto& l = levels[i];
      const auto previous =
          i == 0 ? static_cast<std::size_t>(vdim) : levels[i - 1].size();
      const auto subset =
          (i == 0) || std::ranges::all_of(l, [&levels, i](const auto c) {
            return std::ranges::find(levels[i - 1], c) != levels[i - 1].end();
          });
      if (l.empty() || (l.size() >= previous) || (!subset)) {
        return ctx.registerErrorMessage(
            std::string{"invalid option 'CoarseComponents' of "} + n +
            ": each level must keep a non empty strict subset of the "
            "components of the previous one");
      }
    }
    return levels;
  }  // end of getMGRCoarseComponents

  /*!
   * \return if the values of the integer options of the multigrid reduction
   * of hypre are among the ones documented by hypre (hypre 3.2.0) and, when
   * hypre computes on a GPU, among the ones it computes there: the other
   * ones apply the relaxations, interpolations or restrictions of the host
   * to the matrices on the GPU, and stop the program (Gauss-Seidel
   * relaxations, classical and approximate inverse interpolations and
   * restrictions), or the V-cycle relaxation (1) frees memory which is not
   * its own
   * \param[in, out] ctx: execution context
   * \param[in] opts: options
   */
  [[nodiscard]] static bool checkHypreMGRValues(Context& ctx,
                                                const Parameters& opts) {
    struct Values {
      //! \brief name of the option
      const char* name;
      //! \brief values documented by hypre
      std::vector<int> documented;
      //! \brief values computed by hypre on a GPU
      std::vector<int> gpu;
    };
    const auto values = std::array<Values, 3>{
        {{"RelaxationType",
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 18, 99, 199},
          {0, 2, 7, 9, 18, 99, 199}},
         {"InterpolationType", {0, 1, 2, 3, 4, 12, 13, 14}, {0, 1, 2, 12}},
         {"RestrictionType",
          {0, 1, 2, 3, 4, 5, 12, 13, 14, 15},
          {0, 1, 2, 12}}}};
    const auto gpu = isHypreUsingGPU();
    for (const auto& v : values) {
      if (!contains(opts, v.name)) {
        continue;
      }
      const auto ov = get<int>(ctx, opts, v.name);
      if (isInvalid(ov)) {
        return false;
      }
      const auto& accepted = gpu ? v.gpu : v.documented;
      if (std::ranges::find(accepted, *ov) == accepted.end()) {
        auto msg = "invalid value " + std::to_string(*ov) + " for option '" +
                   std::string{v.name} +
                   "' of the multigrid reduction of hypre" +
                   (gpu ? " computing on the GPU" : "") + " (expected ";
        for (std::size_t i = 0; i != accepted.size(); ++i) {
          msg += (i == 0 ? "" : (i + 1 == accepted.size() ? " or " : ", ")) +
                 std::to_string(accepted[i]);
        }
        return ctx.registerErrorMessage(msg + ")");
      }
    }
    return true;
  }  // end of checkHypreMGRValues

  /*!
   * \return if the block interpolations and restrictions of the multigrid
   * reduction of hypre can be used with the given levels: they group
   * consecutive fine unknowns, which are the components eliminated at a
   * node only with the ordering by vector dimension, and, with hypre 3.2.0
   * (tested on a beam in linear elasticity), the interpolation 13 gives
   * wrong results, the restriction 12 stops the program when a level
   * eliminates several components, the restrictions 14 and 15 give invalid
   * values when a level eliminates one component, and, when hypre computes
   * on a GPU, the interpolation 12 accesses invalid memory when a level
   * eliminates several components
   * \param[in, out] ctx: execution context
   * \param[in] eliminated: number of components eliminated at each level
   * \param[in] interpolation: interpolation
   * \param[in] restriction: restriction
   * \param[in] by_nodes: if true, the unknowns are ordered by nodes
   */
  [[nodiscard]] static bool checkHypreMGRBlockOptions(
      Context& ctx,
      const std::vector<std::size_t>& eliminated,
      const int interpolation,
      const int restriction,
      const bool by_nodes) {
    const auto n = std::string{"the multigrid reduction of hypre"};
    const auto several =
        std::ranges::any_of(eliminated, [](const auto e) { return e > 1; });
    const auto one =
        std::ranges::any_of(eliminated, [](const auto e) { return e == 1; });
    const auto blocks = (interpolation >= 12) || (restriction >= 12);
    if (interpolation == 13) {
      return ctx.registerErrorMessage(
          "the interpolation 13 of " + n +
          " is refused: it gives wrong results (hypre 3.2.0)");
    }
    if (by_nodes && blocks && several) {
      return ctx.registerErrorMessage(
          "the block interpolations and restrictions of " + n +
          " (12 to 15) group consecutive fine unknowns: they require "
          "unknowns ordered by vector dimension when a level eliminates "
          "several components");
    }
    if ((restriction == 12) && several) {
      return ctx.registerErrorMessage(
          "the restriction 12 of " + n +
          " is refused when a level eliminates several components: it stops "
          "the program (hypre 3.2.0)");
    }
    if (((restriction == 14) || (restriction == 15)) && one) {
      return ctx.registerErrorMessage(
          "the restrictions 14 and 15 of " + n +
          " are refused when a level eliminates one component: they give "
          "invalid values (hypre 3.2.0)");
    }
    if (isHypreUsingGPU() && (interpolation == 12) && several) {
      return ctx.registerErrorMessage(
          "the interpolation 12 of " + n +
          " computing on the GPU is refused when a level eliminates several "
          "components: it accesses invalid memory (hypre 3.2.0)");
    }
    return true;
  }  // end of checkHypreMGRBlockOptions

  /*!
   * \return the multigrid reduction of hypre (MGR)
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options: `CoarseComponents` (see
   * `getMGRCoarseComponents`), and the options of hypre (see its
   * documentation): `RelaxationType`, `NumberOfRelaxationSweeps`,
   * `InterpolationType`, `RestrictionType` and `VerbosityLevel` (0 to 7,
   * information on the setup, the resolution and the parameters: the other
   * bits of the print level of MGR write files)
   *
   * The block interpolations and restrictions are checked by
   * `checkHypreMGRBlockOptions`. When hypre
   * computes on a GPU, with the Jacobi relaxation (0, by default) and the
   * interpolation 12, the levels eliminating one component use the Jacobi
   * relaxation of the GPU (7): hypre 3.2.0 replaces their interpolation 12
   * by 2 after having replaced the relaxation 0 by 7 only for the
   * interpolations other than 12, and applies the Jacobi relaxation of the
   * host to the matrices on the GPU.
   */
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  setHypreMGRPreconditioner(Context& ctx,
                            FiniteElementSpace<true>& fespace,
                            const Parameters& opts) {
    using Problem = AbstractNonLinearEvolutionProblem;
    constexpr auto m = hypre_option_maximum;
    const auto integers = std::array<HypreOption<HYPRE_Int>, 5>{
        {{"RelaxationType", 0, 199, HYPRE_MGRSetRelaxType},
         {"NumberOfRelaxationSweeps", 1, m, HYPRE_MGRSetNumRelaxSweeps},
         {"InterpolationType", 0, 14, HYPRE_MGRSetInterpType},
         {"RestrictionType", 0, 15, HYPRE_MGRSetRestrictType},
         {Problem::SolverVerbosityLevel, 0, 7, HYPRE_MGRSetPrintLevel}}};
    auto allowed = std::vector<std::string>{"CoarseComponents"};
    for (const auto& o : integers) {
      allowed.push_back(o.name);
    }
    if (!checkParameters(ctx, opts, allowed)) {
      return {};
    }
    if (!checkHypreMGRValues(ctx, opts)) {
      return {};
    }
    const auto vdim = fespace.GetVDim();
    auto olevels = getMGRCoarseComponents(ctx, opts, vdim);
    if (isInvalid(olevels)) {
      return {};
    }
    auto settings = std::vector<internals::HypreSetting>{};
    if (!addHypreSettings(ctx, settings, opts, integers,
                          "the multigrid reduction of hypre")) {
      return {};
    }
    // number of components eliminated at each level
    auto eliminated = std::vector<std::size_t>{};
    for (std::size_t l = 0; l != olevels->size(); ++l) {
      const auto previous =
          (l == 0) ? static_cast<std::size_t>(vdim) : (*olevels)[l - 1].size();
      eliminated.push_back(previous - (*olevels)[l].size());
    }
    const auto ointerpolation = get_if<int>(ctx, opts, "InterpolationType", 2);
    const auto orestriction = get_if<int>(ctx, opts, "RestrictionType", 0);
    const auto orelaxation = get_if<int>(ctx, opts, "RelaxationType", 0);
    if (isInvalid(ointerpolation) || isInvalid(orestriction) ||
        isInvalid(orelaxation)) {
      return {};
    }
    const auto by_nodes = fespace.GetOrdering() == mfem::Ordering::byNODES;
    if (!checkHypreMGRBlockOptions(ctx, eliminated, *ointerpolation,
                                   *orestriction, by_nodes)) {
      return {};
    }
    auto level_relaxation_types = std::vector<HYPRE_Int>{};
    if ((isHypreUsingGPU()) && (*orelaxation == 0) && (*ointerpolation == 12)) {
      for (const auto e : eliminated) {
        level_relaxation_types.push_back(e == 1 ? 7 : 0);
      }
    }
    return std::make_unique<internals::HypreMGR>(
        std::move(settings), std::move(*olevels), vdim, by_nodes,
        std::move(level_relaxation_types));
  }  // end of setHypreMGRPreconditioner

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
            "' only accepts the preconditioners applied by hypre (the "
            "'HypreSmoother' preconditioner is applied by MFEM)");
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

  /*!
   * \brief set the absolute tolerance of a Krylov solver of hypre
   * \param[in, out] s: solver
   * \param[in] tol: tolerance
   */
  template <typename SolverType>
  static void setHypreAbsoluteTolerance(SolverType& s, const real tol) {
    if constexpr (std::is_same_v<SolverType, mfem::HypreFGMRES>) {
      // mfem::HypreFGMRES has no `SetAbsTol` method
      HYPRE_ParCSRFlexGMRESSetAbsoluteTol(s, tol);
    } else {
      s.SetAbsTol(tol);
    }
  }  // end of setHypreAbsoluteTolerance

  /*!
   * \brief set the dimension of the Krylov space of a GMRES solver of hypre
   * \param[in, out] s: solver
   * \param[in] n: dimension
   */
  template <typename SolverType>
  static void setHypreKrylovDimension(SolverType& s, const int n) {
    if constexpr (std::is_same_v<SolverType, internals::HypreKrylovSolver>) {
      std::ignore = s.SetKDim(n);
    } else if constexpr (requires { s.SetKDim(n); }) {
      s.SetKDim(n);
    }
  }  // end of setHypreKrylovDimension

  /*!
   * \return the generator of a Krylov solver of hypre, whose options are its
   * verbosity level, its relative (`Tolerance`) and absolute tolerances, its
   * maximum number of iterations, its preconditioner (a preconditioner of
   * hypre) and, for the GMRES solvers, the dimension of the Krylov space
   * (`KDim`)
   * \param[in] n: name of the solver
   * \param[in] d: description of the solver
   * \param[in] make: function building the solver from the communicator
   */
  template <typename SolverType>
  [[nodiscard]] static LinearSolverFactory<true>::Generator
  buildHypreKrylovSolverGenerator(
      const std::string n,
      const LinearSolverDescription d,
      std::function<std::unique_ptr<SolverType>(MPI_Comm)> make) {
    return [n, d, make](Context& ctx, FiniteElementSpace<true>& fespace,
                        const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      if (!checkParameters(ctx, params, d.options)) {
        return {};
      }
      // the solvers which only compute on the host can't use the matrices
      // and the vectors of hypre stored on a GPU
      if ((d.placement == LinearSolverPlacement::HOST) && (isHypreUsingGPU())) {
        std::ignore = ctx.registerErrorMessage(
            "the linear solver '" + n +
            "' only computes on the host, which hypre does not");
        return {};
      }
      auto s = make(fespace.GetComm());
      s->iterative_mode = false;
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
        setHypreAbsoluteTolerance(*s, *otol);
      }
      if (contains(params, "Tolerance")) {
        const auto otol = get<double>(ctx, params, "Tolerance");
        if (isInvalid(otol)) {
          return {};
        }
        s->SetTol(boundRelativeTolerance(*otol));
      }
      if (contains(params, "KDim")) {
        const auto odim = get<int>(ctx, params, "KDim");
        if (isInvalid(odim)) {
          return {};
        }
        if (*odim < 1) {
          std::ignore = ctx.registerErrorMessage(
              "invalid Krylov dimension (must be strictly positive)");
          return {};
        }
        setHypreKrylovDimension(*s, *odim);
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        const auto oiterMax =
            get<int>(ctx, params, Problem::SolverMaximumNumberOfIterations);
        if (isInvalid(oiterMax)) {
          return {};
        }
        s->SetMaxIter(*oiterMax);
      }
      auto oprec =
          setHypreKrylovSolverPreconditioner(ctx, fespace, n, d, *s, params);
      if (isInvalid(oprec)) {
        return {};
      }
      return LinearSolverHandler{std::move(s), std::move(*oprec)};
    };
  }  // end of buildHypreKrylovSolverGenerator

  /*!
   * \return the generator of the hybrid solver of hypre, whose options are
   * its verbosity level, its relative tolerance (`Tolerance`), its maximum
   * number of iterations (of each phase), its Krylov solver (`SolverType`:
   * `PCG`, by default, `GMRES` or `BiCGSTAB`), the dimension of the Krylov
   * space of the GMRES solver (`KDim`), the convergence rate below which
   * the diagonal scaling is replaced by BoomerAMG (`ConvergenceRate`, 0.9
   * by default) and options of BoomerAMG (`StrongThreshold`,
   * `CoarseningType`, `InterpolationType`, `RelaxationType`,
   * `AggressiveCoarseningLevels`, `MaximumNumberOfLevels`)
   * \param[in] d: description of the solver
   */
  [[nodiscard]] static LinearSolverFactory<true>::Generator
  buildHypreHybridSolverGenerator(const LinearSolverDescription d) {
    return [d](Context& ctx, FiniteElementSpace<true>&,
               const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      using Setting = internals::HypreSetting;
      if (!checkParameters(ctx, params, d.options)) {
        return {};
      }
      auto settings = std::vector<Setting>{};
      const auto otol = get_if<double>(ctx, params, "Tolerance", 1e-6);
      const auto omaxit = get_if<int>(
          ctx, params, Problem::SolverMaximumNumberOfIterations, 1000);
      if (isInvalid(otol) || isInvalid(omaxit)) {
        return {};
      }
      // the ranges accepted by hypre, which ignores the other values
      if ((*otol < 0) || (*otol > 1) || (*omaxit < 0)) {
        std::ignore = ctx.registerErrorMessage(
            "invalid tolerance or maximum number of iterations of the "
            "hybrid solver of hypre");
        return {};
      }
      if (contains(params, "SolverType")) {
        const auto otype = get<std::string>(ctx, params, "SolverType");
        if (isInvalid(otype)) {
          return {};
        }
        const auto types = std::array<std::pair<const char*, int>, 3>{
            {{"PCG", 1}, {"GMRES", 2}, {"BiCGSTAB", 3}}};
        const auto p =
            std::find_if(types.begin(), types.end(),
                         [&otype](const auto& t) { return *otype == t.first; });
        if (p == types.end()) {
          std::ignore = ctx.registerErrorMessage(
              "invalid type '" + *otype +
              "' of the Krylov solver of the hybrid solver of hypre");
          return {};
        }
        settings.push_back([v = p->second](HYPRE_Solver s) {
          HYPRE_ParCSRHybridSetSolverType(s, v);
        });
      }
      constexpr auto m = hypre_option_maximum;
      const auto integers = std::array<HypreOption<HYPRE_Int>, 7>{
          {{Problem::SolverVerbosityLevel, 0, m,
            HYPRE_ParCSRHybridSetPrintLevel},
           {"KDim", 1, m, HYPRE_ParCSRHybridSetKDim},
           {"CoarseningType", 0, m, HYPRE_ParCSRHybridSetCoarsenType},
           {"InterpolationType", 0, m, HYPRE_ParCSRHybridSetInterpType},
           {"RelaxationType", 0, m, HYPRE_ParCSRHybridSetRelaxType},
           {"AggressiveCoarseningLevels", 0, m,
            HYPRE_ParCSRHybridSetAggNumLevels},
           {"MaximumNumberOfLevels", 1, m, HYPRE_ParCSRHybridSetMaxLevels}}};
      const auto reals = std::array<HypreOption<HYPRE_Real>, 2>{
          {{"ConvergenceRate", 0, 1, HYPRE_ParCSRHybridSetConvergenceTol},
           {"StrongThreshold", 0, 1, HYPRE_ParCSRHybridSetStrongThreshold}}};

      const auto n = "the hybrid solver of hypre";
      if ((!addHypreSettings(ctx, settings, params, integers, n)) ||
          (!addHypreSettings(ctx, settings, params, reals, n))) {
        return {};
      }
      return LinearSolverHandler{
          std::make_unique<internals::HypreHybridSolver>(
              std::move(settings), boundRelativeTolerance(*otol), *omaxit),
          std::unique_ptr<LinearSolverPreconditioner>{}};
    };
  }  // end of buildHypreHybridSolverGenerator
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
    using Problem = AbstractNonLinearEvolutionProblem;
    const auto krylov = [](const bool gmres, std::string description) {
      return LinearSolverDescription{
          .library = "hypre",
          .kind = LinearSolverKind::ITERATIVE,
          .placement = LinearSolverPlacement::HYPRE,
          .preconditioners = AcceptedPreconditioners::SAME_LIBRARY,
          .options = getHypreKrylovSolverOptions(gmres),
          .description = std::move(description)};
    };
    // Krylov solvers wrapped by MFEM
    const auto add = [&f]<typename SolverType>(
                         const std::string& n,
                         const LinearSolverDescription& d) {
      declare(
          f, n, d,
          buildHypreKrylovSolverGenerator<SolverType>(n, d, [](MPI_Comm comm) {
            return std::make_unique<SolverType>(comm);
          }));
    };
    add.template operator()<mfem::HyprePCG>(
        "HyprePCG", krylov(false, "conjugate gradient of hypre"));
    add.template operator()<mfem::HypreGMRES>("HypreGMRES",
                                              krylov(true, "GMRES of hypre"));
    add.template operator()<mfem::HypreFGMRES>(
        "HypreFGMRES", krylov(true, "flexible GMRES of hypre"));
    // Krylov solvers of hypre which are not wrapped by MFEM
    using internals::HypreKrylovSolver;
    using internals::HypreKrylovSolverFunctions;
    const auto add_native = [&f](const std::string& n,
                                 const LinearSolverDescription& d,
                                 const HypreKrylovSolverFunctions& functions) {
      declare(f, n, d,
              buildHypreKrylovSolverGenerator<HypreKrylovSolver>(
                  n, d, [functions](MPI_Comm comm) {
                    return std::make_unique<HypreKrylovSolver>(comm, functions);
                  }));
    };
    add_native("HypreBiCGSTAB", krylov(false, "BiCGSTAB of hypre"),
               internals::getHypreBiCGSTABFunctions());
    // the products of COGMRES only compute on the host (segmentation fault
    // in hypre_SeqVectorMassInnerProd with hypre 3.2.0 on a GPU)
    auto cogmres = krylov(true,
                          "GMRES of hypre with classical Gram-Schmidt "
                          "orthogonalization (communication optimal)");
    cogmres.placement = LinearSolverPlacement::HOST;
    add_native("HypreCOGMRES", cogmres, internals::getHypreCOGMRESFunctions());
    add_native("HypreLGMRES",
               krylov(true,
                      "GMRES of hypre augmented with approximations of "
                      "the error (LGMRES)"),
               internals::getHypreLGMRESFunctions());
    const auto hybrid = LinearSolverDescription{
        .library = "hypre",
        .kind = LinearSolverKind::ITERATIVE,
        .placement = LinearSolverPlacement::HYPRE,
        .preconditioners = AcceptedPreconditioners::NONE,
        // no absolute tolerance: hypre bounds it by one, and the convergence
        // is checked on the relative residual (see `HypreHybridSolver`)
        .options = {Problem::SolverVerbosityLevel, "Tolerance",
                    Problem::SolverMaximumNumberOfIterations, "SolverType",
                    "KDim", "ConvergenceRate", "StrongThreshold",
                    "CoarseningType", "InterpolationType", "RelaxationType",
                    "AggressiveCoarseningLevels", "MaximumNumberOfLevels"},
        .description =
            "hybrid solver of hypre: Krylov solver preconditioned by the "
            "diagonal scaling, then by BoomerAMG if the convergence is slow"};
    declare(f, "HypreHybrid", hybrid, buildHypreHybridSolverGenerator(hybrid));
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
    declare(
        f, "HypreSmoother",
        preconditioner(LinearSolverPlacement::HYPRE,
                       {"Type", "RelaxationSweeps", "RelaxationWeight", "Omega",
                        "PolynomialOrder", "PolynomialFraction",
                        "EigenvalueEstimationIterations", "TaubinLambda",
                        "TaubinMu", "TaubinIterations", "Window"},
                       "smoothers of hypre (Jacobi, l1-Jacobi, Gauss-Seidel, "
                       "polynomial smoothers)"),
        setHypreSmootherPreconditioner);
    declare(f, "HypreFSAI",
            preconditioner(
                LinearSolverPlacement::HYPRE,
                {"AlgorithmType", "LocalSolveType", "MaximumSteps",
                 "MaximumStepSize", "MaximumNonZerosPerRow", "NumberOfLevels",
                 "Threshold", "KapTolerance", "EigenvalueEstimationIterations",
                 Problem::SolverVerbosityLevel},
                "factorized sparse approximate inverse of hypre"),
            setHypreFSAIPreconditioner);
    declare(f, "HypreMGR",
            preconditioner(LinearSolverPlacement::HYPRE,
                           {"CoarseComponents", "RelaxationType",
                            "NumberOfRelaxationSweeps", "InterpolationType",
                            "RestrictionType", Problem::SolverVerbosityLevel},
                           "multigrid reduction of hypre (MGR), the unknowns "
                           "being split by their components"),
            setHypreMGRPreconditioner);
  }  // end of declareHyprePreconditioners

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_USE_MPI */
