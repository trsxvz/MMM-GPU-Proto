/*!
 * \file   src/SuperLUSolvers.cxx
 * \brief  This file implements the module of the direct solver of
 * SuperLU_DIST (see `LinearSolverModules.hxx`)
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_SUPERLU)

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <optional>
#include <algorithm>
#include <functional>
#include "mfem/linalg/superlu.hpp"
#include "mfem/fem/pfespace.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief direct solver of SuperLU_DIST. The matrices of hypre are copied
   * to the format of SuperLU_DIST when they are set, and given to the solver
   * of MFEM just before the resolution, so that the solver of MFEM never
   * holds an operator which is not factorized (it then refuses a new one).
   * The solver of MFEM is built, with the settings, at the first resolution
   * and when the global size of the operator changes; otherwise, the next
   * operators are factorized from scratch, or with the column permutation
   * of the previous factorization if the analysis is reused
   * (`SamePattern`, which requires the same sparsity pattern).
   */
  struct SuperLUDirectSolver final : public mfem::Solver {
    //! \brief a setting of the solver
    using Setting = std::function<void(mfem::SuperLUSolver&)>;
    /*!
     * \param[in] c: communicator
     * \param[in] s: settings
     * \param[in] r: if true, the analysis is reused
     */
    SuperLUDirectSolver(MPI_Comm c, std::vector<Setting> s, const bool r)
        : communicator(c), settings(std::move(s)), reuse_analysis(r) {}
    void SetOperator(const mfem::Operator& op) override {
      if (dynamic_cast<const mfem::HypreParMatrix*>(&op) == nullptr) {
        mgis::raise(
            "SuperLUDirectSolver::SetOperator: the operator is not a matrix "
            "of hypre");
      }
      this->pending = std::make_unique<mfem::SuperLURowLocMatrix>(op);
      this->height = op.Height();
      this->width = op.Width();
    }
    void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
      if (this->pending != nullptr) {
        // the global size is the same on all processes, which take the
        // same decision (the creation of the solver of MFEM is collective);
        // an operator of the same global size whose rows are distributed
        // differently is not supported (mfem-mgis never repartitions)
        const auto n = this->pending->GetGlobalNumRows();
        if ((this->solver == nullptr) || (n != this->global_size)) {
          this->solver.reset();
          this->solver =
              std::make_unique<mfem::SuperLUSolver>(this->communicator);
          for (const auto& setting : this->settings) {
            setting(*(this->solver));
          }
        } else if (this->reuse_analysis) {
          this->solver->SetFact(mfem::superlu::SamePattern);
        }
        this->solver->SetOperator(*(this->pending));
        // the previous matrix is referenced by the solver until here
        this->matrix = std::move(this->pending);
        this->global_size = n;
      }
      if (this->solver == nullptr) {
        mgis::raise("SuperLUDirectSolver::Mult: the operator is not set");
      }
      this->solver->Mult(x, y);
    }

   private:
    //! \brief communicator
    const MPI_Comm communicator;
    //! \brief settings
    const std::vector<Setting> settings;
    //! \brief if true, the analysis is reused
    const bool reuse_analysis;
    //! \brief matrix set and not yet given to the solver of MFEM
    mutable std::unique_ptr<mfem::SuperLURowLocMatrix> pending;
    /*!
     * \brief matrix of the solver of MFEM, declared before the solver: the
     * destructor of the solver of MFEM reads it
     */
    mutable std::unique_ptr<mfem::SuperLURowLocMatrix> matrix;
    //! \brief solver of MFEM
    mutable std::unique_ptr<mfem::SuperLUSolver> solver;
    //! \brief global size of the matrix of the solver of MFEM
    mutable HYPRE_BigInt global_size = 0;
  };  // end of struct SuperLUDirectSolver

  /*!
   * \brief look for the value of an option given by a name
   * \return the value, or an empty value on error
   * \param[in, out] ctx: execution context
   * \param[in] params: parameters
   * \param[in] n: name of the option
   * \param[in] values: names and values
   */
  template <typename ValueType, std::size_t N>
  [[nodiscard]] static std::optional<std::optional<ValueType>>
  getEnumerationOption(
      Context& ctx,
      const Parameters& params,
      const char* const n,
      const std::array<std::pair<const char*, ValueType>, N>& values) {
    if (!contains(params, n)) {
      return std::optional<ValueType>{};
    }
    const auto ov = get<std::string>(ctx, params, n);
    if (isInvalid(ov)) {
      return {};
    }
    const auto p =
        std::find_if(values.begin(), values.end(),
                     [&ov](const auto& v) { return *ov == v.first; });
    if (p == values.end()) {
      std::ignore = ctx.registerErrorMessage(
          "invalid value '" + *ov + "' of the option '" + std::string{n} +
          "' of SuperLU_DIST");
      return {};
    }
    return std::optional<ValueType>{p->second};
  }  // end of getEnumerationOption

  //! \return the options of the direct solver of SuperLU_DIST
  [[nodiscard]] static std::vector<std::string> getSuperLUSolverOptions() {
    using Problem = AbstractNonLinearEvolutionProblem;
    return {"ReuseAnalysis",
            "ColumnPermutation",
            "RowPermutation",
            "IterativeRefinement",
            "Equilibrate",
            "ReplaceTinyPivot",
            "SymmetricPattern",
            "ParallelSymbolicFactorization",
            "NumberOfLookAheads",
            "DeviceOffload",
            Problem::SolverVerbosityLevel};
  }  // end of getSuperLUSolverOptions

  /*!
   * \return the direct solver of SuperLU_DIST
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] params: options (see `declareSuperLULinearSolvers`)
   */
  [[nodiscard]] static LinearSolverHandler makeSuperLUDirectSolver(
      Context& ctx,
      FiniteElementSpace<true>& fespace,
      const Parameters& params) {
    using Problem = AbstractNonLinearEvolutionProblem;
    using Setting = SuperLUDirectSolver::Setting;
    namespace superlu = mfem::superlu;
    if (!checkParameters(ctx, params, getSuperLUSolverOptions())) {
      return {};
    }
    auto settings = std::vector<Setting>{};
    // the statistics are printed by default by MFEM
    const auto olvl =
        get_if<int>(ctx, params, Problem::SolverVerbosityLevel, 0);
    if (isInvalid(olvl)) {
      return {};
    }
    settings.push_back(
        [v = *olvl > 0](mfem::SuperLUSolver& s) { s.SetPrintStatistics(v); });
    const auto ocolumns = getEnumerationOption(
        ctx, params, "ColumnPermutation",
        std::array<std::pair<const char*, superlu::ColPerm>, 6>{
            {{"Natural", superlu::NATURAL},
             {"MMD_ATA", superlu::MMD_ATA},
             {"MMD_AT_PLUS_A", superlu::MMD_AT_PLUS_A},
             {"COLAMD", superlu::COLAMD},
             {"METIS_AT_PLUS_A", superlu::METIS_AT_PLUS_A},
             {"ParMETIS", superlu::PARMETIS}}});
    const auto orows = getEnumerationOption(
        ctx, params, "RowPermutation",
        std::array<std::pair<const char*, superlu::RowPerm>, 3>{
            {{"None", superlu::NOROWPERM},
             {"LargeDiag_MC64", superlu::LargeDiag_MC64},
             {"LargeDiag_HWPM", superlu::LargeDiag_HWPM}}});
    const auto orefinement = getEnumerationOption(
        ctx, params, "IterativeRefinement",
        std::array<std::pair<const char*, superlu::IterRefine>, 4>{
            {{"None", superlu::NOREFINE},
             {"Single", superlu::SLU_SINGLE},
             {"Double", superlu::SLU_DOUBLE},
             {"Extra", superlu::SLU_EXTRA}}});
    if (isInvalid(ocolumns) || isInvalid(orows) || isInvalid(orefinement)) {
      return {};
    }
    if (ocolumns->has_value()) {
      settings.push_back([v = **ocolumns](mfem::SuperLUSolver& s) {
        s.SetColumnPermutation(v);
      });
    }
    if (orows->has_value()) {
      settings.push_back(
          [v = **orows](mfem::SuperLUSolver& s) { s.SetRowPermutation(v); });
    }
    if (orefinement->has_value()) {
      settings.push_back([v = **orefinement](mfem::SuperLUSolver& s) {
        s.SetIterativeRefine(v);
      });
    }
    const auto booleans =
        std::array<std::pair<const char*, void (mfem::SuperLUSolver::*)(bool)>,
                   5>{
            {{"Equilibrate", &mfem::SuperLUSolver::SetEquilibriate},
             {"ReplaceTinyPivot", &mfem::SuperLUSolver::SetReplaceTinyPivot},
             {"SymmetricPattern", &mfem::SuperLUSolver::SetSymmetricPattern},
             {"ParallelSymbolicFactorization",
              &mfem::SuperLUSolver::SetParSymbFact},
             {"DeviceOffload", &mfem::SuperLUSolver::SetDeviceOffload}}};
    for (const auto& [n, f] : booleans) {
      if (!contains(params, n)) {
        continue;
      }
      const auto ov = get<bool>(ctx, params, n);
      if (isInvalid(ov)) {
        return {};
      }
      settings.push_back([f, v = *ov](mfem::SuperLUSolver& s) { (s.*f)(v); });
    }
    if (contains(params, "NumberOfLookAheads")) {
      const auto ov = get<int>(ctx, params, "NumberOfLookAheads");
      if (isInvalid(ov)) {
        return {};
      }
      if (*ov < 0) {
        std::ignore = ctx.registerErrorMessage(
            "invalid number of look-aheads of SuperLU_DIST");
        return {};
      }
      settings.push_back(
          [v = *ov](mfem::SuperLUSolver& s) { s.SetNumLookAheads(v); });
    }
    const auto oreuse = get_if<bool>(ctx, params, "ReuseAnalysis", false);
    const auto oparallel_symbolic =
        get_if<bool>(ctx, params, "ParallelSymbolicFactorization", false);
    if (isInvalid(oreuse) || isInvalid(oparallel_symbolic)) {
      return {};
    }
    // the parallel symbolic factorization of SuperLU_DIST, selected by the
    // ParMETIS ordering (unless explicitly disabled) and selecting it, fails
    // (segmentation fault in symbfact_dist with SuperLU_DIST 9.2.1) on one
    // process and when the analysis is reused
    const auto parmetis =
        (ocolumns->has_value()) && (**ocolumns == superlu::PARMETIS);
    const auto parallel_symbolic =
        *oparallel_symbolic ||
        (parmetis && (!contains(params, "ParallelSymbolicFactorization")));
    if (parallel_symbolic) {
      if (*oreuse) {
        std::ignore = ctx.registerErrorMessage(
            "the analysis of SuperLU_DIST can't be reused with its parallel "
            "symbolic factorization (column permutation 'ParMETIS' or option "
            "'ParallelSymbolicFactorization')");
        return {};
      }
      if (fespace.GetNRanks() < 2) {
        std::ignore = ctx.registerErrorMessage(
            "the parallel symbolic factorization of SuperLU_DIST (column "
            "permutation 'ParMETIS' or option "
            "'ParallelSymbolicFactorization') requires at least two "
            "processes");
        return {};
      }
    }
    return LinearSolverHandler{
        std::make_unique<SuperLUDirectSolver>(fespace.GetComm(),
                                              std::move(settings), *oreuse),
        std::unique_ptr<LinearSolverPreconditioner>{}};
  }  // end of makeSuperLUDirectSolver

  void declareSuperLULinearSolvers(LinearSolverFactory<true>& f) {
    declare(f, "SuperLUSolver",
            LinearSolverDescription{
                .library = "SuperLU_DIST",
                .kind = LinearSolverKind::DIRECT,
                .placement = LinearSolverPlacement::HOST,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = getSuperLUSolverOptions(),
                .description = "direct solver of SuperLU_DIST"},
            makeSuperLUDirectSolver);
  }  // end of declareSuperLULinearSolvers

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_SUPERLU) */
