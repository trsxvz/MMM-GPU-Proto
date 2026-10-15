/*!
 * \file   src/MUMPSSolvers.cxx
 * \brief  This file implements the module of the direct solver of MUMPS
 * (see `LinearSolverModules.hxx`), and the MUMPS solver which reuses its
 * analysis for the matrices of the same pattern
 * \date   12/10/2026
 */

#include "MUMPSSolvers.hxx"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS)

#include <array>
#include <utility>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include "mfem/fem/pfespace.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis::internals {

  MUMPSSolverReusingAnalysis::MUMPSSolverReusingAnalysis(const MPI_Comm c)
      : mfem::MUMPSSolver(c), communicator(c) {}

  void MUMPSSolverReusingAnalysis::SetOperator(const mfem::Operator& op) {
    const auto* const A = dynamic_cast<const mfem::HypreParMatrix*>(&op);
    if (A == nullptr) {
      // error reported by MFEM
      mfem::MUMPSSolver::SetOperator(op);
      return;
    }
    const auto same =
        this->pattern.isSameOnAllProcesses(*A, this->communicator);
    this->SetReorderingReuse(same);
    // the pattern is forgotten until the factorization succeeds: if it
    // fails by an exception, the next matrix is analysed again
    auto previous = std::exchange(this->pattern, HypreMatrixPatternRecord{});
    mfem::MUMPSSolver::SetOperator(op);
    if (same) {
      this->pattern = std::move(previous);
      return;
    }
    ++(this->number_of_analyses);
    this->pattern.record(*A);
  }  // end of SetOperator

  size_type MUMPSSolverReusingAnalysis::getNumberOfAnalyses() const noexcept {
    return this->number_of_analyses;
  }  // end of getNumberOfAnalyses

  //! \return the options of the MUMPS solver
  [[nodiscard]] static std::vector<std::string> getMUMPSSolverOptions() {
    using Problem = AbstractNonLinearEvolutionProblem;
    return {"Symmetric",     "PositiveDefinite",
            "ReuseAnalysis", "ReorderingStrategy",
            "BLRTolerance",  Problem::SolverVerbosityLevel};
  }  // end of getMUMPSSolverOptions

  /*!
   * \return the generator of the MUMPS solver, whose options are:
   * - `Symmetric` and `PositiveDefinite`: type of the matrix;
   * - `ReuseAnalysis`: if true, the analysis is reused for the matrices of
   *   the same pattern;
   * - `ReorderingStrategy`: `Automatic` (by default), `AMD`, `AMF`, `PORD`,
   *   `METIS`, `ParMETIS`, `Scotch` or `PTScotch`;
   * - `BLRTolerance`: tolerance of the block low-rank factorization (zero,
   *   an exact factorization, by default; an approximate factorization when
   *   positive, whose result is not checked);
   * - `VerbosityLevel`: print level of MUMPS (`ICNTL(4)`, one, errors only,
   *   by default).
   */
  [[nodiscard]] static std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildMUMPSSolverGenerator() {
    return [](Context& ctx, FiniteElementSpace<true>& fespace,
              const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      if (!checkParameters(ctx, params, getMUMPSSolverOptions())) {
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
      // errors only, by default
      const auto olvl =
          get_if<int>(ctx, params, Problem::SolverVerbosityLevel, 1);
      if (isInvalid(olvl)) {
        return {};
      }
      s->SetPrintLevel(*olvl);
      if (contains(params, "ReorderingStrategy")) {
        using MUMPSSolver = mfem::MUMPSSolver;
        const auto strategies =
            std::array<std::pair<const char*, MUMPSSolver::ReorderingStrategy>,
                       8>{{{"Automatic", MUMPSSolver::AUTOMATIC},
                           {"AMD", MUMPSSolver::AMD},
                           {"AMF", MUMPSSolver::AMF},
                           {"PORD", MUMPSSolver::PORD},
                           {"METIS", MUMPSSolver::METIS},
                           {"ParMETIS", MUMPSSolver::PARMETIS},
                           {"Scotch", MUMPSSolver::SCOTCH},
                           {"PTScotch", MUMPSSolver::PTSCOTCH}}};
        const auto ostrategy =
            get<std::string>(ctx, params, "ReorderingStrategy");
        if (isInvalid(ostrategy)) {
          return {};
        }
        const auto p = std::find_if(
            strategies.begin(), strategies.end(),
            [&ostrategy](const auto& v) { return *ostrategy == v.first; });
        if (p == strategies.end()) {
          return ctx.registerErrorMessage("invalid reordering strategy '" +
                                          *ostrategy + "' of MUMPS");
        }
        s->SetReorderingStrategy(p->second);
      }
      if (contains(params, "BLRTolerance")) {
#if MFEM_MUMPS_VERSION >= 510
        const auto otol = get<double>(ctx, params, "BLRTolerance");
        if (isInvalid(otol)) {
          return {};
        }
        if (*otol < 0) {
          return ctx.registerErrorMessage(
              "invalid tolerance of the block low-rank factorization of "
              "MUMPS");
        }
        s->SetBLRTol(*otol);
#else  /* MFEM_MUMPS_VERSION >= 510 */
        return ctx.registerErrorMessage(
            "the block low-rank factorization of MUMPS requires MUMPS 5.1");
#endif /* MFEM_MUMPS_VERSION >= 510 */
      }
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
  }  // end of buildMUMPSSolverGenerator

  void declareMUMPSLinearSolvers(LinearSolverFactory<true>& f) {
    declare(f, "MUMPSSolver",
            LinearSolverDescription{
                .library = "MUMPS",
                .kind = LinearSolverKind::DIRECT,
                .placement = LinearSolverPlacement::HOST,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = getMUMPSSolverOptions(),
                .description = "MUMPS parallel direct solver"},
            buildMUMPSSolverGenerator());
  }  // end of declareMUMPSLinearSolvers

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS) */
