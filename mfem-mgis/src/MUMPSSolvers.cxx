/*!
 * \file   src/MUMPSSolvers.cxx
 * \brief  This file implements the MUMPS solver which reuses its analysis
 * for the matrices of the same pattern
 * \date   12/10/2026
 */

#include "MUMPSSolvers.hxx"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS)

#include <utility>
#include <memory>
#include <string>
#include <vector>
#include "mfem/fem/pfespace.hpp"
#include "MFEMMGIS/Parameters.hxx"
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
  }  // end of buildMUMPSSolverGenerator

  void declareMUMPSLinearSolvers(LinearSolverFactory<true>& f) {
    declare(f, "MUMPSSolver",
            LinearSolverDescription{
                .library = "MUMPS",
                .kind = LinearSolverKind::DIRECT,
                .placement = LinearSolverPlacement::HOST,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Symmetric", "PositiveDefinite", "ReuseAnalysis"},
                .description = "MUMPS parallel direct solver"},
            buildMUMPSSolverGenerator());
  }  // end of declareMUMPSLinearSolvers

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS) */
