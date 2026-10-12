/*!
 * \file   src/MUMPSSolvers.cxx
 * \brief  This file implements the MUMPS solver which reuses its analysis
 * for the matrices of the same pattern
 * \date   12/10/2026
 */

#include "MUMPSSolvers.hxx"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS)

#include <utility>

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

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS) */
