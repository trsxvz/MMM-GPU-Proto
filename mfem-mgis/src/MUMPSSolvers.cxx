/*!
 * \file   src/MUMPSSolvers.cxx
 * \brief  This file implements the MUMPS solver which reuses its analysis
 * for the matrices of the same pattern
 * \date   12/10/2026
 */

#include "MUMPSSolvers.hxx"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS)

#include <utility>
#include <algorithm>

namespace mfem_mgis::internals {

  MUMPSSolverReusingAnalysis::MUMPSSolverReusingAnalysis(const MPI_Comm c)
      : mfem::MUMPSSolver(c),
        communicator(c) {}  // end of MUMPSSolverReusingAnalysis

  void MUMPSSolverReusingAnalysis::SetOperator(const mfem::Operator& op) {
    const auto* const A = dynamic_cast<const mfem::HypreParMatrix*>(&op);
    if (A == nullptr) {
      // error reported by MFEM
      mfem::MUMPSSolver::SetOperator(op);
      return;
    }
    // the comparison reads the matrix in place: nothing is allocated
    // before the collective reduction
    A->HostRead();
    auto same = static_cast<int>(this->pattern.has_value() &&
                                 equal(*(this->pattern), view(*A)));
    A->HypreRead();
    MPI_Allreduce(MPI_IN_PLACE, &same, 1, MPI_INT, MPI_LAND,
                  this->communicator);
    this->SetReorderingReuse(same != 0);
    // the pattern is forgotten until the factorization succeeds: if it
    // fails by an exception, the next matrix is analysed again
    auto previous = std::exchange(this->pattern, std::nullopt);
    mfem::MUMPSSolver::SetOperator(op);
    if (same != 0) {
      this->pattern = std::move(previous);
      return;
    }
    ++(this->number_of_analyses);
    A->HostRead();
    this->pattern = copy(view(*A));
    A->HypreRead();
  }  // end of SetOperator

  size_type MUMPSSolverReusingAnalysis::getNumberOfAnalyses() const noexcept {
    return this->number_of_analyses;
  }  // end of getNumberOfAnalyses

  MUMPSSolverReusingAnalysis::Pattern<MUMPSSolverReusingAnalysis::View>
  MUMPSSolverReusingAnalysis::view(const mfem::HypreParMatrix& A) {
    const auto* const m = static_cast<hypre_ParCSRMatrix*>(A);
    const auto* const diagonal = hypre_ParCSRMatrixDiag(m);
    const auto* const off_diagonal = hypre_ParCSRMatrixOffd(m);
    const auto n = hypre_CSRMatrixNumRows(diagonal);
    auto make_view = []<typename ValueType>(const ValueType* const values,
                                            const HYPRE_Int s) {
      if ((values == nullptr) || (s <= 0)) {
        return View<const ValueType>{};
      }
      return View<const ValueType>{values, static_cast<std::size_t>(s)};
    };
    auto p = Pattern<View>{};
    p.global_number_of_rows = hypre_ParCSRMatrixGlobalNumRows(m);
    p.first_row = hypre_ParCSRMatrixFirstRowIndex(m);
    p.first_column = hypre_ParCSRMatrixFirstColDiag(m);
    p.diagonal_rows = make_view(hypre_CSRMatrixI(diagonal), n + 1);
    p.off_diagonal_rows = make_view(hypre_CSRMatrixI(off_diagonal), n + 1);
    auto number_of_entries = [n](const View<const HYPRE_Int> rows) {
      return rows.empty() ? HYPRE_Int{} : rows[static_cast<std::size_t>(n)];
    };
    p.diagonal_columns = make_view(hypre_CSRMatrixJ(diagonal),
                                   number_of_entries(p.diagonal_rows));
    p.off_diagonal_columns = make_view(hypre_CSRMatrixJ(off_diagonal),
                                       number_of_entries(p.off_diagonal_rows));
    p.off_diagonal_column_map = make_view(hypre_ParCSRMatrixColMapOffd(m),
                                          hypre_CSRMatrixNumCols(off_diagonal));
    return p;
  }  // end of view

  MUMPSSolverReusingAnalysis::Pattern<MUMPSSolverReusingAnalysis::Copy>
  MUMPSSolverReusingAnalysis::copy(const Pattern<View>& v) {
    auto make_copy = []<typename ValueType>(const View<ValueType> values) {
      return Copy<ValueType>(values.begin(), values.end());
    };
    return {.global_number_of_rows = v.global_number_of_rows,
            .first_row = v.first_row,
            .first_column = v.first_column,
            .diagonal_rows = make_copy(v.diagonal_rows),
            .diagonal_columns = make_copy(v.diagonal_columns),
            .off_diagonal_rows = make_copy(v.off_diagonal_rows),
            .off_diagonal_columns = make_copy(v.off_diagonal_columns),
            .off_diagonal_column_map = make_copy(v.off_diagonal_column_map)};
  }  // end of copy

  bool MUMPSSolverReusingAnalysis::equal(const Pattern<Copy>& c,
                                         const Pattern<View>& v) noexcept {
    return (c.global_number_of_rows == v.global_number_of_rows) &&
           (c.first_row == v.first_row) && (c.first_column == v.first_column) &&
           std::ranges::equal(c.diagonal_rows, v.diagonal_rows) &&
           std::ranges::equal(c.diagonal_columns, v.diagonal_columns) &&
           std::ranges::equal(c.off_diagonal_rows, v.off_diagonal_rows) &&
           std::ranges::equal(c.off_diagonal_columns, v.off_diagonal_columns) &&
           std::ranges::equal(c.off_diagonal_column_map,
                              v.off_diagonal_column_map);
  }  // end of equal

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS) */
