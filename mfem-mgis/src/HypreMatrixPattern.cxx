/*!
 * \file   src/HypreMatrixPattern.cxx
 * \brief  This file implements the pattern of a parallel matrix of hypre
 * \date   12/10/2026
 */

#include "HypreMatrixPattern.hxx"

#ifdef MFEM_USE_MPI

#include <span>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace mfem_mgis::internals {

  namespace {

    /*!
     * \brief index arrays of a block of the local part of a matrix, where
     * hypre stores them, and their sizes
     */
    struct BlockIndices {
      //! \brief offsets of the rows, null if not allocated
      const HYPRE_Int* rows = nullptr;
      //! \brief column indices, null if not allocated
      const HYPRE_Int* columns = nullptr;
      //! \brief size of the array of the offsets of the rows
      std::size_t number_of_row_offsets = 0;
      //! \brief size of the array of the column indices
      std::size_t number_of_entries = 0;
    };

    /*!
     * \return the index arrays of a block
     * \param[in] b: block
     * \param[in] host: if true, the arrays are readable on the host and
     * the number of entries is read in the offsets of the rows
     */
    BlockIndices getBlockIndices(const hypre_CSRMatrix* const b,
                                 const bool host) {
      auto r = BlockIndices{};
      const auto n = static_cast<std::size_t>(hypre_CSRMatrixNumRows(b));
      r.rows = hypre_CSRMatrixI(b);
      r.columns = hypre_CSRMatrixJ(b);
      if (r.rows == nullptr) {
        return r;
      }
      r.number_of_row_offsets = n + 1;
      const auto nnz = host ? r.rows[n] : hypre_CSRMatrixNumNonzeros(b);
      if ((r.columns != nullptr) && (nnz > 0)) {
        r.number_of_entries = static_cast<std::size_t>(nnz);
      }
      return r;
    }  // end of getBlockIndices

    //! \brief index arrays of the local part of a matrix
    struct MatrixIndices {
      //! \brief if true, the arrays are readable on the host
      bool host = true;
      //! \brief memory location of the arrays of the blocks
      HYPRE_MemoryLocation location = HYPRE_MEMORY_HOST;
      //! \brief diagonal block
      BlockIndices diagonal;
      //! \brief off-diagonal block
      BlockIndices off_diagonal;
      //! \brief global indices of the off-diagonal columns (on the host)
      std::span<const HYPRE_BigInt> off_diagonal_column_map;
    };

    /*!
     * \return the index arrays of the local part of a matrix
     * \param[in] A: matrix
     */
    MatrixIndices getMatrixIndices(const mfem::HypreParMatrix& A) {
      // arrays where hypre reads them
      A.HypreRead();
      const auto* const m = static_cast<hypre_ParCSRMatrix*>(A);
      const auto* const diagonal = hypre_ParCSRMatrixDiag(m);
      const auto* const off_diagonal = hypre_ParCSRMatrixOffd(m);
      auto r = MatrixIndices{};
      r.location = hypre_CSRMatrixMemoryLocation(diagonal);
      const auto l = hypre_GetActualMemLocation(r.location);
      r.host = (l == hypre_MEMORY_HOST) || (l == hypre_MEMORY_HOST_PINNED) ||
               (l == hypre_MEMORY_UNIFIED);
      r.diagonal = getBlockIndices(diagonal, r.host);
      r.off_diagonal = getBlockIndices(off_diagonal, r.host);
      const auto* const column_map = hypre_ParCSRMatrixColMapOffd(m);
      const auto number_of_columns = hypre_CSRMatrixNumCols(off_diagonal);
      if ((column_map != nullptr) && (number_of_columns > 0)) {
        r.off_diagonal_column_map = std::span<const HYPRE_BigInt>(
            column_map, static_cast<std::size_t>(number_of_columns));
      }
      return r;
    }  // end of getMatrixIndices

    /*!
     * \brief copy an index array to the host, into an array of its size
     * \param[out] d: destination
     * \param[in] s: source
     * \param[in] l: memory location of the source
     */
    void copyToHost(std::vector<HYPRE_Int>& d,
                    const HYPRE_Int* const s,
                    const HYPRE_MemoryLocation l) {
      if (!d.empty()) {
        hypre_TMemcpy(d.data(), s, HYPRE_Int, d.size(), HYPRE_MEMORY_HOST, l);
      }
    }  // end of copyToHost

    /*!
     * \return true if an array of the recorded pattern has the given size
     * and values
     * \param[in] recorded: recorded array
     * \param[in] values: values, readable on the host
     * \param[in] n: number of values
     */
    template <typename ValueType>
    bool equal(const std::vector<ValueType>& recorded,
               const ValueType* const values,
               const std::size_t n) {
      if (recorded.size() != n) {
        return false;
      }
      return (n == 0) || std::equal(recorded.begin(), recorded.end(), values);
    }  // end of equal

  }  // end of namespace

  bool HypreMatrixPatternRecord::isSame(const mfem::HypreParMatrix& A) const {
    if (!this->pattern.has_value()) {
      return false;
    }
    const auto& p = *(this->pattern);
    const auto* const m = static_cast<hypre_ParCSRMatrix*>(A);
    if ((p.global_number_of_rows != hypre_ParCSRMatrixGlobalNumRows(m)) ||
        (p.first_row != hypre_ParCSRMatrixFirstRowIndex(m)) ||
        (p.first_column != hypre_ParCSRMatrixFirstColDiag(m))) {
      return false;
    }
    const auto i = getMatrixIndices(A);
    if (!equal(p.off_diagonal_column_map, i.off_diagonal_column_map.data(),
               i.off_diagonal_column_map.size())) {
      return false;
    }
    // sizes of the index arrays
    if ((p.diagonal_rows.size() != i.diagonal.number_of_row_offsets) ||
        (p.diagonal_columns.size() != i.diagonal.number_of_entries) ||
        (p.off_diagonal_rows.size() != i.off_diagonal.number_of_row_offsets) ||
        (p.off_diagonal_columns.size() != i.off_diagonal.number_of_entries)) {
      return false;
    }
    if (!i.host) {
      // copies of the index arrays into the buffers, which have the sizes
      // of the recorded arrays when the pattern was recorded from the device
      if ((this->buffers.diagonal_rows.size() != p.diagonal_rows.size()) ||
          (this->buffers.diagonal_columns.size() !=
           p.diagonal_columns.size()) ||
          (this->buffers.off_diagonal_rows.size() !=
           p.off_diagonal_rows.size()) ||
          (this->buffers.off_diagonal_columns.size() !=
           p.off_diagonal_columns.size())) {
        // pattern recorded from the host: a new analysis, conservatively
        return false;
      }
      copyToHost(this->buffers.diagonal_rows, i.diagonal.rows, i.location);
      copyToHost(this->buffers.diagonal_columns, i.diagonal.columns,
                 i.location);
      copyToHost(this->buffers.off_diagonal_rows, i.off_diagonal.rows,
                 i.location);
      copyToHost(this->buffers.off_diagonal_columns, i.off_diagonal.columns,
                 i.location);
      return (p.diagonal_rows == this->buffers.diagonal_rows) &&
             (p.diagonal_columns == this->buffers.diagonal_columns) &&
             (p.off_diagonal_rows == this->buffers.off_diagonal_rows) &&
             (p.off_diagonal_columns == this->buffers.off_diagonal_columns);
    }
    return equal(p.diagonal_rows, i.diagonal.rows,
                 i.diagonal.number_of_row_offsets) &&
           equal(p.diagonal_columns, i.diagonal.columns,
                 i.diagonal.number_of_entries) &&
           equal(p.off_diagonal_rows, i.off_diagonal.rows,
                 i.off_diagonal.number_of_row_offsets) &&
           equal(p.off_diagonal_columns, i.off_diagonal.columns,
                 i.off_diagonal.number_of_entries);
  }  // end of isSame

  bool isTrueOnAllProcesses(const bool b, const MPI_Comm c) {
    auto r = static_cast<int>(b);
    MPI_Allreduce(MPI_IN_PLACE, &r, 1, MPI_INT, MPI_LAND, c);
    return r != 0;
  }  // end of isTrueOnAllProcesses

  bool HypreMatrixPatternRecord::isSameOnAllProcesses(
      const mfem::HypreParMatrix& A, const MPI_Comm c, const bool b) const {
    return isTrueOnAllProcesses(b && this->isSame(A), c);
  }  // end of isSameOnAllProcesses

  void HypreMatrixPatternRecord::record(const mfem::HypreParMatrix& A) {
    const auto i = getMatrixIndices(A);
    const auto* const m = static_cast<hypre_ParCSRMatrix*>(A);
    auto copy = [&i](std::vector<HYPRE_Int>& d, const HYPRE_Int* const s,
                     const std::size_t n) {
      d.resize(n);
      if (i.host) {
        std::copy(s, s + n, d.begin());
      } else {
        copyToHost(d, s, i.location);
      }
    };
    auto p = HypreMatrixPattern{};
    p.global_number_of_rows = hypre_ParCSRMatrixGlobalNumRows(m);
    p.first_row = hypre_ParCSRMatrixFirstRowIndex(m);
    p.first_column = hypre_ParCSRMatrixFirstColDiag(m);
    copy(p.diagonal_rows, i.diagonal.rows, i.diagonal.number_of_row_offsets);
    copy(p.diagonal_columns, i.diagonal.columns, i.diagonal.number_of_entries);
    copy(p.off_diagonal_rows, i.off_diagonal.rows,
         i.off_diagonal.number_of_row_offsets);
    copy(p.off_diagonal_columns, i.off_diagonal.columns,
         i.off_diagonal.number_of_entries);
    p.off_diagonal_column_map.assign(i.off_diagonal_column_map.begin(),
                                     i.off_diagonal_column_map.end());
    // buffers of the comparisons, only needed for arrays on the device
    this->buffers = HypreMatrixPattern{};
    if (!i.host) {
      this->buffers.diagonal_rows.resize(p.diagonal_rows.size());
      this->buffers.diagonal_columns.resize(p.diagonal_columns.size());
      this->buffers.off_diagonal_rows.resize(p.off_diagonal_rows.size());
      this->buffers.off_diagonal_columns.resize(p.off_diagonal_columns.size());
    }
    this->pattern = std::move(p);
  }  // end of record

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_USE_MPI */
