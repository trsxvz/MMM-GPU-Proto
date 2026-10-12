/*!
 * \file   src/HypreMatrixPattern.hxx
 * \brief  This file declares the pattern of a parallel matrix of hypre,
 * compared by the direct solvers which reuse their analysis (internal
 * header)
 * \date   12/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_HYPREMATRIXPATTERN_HXX
#define LIB_MFEM_MGIS_SRC_HYPREMATRIXPATTERN_HXX

#include "mfem/config/config.hpp"

#ifdef MFEM_USE_MPI

#include <vector>
#include <optional>
#include <mpi.h>
#include "mfem/linalg/hypre.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief pattern of the local part of a parallel matrix of hypre: global
   * number of rows, first local row and column, rows and column indices of
   * the diagonal and off-diagonal blocks, global indices of the
   * off-diagonal columns. Two matrices with the same pattern on all
   * processes give the same coordinates of their coefficients to a direct
   * solver.
   */
  struct HypreMatrixPattern {
    //! \brief number of rows of the matrix
    HYPRE_BigInt global_number_of_rows = 0;
    //! \brief global index of the first local row
    HYPRE_BigInt first_row = 0;
    //! \brief global index of the first column of the diagonal block
    HYPRE_BigInt first_column = 0;
    //! \brief offsets of the rows of the diagonal block
    std::vector<HYPRE_Int> diagonal_rows;
    //! \brief column indices of the diagonal block
    std::vector<HYPRE_Int> diagonal_columns;
    //! \brief offsets of the rows of the off-diagonal block
    std::vector<HYPRE_Int> off_diagonal_rows;
    //! \brief column indices of the off-diagonal block
    std::vector<HYPRE_Int> off_diagonal_columns;
    //! \brief global indices of the columns of the off-diagonal block
    std::vector<HYPRE_BigInt> off_diagonal_column_map;
  };

  /*!
   * \brief pattern of the last matrix analysed by a direct solver, compared
   * to the ones of the next matrices.
   *
   * The index arrays of the matrices are read where hypre stores them: in
   * place on the host, or copied from the device, without their values,
   * into buffers allocated with the recorded pattern. Nothing is allocated
   * before the collective reduction of the comparison.
   */
  struct HypreMatrixPatternRecord {
    /*!
     * \return true if, on all the processes of a communicator, a pattern
     * is recorded, the local part of the given matrix has this pattern and
     * a local condition holds. Collective: all the processes take part in
     * the reduction, whatever their local results.
     * \param[in] A: matrix
     * \param[in] c: communicator
     * \param[in] b: local condition
     */
    [[nodiscard]] bool isSameOnAllProcesses(const mfem::HypreParMatrix&,
                                            MPI_Comm,
                                            const bool = true) const;
    /*!
     * \brief record the pattern of the local part of a matrix
     * \param[in] A: matrix
     */
    void record(const mfem::HypreParMatrix&);

   private:
    //! \return true if the local part of a matrix has the recorded pattern
    [[nodiscard]] bool isSame(const mfem::HypreParMatrix&) const;
    //! \brief recorded pattern, if any
    std::optional<HypreMatrixPattern> pattern;
    /*!
     * \brief host copies of the index arrays of a matrix stored on the
     * device, with the sizes of the recorded pattern
     */
    mutable HypreMatrixPattern buffers;
  };

  /*!
   * \return true if a condition holds on all the processes of a
   * communicator. Collective.
   * \param[in] b: local condition
   * \param[in] c: communicator
   */
  [[nodiscard]] bool isTrueOnAllProcesses(const bool, MPI_Comm);

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_USE_MPI */

#endif /* LIB_MFEM_MGIS_SRC_HYPREMATRIXPATTERN_HXX */
