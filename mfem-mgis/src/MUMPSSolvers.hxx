/*!
 * \file   src/MUMPSSolvers.hxx
 * \brief  This file declares the MUMPS solver which reuses its analysis for
 * the matrices of the same pattern, built by the factory of linear solvers
 * (internal header)
 * \date   12/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_MUMPSSOLVERS_HXX
#define LIB_MFEM_MGIS_SRC_MUMPSSOLVERS_HXX

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS)

#include <span>
#include <vector>
#include <optional>
#include <type_traits>
#include <mpi.h>
#include "mfem/linalg/hypre.hpp"
#include "mfem/linalg/mumps.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief MUMPS solver which reuses the analysis (ordering and symbolic
   * factorization) of the previous matrix when the new one has the same
   * pattern on all processes, as the matrices of the successive iterations
   * of the non linear solvers.
   *
   * `mfem::MUMPSSolver::SetReorderingReuse` (MFEM 4.10) does not check the
   * pattern of the new matrix, which would be factorized with the analysis
   * of the previous one. The pattern of the local part of the last
   * analysed matrix (rows, column indices of the diagonal and off-diagonal
   * blocks, global indices of the off-diagonal columns) is kept and
   * compared to the one of each new matrix. The analysis is collective: a
   * new analysis is made on all processes if the pattern changed on one of
   * them.
   */
  struct MUMPSSolverReusingAnalysis final : public mfem::MUMPSSolver {
    /*!
     * \brief constructor
     * \param[in] c: communicator
     */
    explicit MUMPSSolverReusingAnalysis(MPI_Comm);
    //! \brief set the operator
    void SetOperator(const mfem::Operator&) override;
    //! \return the number of analyses made
    [[nodiscard]] size_type getNumberOfAnalyses() const noexcept;

   private:
    //! \brief pattern of the local part of a matrix, viewed or copied
    template <template <typename> typename Storage>
    struct Pattern {
      //! \brief number of rows of the matrix
      HYPRE_BigInt global_number_of_rows = 0;
      //! \brief global index of the first local row
      HYPRE_BigInt first_row = 0;
      //! \brief global index of the first column of the diagonal block
      HYPRE_BigInt first_column = 0;
      //! \brief offsets of the rows of the diagonal block
      Storage<const HYPRE_Int> diagonal_rows;
      //! \brief column indices of the diagonal block
      Storage<const HYPRE_Int> diagonal_columns;
      //! \brief offsets of the rows of the off-diagonal block
      Storage<const HYPRE_Int> off_diagonal_rows;
      //! \brief column indices of the off-diagonal block
      Storage<const HYPRE_Int> off_diagonal_columns;
      //! \brief global indices of the columns of the off-diagonal block
      Storage<const HYPRE_BigInt> off_diagonal_column_map;
    };
    //! \brief storage of a view
    template <typename ValueType>
    using View = std::span<ValueType>;
    //! \brief storage of a copy
    template <typename ValueType>
    using Copy = std::vector<std::remove_const_t<ValueType>>;
    /*!
     * \return a view of the pattern of the local part of a matrix, whose
     * values must be readable on the host
     * \param[in] A: matrix
     */
    static Pattern<View> view(const mfem::HypreParMatrix&);
    //! \return a copy of a pattern
    static Pattern<Copy> copy(const Pattern<View>&);
    //! \return true if a copied pattern is equal to a viewed one
    static bool equal(const Pattern<Copy>&, const Pattern<View>&) noexcept;
    //! \brief communicator
    const MPI_Comm communicator;
    //! \brief pattern of the last analysed matrix, empty if none
    std::optional<Pattern<Copy>> pattern;
    //! \brief number of analyses made
    size_type number_of_analyses = 0;
  };  // end of struct MUMPSSolverReusingAnalysis

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS) */

#endif /* LIB_MFEM_MGIS_SRC_MUMPSSOLVERS_HXX */
