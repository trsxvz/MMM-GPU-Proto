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

#include <mpi.h>
#include "mfem/linalg/mumps.hpp"
#include "MFEMMGIS/Config.hxx"
#include "HypreMatrixPattern.hxx"

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
   * analysed matrix is kept and compared to the one of each new matrix. The
   * analysis is collective: a new analysis is made on all processes if the
   * pattern changed on one of them.
   */
  struct MUMPSSolverReusingAnalysis final : public mfem::MUMPSSolver {
    /*!
     * \brief constructor
     * \param[in] c: communicator
     */
    explicit MUMPSSolverReusingAnalysis(MPI_Comm c);
    //! \brief set the operator
    void SetOperator(const mfem::Operator&) override;
    //! \return the number of analyses made
    [[nodiscard]] size_type getNumberOfAnalyses() const noexcept;

   private:
    //! \brief communicator
    const MPI_Comm communicator;
    //! \brief pattern of the last analysed matrix, if any
    HypreMatrixPatternRecord pattern;
    //! \brief number of analyses made
    size_type number_of_analyses = 0;
  };  // end of struct MUMPSSolverReusingAnalysis

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS) */

#endif /* LIB_MFEM_MGIS_SRC_MUMPSSOLVERS_HXX */
