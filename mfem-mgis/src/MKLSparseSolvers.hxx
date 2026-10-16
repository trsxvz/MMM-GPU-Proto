/*!
 * \file   src/MKLSparseSolvers.hxx
 * \brief  This file declares the interface of the sparse direct solvers of
 * oneMKL (see `MKLSparseSolvers.cxx`)
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_MKLSPARSESOLVERS_HXX
#define LIB_MFEM_MGIS_MKLSPARSESOLVERS_HXX

#include "mfem/config/config.hpp"

#if defined(MFEM_MGIS_HAVE_MKL_PARDISO) || \
    (defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_MKL_CPARDISO))

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::internals {

  //! \brief interface of the sparse direct solvers of oneMKL
  struct MKLSparseDirectSolverBase {
    //! \return the number of analyses of the matrices
    [[nodiscard]] virtual size_type getNumberOfAnalyses() const noexcept = 0;
    //! \brief destructor
    virtual ~MKLSparseDirectSolverBase();
  };  // end of struct MKLSparseDirectSolverBase

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_MGIS_HAVE_MKL_PARDISO) || (defined(MFEM_USE_MPI) && \
          defined(MFEM_MGIS_HAVE_MKL_CPARDISO)) */

#endif /* LIB_MFEM_MGIS_MKLSPARSESOLVERS_HXX */
