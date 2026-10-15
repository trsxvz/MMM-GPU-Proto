/*!
 * \file   src/PETScSupport.hxx
 * \brief  This file declares the functions supporting the solvers and the
 * preconditioners of PETSc (internal header)
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_PETSCSUPPORT_HXX
#define LIB_MFEM_MGIS_SRC_PETSCSUPPORT_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief initialize PETSc, if not done yet by mfem-mgis, for its solvers
   * and its preconditioners, PETSc being finalized at the end of the
   * execution of mfem-mgis. This does not select the non linear solver of
   * PETSc (see `usePETSc`). This function shall not be called if the program
   * initialized PETSc itself.
   */
  void initializePETScOnDemand();

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_PETSCSUPPORT_HXX */
