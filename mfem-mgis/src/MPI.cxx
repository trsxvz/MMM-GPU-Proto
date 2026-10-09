/*!
 * \file   src/MPI.cxx
 * \brief  This file implements the functions declared in `MFEMMGIS/MPI.hxx`
 * \author Thomas Helfer
 * \date   25/02/2026
 */

#include <string>
#include "mfem/config/config.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/MPI.hxx"

namespace mfem_mgis {

  bool isTrueOnAllProcesses(const MeshDiscretization& m,
                            const bool b) noexcept {
    if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto r = b;
      MPI_Allreduce(MPI_IN_PLACE, &r, 1, MPI_CXX_BOOL, MPI_LAND,
                    getMPICommunicator(m));
      return r;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return b;
  }  // end of isTrueOnAllProcesses

  void raiseOnFailureOnAnyProcess(const MeshDiscretization& m,
                                  Context& ctx,
                                  const bool success,
                                  std::string_view what) {
    if (isTrueOnAllProcesses(m, success)) {
      return;
    }
    if (!success) {
      raise(std::string{what} + ": " + ctx.getErrorMessage());
    }
    raise(std::string{what} + ": failure on another process");
  }  // end of raiseOnFailureOnAnyProcess

}  // end of namespace mfem_mgis
