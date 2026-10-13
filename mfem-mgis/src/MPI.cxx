/*!
 * \file   src/MPI.cxx
 * \brief  This file implements the functions declared in `MFEMMGIS/MPI.hxx`
 * \author Thomas Helfer
 * \date   25/02/2026
 */

#include <string>
#include <tuple>
#include "mfem/config/config.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/MPI.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI
  bool isTrueOnAllProcesses(const MPI_Comm& c, const bool b) noexcept {
    auto r = b;
    MPI_Allreduce(MPI_IN_PLACE, &r, 1, MPI_CXX_BOOL, MPI_LAND, c);
    return r;
  }

  bool isMPIMainProcess(const MPI_Comm& c) noexcept {
    return getMPIRank(c) == 0;
  }  // end of isMPIMainProcess

  int getMPIRank(const MPI_Comm& c) noexcept {
    int rank = 0;
#ifdef MFEM_USE_MPI
    MPI_Comm_rank(c, &rank);
#endif /* MFEM_USE_MPI */
    return rank;
  }  // end of getMPIRank

  int getMPISize(const MPI_Comm& c) noexcept {
    int size = 1;
#ifdef MFEM_USE_MPI
    MPI_Comm_size(c, &size);
#endif /* MFEM_USE_MPI */
    return size;
  }  // end of getMPISize

#endif /* MFEM_USE_MPI */

  bool isTrueOnAllProcesses(const MeshDiscretization& m,
                            const bool b) noexcept {
    if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return isTrueOnAllProcesses(getMPICommunicator(m), b);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return b;
  }  // end of isTrueOnAllProcesses

  bool isMPIMainProcess(const MeshDiscretization& m) noexcept {
    if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return isMPIMainProcess(getMPICommunicator(m));
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return true;
  }  // end of isMPIMainProcess

  int getMPIRank(const MeshDiscretization& m) noexcept {
    if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return getMPIRank(getMPICommunicator(m));
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return 0;
  }  // end of getMPIRank

  int getMPISize(const MeshDiscretization& m) noexcept {
    if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return getMPISize(getMPICommunicator(m));
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return 1;
  }  // end of getMPISize

  void raiseOnFailureOnAnyProcess(const MeshDiscretization& m,
                                  Context& ctx,
                                  const bool success,
                                  std::string_view what) {
    if (isTrueOnAllProcesses(m, success)) {
      return;
    }
    // the error messages of the processes which failed, given to all the
    // processes (in a new context: the context of a process which succeeded
    // may hold messages of failures already handled)
    auto failures = Context{};
    if (!success) {
      std::ignore = failures.registerErrorMessageWithoutSourceLocation(
          ctx.getErrorMessage());
    }
    synchronize(failures, m);
    const auto message = failures.getErrorMessage();
    throw(CollectiveFailure(
        std::string{what} + ": " +
        (message.empty() ? "failure without message" : message)));
  }  // end of raiseOnFailureOnAnyProcess

  CollectiveFailure::~CollectiveFailure() noexcept = default;

}  // end of namespace mfem_mgis
