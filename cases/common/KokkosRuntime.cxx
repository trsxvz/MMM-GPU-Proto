/*!
 * \file   KokkosRuntime.cxx
 * \brief  Initialization of Kokkos (see KokkosRuntime.hxx), compiled by the
 * compiler of Kokkos (target Kokkos::kokkos).
 */

#include <Kokkos_Core.hpp>
#include "KokkosRuntime.hxx"

namespace proto {

  KokkosRuntime::KokkosRuntime() {
    if (!Kokkos::is_initialized()) {
      Kokkos::initialize();
      this->owner = true;
    }
  }  // end of KokkosRuntime

  KokkosRuntime::~KokkosRuntime() {
    if (this->owner && Kokkos::is_initialized() && !Kokkos::is_finalized()) {
      Kokkos::finalize();
    }
  }  // end of ~KokkosRuntime

}  // end of namespace proto
