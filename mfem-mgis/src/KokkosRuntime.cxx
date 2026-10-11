/*!
 * \file   src/KokkosRuntime.cxx
 * \brief  This file implements the initialization and the finalization of
 * Kokkos (see `KokkosRuntime.hxx`)
 * \date   11/10/2026
 * \note   the only source of mfem-mgis including Kokkos: it is compiled by
 * the compiler of Kokkos (target `MFEMMGISKokkos`)
 */

#include <stdexcept>
#include <Kokkos_Core.hpp>
#include "KokkosRuntime.hxx"

namespace mfem_mgis::internals {

  bool initializeKokkos(const int gpu) {
    if (Kokkos::is_finalized()) {
      throw std::runtime_error{
          "initialize: Kokkos is finalized and can't be initialized again"};
    }
    if (Kokkos::is_initialized()) {
      return false;
    }
    auto settings = Kokkos::InitializationSettings{};
    if (gpu >= 0) {
      settings.set_device_id(gpu);
    }
    Kokkos::initialize(settings);
    return true;
  }  // end of initializeKokkos

  void finalizeKokkos() {
    if (Kokkos::is_initialized() && !Kokkos::is_finalized()) {
      Kokkos::finalize();
    }
  }  // end of finalizeKokkos

  bool isKokkosInitialized() noexcept {
    return Kokkos::is_initialized() && !Kokkos::is_finalized();
  }  // end of isKokkosInitialized

}  // end of namespace mfem_mgis::internals
