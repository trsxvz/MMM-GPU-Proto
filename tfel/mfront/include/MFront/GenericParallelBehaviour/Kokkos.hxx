/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/Kokkos.hxx
 * \brief  This file declares the functions of the `kokkos` backend of the
 * `generic-parallel` behaviour interface which are specific to the
 * behaviours. The primitive of the backend is declared in
 * `MFront/GenericParallel/Kokkos.hxx`.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_KOKKOS_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_KOKKOS_HXX

#include <stdexcept>
#include <type_traits>
#include "MFront/GenericParallel/Kokkos.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::kokkos {

  /*!
   * \return a copy of a data in the memory space of an execution space,
   * freed when Kokkos is finalized
   * \tparam ExecutionSpace: execution space
   * \param[in] v: data
   */
  template <typename ExecutionSpace, typename T>
  const T* copyToDevice(const T& v) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "the shared data must be trivially copyable");
    using MemorySpace = typename ExecutionSpace::memory_space;
    using Bytes = Kokkos::View<char*, MemorySpace, Kokkos::MemoryUnmanaged>;
    using HostBytes =
        Kokkos::View<const char*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
    if (!Kokkos::is_initialized()) {
      throw std::runtime_error("Kokkos is not initialized");
    }
    auto* const p = static_cast<T*>(Kokkos::kokkos_malloc<MemorySpace>(
        "mfront::gpb::kokkos::copyToDevice", sizeof(T)));
    Kokkos::deep_copy(Bytes(reinterpret_cast<char*>(p), sizeof(T)),
                      HostBytes(reinterpret_cast<const char*>(&v), sizeof(T)));
    Kokkos::push_finalize_hook([p] { Kokkos::kokkos_free<MemorySpace>(p); });
    return p;
  }  // end of copyToDevice

  /*!
   * \return a copy of a data shared by all the integration points (see the
   * `@SharedData` keyword) in the memory space of an execution space, made
   * at the first call (Kokkos uses one device per process) and freed when
   * Kokkos is finalized.
   * \tparam ExecutionSpace: execution space
   * \tparam get_shared_data: function returning the data
   */
  template <typename ExecutionSpace, auto get_shared_data>
  const auto* getSharedDataCopy() {
    static const auto* const p =
        copyToDevice<ExecutionSpace>(get_shared_data());
    return p;
  }  // end of getSharedDataCopy

}  // end of namespace mfront::gpb::kokkos

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_KOKKOS_HXX */
