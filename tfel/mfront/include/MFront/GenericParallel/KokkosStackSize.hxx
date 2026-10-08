/*!
 * \file   mfront/include/MFront/GenericParallel/KokkosStackSize.hxx
 * \brief  This file declares a function adapting the stack size of the
 * threads of the device to a function applied by `Kokkos`, shared by the
 * `generic-parallel` interfaces for behaviours and for material
 * properties.
 * \date   08/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_KOKKOSSTACKSIZE_HXX
#define LIB_MFRONT_GENERICPARALLEL_KOKKOSSTACKSIZE_HXX

#include <string>
#include <type_traits>
#include <Kokkos_Core.hpp>
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#include "MFront/GenericParallel/GPUStackSize.hxx"
#endif

namespace mfront::generic_parallel::kokkos {

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)

  /*!
   * \brief kernel which is never launched: its attributes give the stack
   * size required to apply a function at one point. Kokkos does not
   * increase the stack size of the threads, which is too small by default
   * for some behaviours.
   */
  template <typename Function>
  __global__ void stackSizeProbe(const Function f) {
    f(0);
  }  // end of stackSizeProbe

#endif /* defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) */

  /*!
   * \brief increase the stack size of the threads of the device if the
   * given function, applied at one point in the given execution space,
   * requires it. Nothing is done if the execution space is not a CUDA or a
   * HIP execution space.
   * \return an error message, empty on success
   * \tparam ExecutionSpace: execution space
   * \param[in] f: function applied at one point
   * \param[in] preserve: if true, the stack size is not modified, and an
   * error is returned if it is too small
   */
  template <typename ExecutionSpace, typename Function>
  std::string ensureStackSize(const Function& f, const bool preserve) {
#if defined(KOKKOS_ENABLE_CUDA)
    constexpr auto device = std::is_same_v<ExecutionSpace, Kokkos::Cuda>;
#elif defined(KOKKOS_ENABLE_HIP)
    constexpr auto device = std::is_same_v<ExecutionSpace, Kokkos::HIP>;
#else
    constexpr auto device = false;
#endif
    static_cast<void>(f);
    static_cast<void>(preserve);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    if constexpr (device) {
      const auto r = gpu::ensureStackSize<stackSizeProbe<Function>>(preserve);
      if (!r.succeeded()) {
        return r.getErrorMessage();
      }
    }
#else
    static_cast<void>(device);
#endif
    return {};
  }  // end of ensureStackSize

}  // end of namespace mfront::generic_parallel::kokkos

#endif /* LIB_MFRONT_GENERICPARALLEL_KOKKOSSTACKSIZE_HXX */
