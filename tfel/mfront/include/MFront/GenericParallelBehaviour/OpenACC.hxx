/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/OpenACC.hxx
 * \brief  This file declares the functions of the `openacc` backend of the
 * `generic-parallel` behaviour interface which are specific to the
 * behaviours. The primitive of the backend is declared in
 * `MFront/GenericParallel/OpenACC.hxx`.
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_OPENACC_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_OPENACC_HXX

#include <utility>
#include <stdexcept>
#include <type_traits>
#include <openacc.h>
#include "MFront/GenericParallel/OpenACC.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::openacc {

  /*!
   * \return a copy of a data shared by all the integration points (see the
   * `@SharedData` keyword) in the memory of the current device of OpenACC
   * (`acc_malloc`), which is also accessible when the memory of the device
   * is not unified with the one of the host. A copy is made at the first
   * call on each device and is never freed, as with the `cuda` and `hip`
   * backends.
   * \tparam get_shared_data: function returning the data
   */
  template <auto get_shared_data>
  const auto* getSharedDataCopy() {
    using T = std::remove_cvref_t<decltype(get_shared_data())>;
    static_assert(std::is_trivially_copyable_v<T>,
                  "the shared data must be trivially copyable");
    // device type and number of the current device
    using Device = std::pair<int, int>;
    static auto copies = SharedDataCopies<T, Device>{};
    const auto type = acc_get_device_type();
    const auto device =
        Device{static_cast<int>(type), acc_get_device_num(type)};
    return copies.get(device, [] {
      auto v = get_shared_data();
      auto* const p = static_cast<T*>(acc_malloc(sizeof(T)));
      if (p == nullptr) {
        throw std::runtime_error(
            "mfront::gpb::openacc::getSharedDataCopy: allocation failed");
      }
      acc_memcpy_to_device(p, &v, sizeof(T));
      return static_cast<const T*>(p);
    });
  }  // end of getSharedDataCopy

}  // end of namespace mfront::gpb::openacc

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_OPENACC_HXX */
