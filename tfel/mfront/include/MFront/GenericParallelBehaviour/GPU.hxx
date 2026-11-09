/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/GPU.hxx
 * \brief  This file declares the functions of the `cuda` and `hip` backends
 * of the `generic-parallel` behaviour interface which are specific to the
 * behaviours. The primitive of the backends is declared in
 * `MFront/GenericParallel/GPU.hxx`, whose requirements apply.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_GPU_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_GPU_HXX

#include <stdexcept>
#include <type_traits>
#include "MFront/GenericParallel/GPU.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::gpu {

  /*!
   * \return a copy of a data in the memory of the current device
   * \param[in] v: data
   */
  template <typename T>
  const T* copyToDevice(const T& v) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "the shared data must be trivially copyable");
    constexpr auto success = MFRONT_GENERIC_PARALLEL_GPU_API(Success);
    T* p = nullptr;
    if (MFRONT_GENERIC_PARALLEL_GPU_API(Malloc)(reinterpret_cast<void**>(&p),
                                                sizeof(T)) != success) {
      throw std::runtime_error("memory allocation failed on the device");
    }
    if (MFRONT_GENERIC_PARALLEL_GPU_API(Memcpy)(
            p, &v, sizeof(T),
            MFRONT_GENERIC_PARALLEL_GPU_API(MemcpyHostToDevice)) != success) {
      static_cast<void>(MFRONT_GENERIC_PARALLEL_GPU_API(Free)(p));
      throw std::runtime_error("copy to the device failed");
    }
    return p;
  }  // end of copyToDevice

  /*!
   * \return a copy of a data shared by all the integration points (see the
   * `@SharedData` keyword) in the memory of the current device. A copy is
   * made at the first call on each device and is never freed.
   * \tparam get_shared_data: function returning the data
   */
  template <auto get_shared_data>
  const auto* getSharedDataCopy() {
    using T = std::remove_cvref_t<decltype(get_shared_data())>;
    static auto copies = SharedDataCopies<T, int>{};
    auto device = int{};
    if (MFRONT_GENERIC_PARALLEL_GPU_API(GetDevice)(&device) !=
        MFRONT_GENERIC_PARALLEL_GPU_API(Success)) {
      throw std::runtime_error("unable to get the current device");
    }
    return copies.get(device, [] { return copyToDevice(get_shared_data()); });
  }  // end of getSharedDataCopy

}  // end of namespace mfront::gpb::gpu

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_GPU_HXX */
