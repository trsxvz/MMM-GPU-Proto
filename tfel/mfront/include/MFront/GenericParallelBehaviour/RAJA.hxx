/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/RAJA.hxx
 * \brief  This file declares the functions of the `raja` backend of the
 * `generic-parallel` behaviour interface which are specific to the
 * behaviours. The primitive of the backend is declared in
 * `MFront/GenericParallel/RAJA.hxx`.
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_RAJA_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_RAJA_HXX

#include <stdexcept>
#include <type_traits>
#include "MFront/GenericParallel/RAJA.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::raja {

  /*!
   * \return a copy of a data in the memory of a resource (memory of the
   * device of the resource for the policies of CUDA and HIP), kept until the
   * end of the process
   * \param[in] resource: resource
   * \param[in] v: data
   */
  template <typename Resource, typename T>
  const T* copyToDevice(Resource& resource, const T& v) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "the shared data must be trivially copyable");
    auto* const p =
        resource.template allocate<T>(1, camp::resources::MemoryAccess::Device);
    if (p == nullptr) {
      throw std::runtime_error(
          "mfront::gpb::raja::copyToDevice: allocation failed");
    }
    resource.memcpy(p, &v, sizeof(T));
    resource.wait();
    return p;
  }  // end of copyToDevice

  /*!
   * \return a copy of a data shared by all the integration points (see the
   * `@SharedData` keyword) in the memory of the default resource of an
   * execution policy (memory of the current device for the policies of CUDA
   * and HIP). A copy is made at the first call on each device and kept until
   * the end of the process, as with the `cuda` and `hip` backends.
   * \tparam ExecutionPolicy: execution policy
   * \tparam get_shared_data: function returning the data
   */
  template <typename ExecutionPolicy, auto get_shared_data>
  const auto* getSharedDataCopy() {
    using T = std::remove_cvref_t<decltype(get_shared_data())>;
    static auto copies = SharedDataCopies<T, int>{};
    auto resource =
        mfront::generic_parallel::raja::getDefaultResource<ExecutionPolicy>();
    return copies.get(
        mfront::generic_parallel::raja::getDevice(resource),
        [&resource] { return copyToDevice(resource, get_shared_data()); });
  }  // end of getSharedDataCopy

}  // end of namespace mfront::gpb::raja

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_RAJA_HXX */
