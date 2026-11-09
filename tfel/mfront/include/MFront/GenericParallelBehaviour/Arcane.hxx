/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/Arcane.hxx
 * \brief  This file declares the functions of the `arcane` backend of the
 * `generic-parallel` behaviour interface which are specific to the
 * behaviours. The primitive of the backend is declared in
 * `MFront/GenericParallel/Arcane.hxx`.
 * \date   23/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_ARCANE_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_ARCANE_HXX

#include <cstring>
#include <stdexcept>
#include <type_traits>
#include "MFront/GenericParallel/Arcane.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::arcane {

  /*!
   * \return a copy of a data shared by all the integration points (see the
   * `@SharedData` keyword) in the unified memory of Arccore, accessible by
   * the host and by the devices whatever the queue. A copy is made at the
   * first call and is never freed.
   * \return the copy, or a null pointer if the policy is invalid (the call
   * reports this invalid argument before any computation)
   * \tparam get_shared_data: function returning the data
   * \param[in] policy: execution policy
   */
  template <auto get_shared_data>
  const auto* getSharedDataCopy(
      const mfront_gp_AbstractBackendExecutionPolicy* const policy) {
    namespace utilities = ::mfront::generic_parallel::arcane;
    using T = std::remove_cvref_t<decltype(get_shared_data())>;
    static_assert(std::is_trivially_copyable_v<T>,
                  "the shared data must be trivially copyable");
    static auto copies = SharedDataCopies<T, ::Arcane::IMemoryAllocator*>{};
    const auto p = utilities::checkExecutionPolicy(policy, true);
    if (p.error != nullptr) {
      return static_cast<const T*>(nullptr);
    }
    auto* const allocator = utilities::getUnifiedMemoryAllocator();
    return copies.get(allocator, [allocator]() -> const T* {
      auto* const c = utilities::allocate<T>(allocator);
      if (c == nullptr) {
        throw std::runtime_error(
            "mfront::gpb::arcane::getSharedDataCopy: allocation failed");
      }
      const auto v = get_shared_data();
      std::memcpy(c, &v, sizeof(T));
      return c;
    });
  }  // end of getSharedDataCopy

}  // end of namespace mfront::gpb::arcane

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_ARCANE_HXX */
