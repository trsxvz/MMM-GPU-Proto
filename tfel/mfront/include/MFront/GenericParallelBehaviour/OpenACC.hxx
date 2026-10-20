/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/OpenACC.hxx
 * \brief  This file implements the integration of a batch of integration
 * points with OpenACC.
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

#include <cstddef>
#include <utility>
#include <stdexcept>
#include <type_traits>
#include <openacc.h>
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallel/OpenACCStackSize.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::openacc {

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * The points are processed by a `parallel loop` construct of OpenACC, on
   * the device selected by the compiler and the OpenACC runtime (for
   * instance, a GPU with `nvc++ -acc=gpu`, the cores of the host with
   * `nvc++ -acc=multicore`). If the computations are performed by a GPU,
   * the fields of the batch must be accessible from it, for instance
   * allocated in unified or managed memory (`nvc++ -gpu=mem:managed`). This
   * function waits for the completion of the computations. The statistics
   * of the points are merged by a reduction, and the statistics given by the
   * execution options, if any, are accumulated on the host.
   *
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] o: execution options. The queue is not used.
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \param[in] device_stack_size: minimal stack size, in bytes, of the
   * threads of the device when the computations are performed by an NVIDIA
   * GPU (see `mfront::generic_parallel::openacc::setMinimalStackSize`). The
   * stack size is not modified if this value is null.
   * \return the status
   */
  template <typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            const ExecutionOptions& o,
            const real rdt0,
            const PointFunction& f,
            const std::size_t device_stack_size = 0) {
    // the statistics are reduced member by member. The minimum reduction of
    // the index of the first failure is performed on `size_type`, an alias
    // of `size_t`: the minimum reductions on `unsigned long long` are wrong
    // on NVIDIA GPUs with `nvc++` 25.11 (see the material properties).
    auto number_of_failures = size_type{};
    auto first_failure_index = no_failure_index;
    auto number_of_unreliable_results = size_type{};
    auto rdt = rdt0;
    if (!mfront::generic_parallel::openacc::setMinimalStackSize(
            device_stack_size, o.preserve_stack_size != 0)) {
      return report(s, -3, mfront::generic_parallel::stack_size_error_message);
    }
#pragma acc parallel loop copyin(f)                                 \
    reduction(+ : number_of_failures, number_of_unreliable_results) \
        reduction(min : first_failure_index, rdt)
    for (size_type i = 0; i < n; ++i) {
      const auto r = f(i);
      number_of_failures += r.number_of_failures;
      first_failure_index = r.first_failure_index < first_failure_index
                                ? r.first_failure_index
                                : first_failure_index;
      number_of_unreliable_results += r.number_of_unreliable_results;
      rdt = r.rdt < rdt ? r.rdt : rdt;
    }
    return report(s, o,
                  BatchStatistics{.number_of_failures = number_of_failures,
                                  .first_failure_index = first_failure_index,
                                  .number_of_unreliable_results =
                                      number_of_unreliable_results,
                                  .rdt = rdt});
  }  // end of apply

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

  /*!
   * \return the primitive of the backend, as expected by
   * `mfront::gpb::integrateBatch` and the similar functions of
   * `MFront/GenericParallelBehaviour/Integrate.hxx`, called by the generated
   * functions
   * \param[in,out] s: output status
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   */
  inline auto makeApply(OutputStatus& s,
                        const std::size_t device_stack_size = 0) {
    return [&s, device_stack_size](const size_type n, const ExecutionOptions& o,
                                   const real rdt0, const auto& f) {
      return openacc::apply(s, n, o, rdt0, f, device_stack_size);
    };
  }  // end of makeApply

}  // end of namespace mfront::gpb::openacc

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_OPENACC_HXX */
