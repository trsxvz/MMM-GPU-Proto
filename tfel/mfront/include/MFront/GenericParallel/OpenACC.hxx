/*!
 * \file   mfront/include/MFront/GenericParallel/OpenACC.hxx
 * \brief  This file declares the primitive of the `openacc` backend of the
 * `generic-parallel` interfaces.
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_OPENACC_HXX
#define LIB_MFRONT_GENERICPARALLEL_OPENACC_HXX

#include <cstddef>
#include <openacc.h>
#include "MFront/GenericParallel/Batch.hxx"
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallel/OpenACCStackSize.hxx"
#include "MFront/GenericParallel/OpenACCBackendExecutionPolicy.hxx"

namespace mfront::generic_parallel::openacc {

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
   * calling code, if any, are accumulated on the host.
   *
   * \tparam StatisticsHandler: handler of the statistics (see
   * `MFront/GenericParallel/Batch.hxx`)
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] e: execution of the batch. The policy is a
   * `OpenACCBackendExecutionPolicy`, or null for the default policy: only its
   * common options are used, the computations being synchronous.
   * \param[in] h: handler of the statistics
   * \param[in] f: function
   * \param[in] device_stack_size: minimal stack size, in bytes, of the
   * threads of the device when the computations are performed by an NVIDIA
   * GPU (see `setMinimalStackSize`). The stack size is not modified if this
   * value is null.
   * \return the status
   */
  template <typename StatisticsHandler, typename PointFunction>
  int apply(typename StatisticsHandler::OutputStatus& s,
            const size_type n,
            const BatchExecution<typename StatisticsHandler::Statistics>& e,
            const StatisticsHandler& h,
            const PointFunction& f,
            const std::size_t device_stack_size = 0) {
    if (n == 0) {
      return s.status;
    }
    const auto p = checkBackendExecutionPolicy<OpenACCBackendExecutionPolicy>(
        e.policy, e.statistics != nullptr);
    if (p.error != nullptr) {
      return report(s, -5, p.error);
    }
    // the statistics are reduced member by member. The minimum reductions
    // of the keys of the first failure and of the first unreliable result
    // are performed on `size_type`, an alias of `size_t` of 64 bits: the
    // minimum reductions on `unsigned long long` are wrong on NVIDIA GPUs
    // with `nvc++` 25.11.
    static_assert(sizeof(size_type) == sizeof(unsigned long long));
    const auto r0 = h.initial();
    auto number_of_failures = r0.number_of_failures;
    auto first_failure = static_cast<size_type>(r0.first_failure);
    auto number_of_unreliable_results = r0.number_of_unreliable_results;
    auto first_unreliable_result =
        static_cast<size_type>(r0.first_unreliable_result);
    auto rdt = r0.rdt;
    if (!setMinimalStackSize(device_stack_size,
                             p.common_options.preserve_stack_size)) {
      return report(s, -7, stack_size_error_message);
    }
#pragma acc parallel loop copyin(f)                                 \
    reduction(+ : number_of_failures, number_of_unreliable_results) \
        reduction(min : first_failure, first_unreliable_result, rdt)
    for (size_type i = 0; i < n; ++i) {
      const auto r = f(i);
      const auto ff = static_cast<size_type>(r.first_failure);
      const auto fu = static_cast<size_type>(r.first_unreliable_result);
      number_of_failures += r.number_of_failures;
      first_failure = ff < first_failure ? ff : first_failure;
      number_of_unreliable_results += r.number_of_unreliable_results;
      first_unreliable_result =
          fu < first_unreliable_result ? fu : first_unreliable_result;
      rdt = r.rdt < rdt ? r.rdt : rdt;
    }
    return h.report(
        s, e.statistics,
        decltype(r0){
            .number_of_failures = number_of_failures,
            .first_failure = static_cast<unsigned long long>(first_failure),
            .number_of_unreliable_results = number_of_unreliable_results,
            .first_unreliable_result =
                static_cast<unsigned long long>(first_unreliable_result),
            .rdt = rdt});
  }  // end of apply

  /*!
   * \return the primitive of the backend, called by the functions of the
   * interfaces (see for instance `mfront::gpb::integrateBatch`)
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   */
  inline auto makeApply(const std::size_t device_stack_size = 0) {
    return [device_stack_size](auto& s, const size_type n, const auto& e,
                               const auto& h, const auto& f) {
      return openacc::apply(s, n, e, h, f, device_stack_size);
    };
  }  // end of makeApply

}  // end of namespace mfront::generic_parallel::openacc

#endif /* LIB_MFRONT_GENERICPARALLEL_OPENACC_HXX */
