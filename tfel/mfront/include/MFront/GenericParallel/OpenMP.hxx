/*!
 * \file   mfront/include/MFront/GenericParallel/OpenMP.hxx
 * \brief  This file declares the primitive of the `openmp` backend of the
 * `generic-parallel` interfaces, on the host or on a device (target
 * offloading).
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_OPENMP_HXX
#define LIB_MFRONT_GENERICPARALLEL_OPENMP_HXX

#include <cstddef>
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallel/OpenMPBackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/Batch.hxx"

namespace mfront::generic_parallel::openmp {

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * If the computations are offloaded, the fields of the batch must be
   * accessible from the device, for instance allocated in unified or
   * managed memory (`nvc++ -gpu=mem:managed`). This function waits for the
   * completion of the computations. The statistics given by the calling
   * code, if any, are accumulated on the host.
   *
   * \tparam offload: if true, the computations are offloaded to the default
   * device (`target` construct). Otherwise, they are performed by the
   * threads of the host (`parallel` construct) and the code does not depend
   * on the support of offloading by the OpenMP runtime.
   * \tparam StatisticsHandler: handler of the statistics (see
   * `MFront/GenericParallel/Batch.hxx`)
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] e: execution of the batch. The policy is a
   * `OpenMPBackendExecutionPolicy`, or null for the default policy: only its
   * common options are used, the computations being synchronous.
   * \param[in] h: handler of the statistics
   * \param[in] f: function
   * \param[in] device_stack_size: minimal stack size, in bytes, of the
   * threads of the device when the computations are offloaded to an NVIDIA
   * GPU through the CUDA runtime (see
   * `mfront::generic_parallel::setMinimalCUDAStackSize`). The stack size is
   * not modified if this value is null.
   * \return the status
   */
  template <bool offload, typename StatisticsHandler, typename PointFunction>
  int apply(typename StatisticsHandler::OutputStatus& s,
            const size_type n,
            const BatchExecution<typename StatisticsHandler::Statistics>& e,
            const StatisticsHandler& h,
            const PointFunction& f,
            const std::size_t device_stack_size = 0) {
    if (n == 0) {
      return s.status;
    }
    const auto p = checkBackendExecutionPolicy<OpenMPBackendExecutionPolicy>(
        e.policy, e.statistics != nullptr);
    if (p.error != nullptr) {
      return report(s, -5, p.error);
    }
    // the statistics are reduced member by member
    const auto r0 = h.initial();
    auto number_of_failures = r0.number_of_failures;
    auto first_failure = r0.first_failure;
    auto number_of_unreliable_results = r0.number_of_unreliable_results;
    auto first_unreliable_result = r0.first_unreliable_result;
    auto rdt = r0.rdt;
    if constexpr (offload) {
      if (!setMinimalCUDAStackSize(device_stack_size,
                                   p.common_options.preserve_stack_size)) {
        return report(s, -7, stack_size_error_message);
      }
#pragma omp target teams distribute parallel for map(to : f)        \
    reduction(+ : number_of_failures, number_of_unreliable_results) \
    reduction(min : first_failure, first_unreliable_result, rdt)
      for (size_type i = 0; i < n; ++i) {
        const auto r = f(i);
        number_of_failures += r.number_of_failures;
        first_failure =
            r.first_failure < first_failure ? r.first_failure : first_failure;
        number_of_unreliable_results += r.number_of_unreliable_results;
        first_unreliable_result =
            r.first_unreliable_result < first_unreliable_result
                ? r.first_unreliable_result
                : first_unreliable_result;
        rdt = r.rdt < rdt ? r.rdt : rdt;
      }
    } else {
      static_cast<void>(device_stack_size);
#pragma omp parallel for reduction(+ : number_of_failures,           \
                                       number_of_unreliable_results) \
    reduction(min : first_failure, first_unreliable_result, rdt)
      for (size_type i = 0; i < n; ++i) {
        const auto r = f(i);
        number_of_failures += r.number_of_failures;
        first_failure =
            r.first_failure < first_failure ? r.first_failure : first_failure;
        number_of_unreliable_results += r.number_of_unreliable_results;
        first_unreliable_result =
            r.first_unreliable_result < first_unreliable_result
                ? r.first_unreliable_result
                : first_unreliable_result;
        rdt = r.rdt < rdt ? r.rdt : rdt;
      }
    }
    return h.report(
        s, e.statistics,
        decltype(r0){
            .number_of_failures = number_of_failures,
            .first_failure = first_failure,
            .number_of_unreliable_results = number_of_unreliable_results,
            .first_unreliable_result = first_unreliable_result,
            .rdt = rdt});
  }  // end of apply

  /*!
   * \return the primitive of the backend, called by the functions of the
   * interfaces (see for instance `mfront::gpb::integrateBatch`)
   * \tparam offload: if true, the computations are offloaded
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   */
  template <bool offload>
  auto makeApply(const std::size_t device_stack_size = 0) {
    return [device_stack_size](auto& s, const size_type n, const auto& e,
                               const auto& h, const auto& f) {
      return openmp::apply<offload>(s, n, e, h, f, device_stack_size);
    };
  }  // end of makeApply

}  // end of namespace mfront::generic_parallel::openmp

#endif /* LIB_MFRONT_GENERICPARALLEL_OPENMP_HXX */
