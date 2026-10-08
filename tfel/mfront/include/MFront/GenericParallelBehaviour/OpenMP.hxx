/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/OpenMP.hxx
 * \brief  This file implements the integration of a batch of integration
 * points with OpenMP, on the host or on a device (target offloading).
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_OPENMP_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_OPENMP_HXX

#include <cstddef>
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::openmp {

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * If the computations are offloaded, the fields of the batch must be
   * accessible from the device, for instance allocated in unified or
   * managed memory (`nvc++ -gpu=mem:managed`). This function waits for the
   * completion of the computations. The statistics given by the execution
   * options, if any, are accumulated on the host.
   *
   * \tparam offload: if true, the computations are offloaded to the default
   * device (`target` construct). Otherwise, they are performed by the
   * threads of the host (`parallel` construct) and the code does not depend
   * on the support of offloading by the OpenMP runtime.
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] o: execution options. The queue is not used.
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \param[in] device_stack_size: minimal stack size, in bytes, of the
   * threads of the device when the computations are offloaded to an NVIDIA
   * GPU through the CUDA runtime (see
   * `mfront::generic_parallel::setMinimalCUDAStackSize`). The stack size is
   * not modified if this value is null.
   * \return the status
   */
  template <bool offload, typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            const ExecutionOptions& o,
            const real rdt0,
            const PointFunction& f,
            const std::size_t device_stack_size = 0) {
    // the statistics are reduced member by member
    auto number_of_failures = size_type{};
    auto first_failure_index = no_failure_index;
    auto number_of_unreliable_results = size_type{};
    auto rdt = rdt0;
    if constexpr (offload) {
      if (!mfront::generic_parallel::setMinimalCUDAStackSize(
              device_stack_size, o.preserve_stack_size != 0)) {
        return report(s, -3,
                      mfront::generic_parallel::stack_size_error_message);
      }
#pragma omp target teams distribute parallel for map(to : f)        \
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
    } else {
      static_cast<void>(device_stack_size);
#pragma omp parallel for reduction(+ : number_of_failures,           \
                                       number_of_unreliable_results) \
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
    }
    return report(s, o,
                  BatchStatistics{.number_of_failures = number_of_failures,
                                  .first_failure_index = first_failure_index,
                                  .number_of_unreliable_results =
                                      number_of_unreliable_results,
                                  .rdt = rdt});
  }  // end of apply

  /*!
   * \return the primitive of the backend, as expected by
   * `mfront::gpb::integrateBatch` and the similar functions of
   * `MFront/GenericParallelBehaviour/Integrate.hxx`, called by the generated
   * functions
   * \tparam offload: if true, the computations are offloaded
   * \param[in,out] s: output status
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   */
  template <bool offload>
  auto makeApply(OutputStatus& s, const std::size_t device_stack_size = 0) {
    return [&s, device_stack_size](const size_type n, const ExecutionOptions& o,
                                   const real rdt0, const auto& f) {
      return openmp::apply<offload>(s, n, o, rdt0, f, device_stack_size);
    };
  }  // end of makeApply

}  // end of namespace mfront::gpb::openmp

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_OPENMP_HXX */
