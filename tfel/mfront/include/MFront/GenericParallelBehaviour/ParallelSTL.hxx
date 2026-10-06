/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/ParallelSTL.hxx
 * \brief  This file implements the integration of a batch of integration
 * points with the parallel algorithms of the C++ standard library.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_PARALLELSTL_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_PARALLELSTL_HXX

#include <cstddef>
#include <iterator>
#include <numeric>
#include <execution>
#include "MFront/GenericParallel/CountingIterator.hxx"
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::stlpar {

  //! \brief a simple alias
  using CountingIterator =
      mfront::generic_parallel::CountingIterator<size_type>;

  static_assert(std::random_access_iterator<CountingIterator>);

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch. This function waits for the completion of the
   * computations. The statistics given by the execution options, if any,
   * are accumulated on the host.
   *
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \tparam ExecutionPolicy: execution policy
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] o: execution options. The queue is not used.
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \param[in] policy: execution policy
   * \param[in] device_stack_size: minimal stack size, in bytes, of the
   * threads of the device when the parallel algorithms are offloaded to an
   * NVIDIA GPU by `nvc++ -stdpar=gpu` (see
   * `mfront::generic_parallel::setMinimalCUDAStackSize`). The stack size is
   * not modified if this value is null.
   * \return the status
   */
  template <typename PointFunction, typename ExecutionPolicy>
  int apply(OutputStatus& s,
            const size_type n,
            const ExecutionOptions& o,
            const real rdt0,
            const PointFunction& f,
            ExecutionPolicy&& policy,
            const std::size_t device_stack_size = 0) {
    if (!mfront::generic_parallel::setMinimalCUDAStackSize(
            device_stack_size, o.preserve_stack_size != 0)) {
      return report(s, -3, mfront::generic_parallel::stack_size_error_message);
    }
    const auto r = std::transform_reduce(
        policy, CountingIterator{0}, CountingIterator{n},
        BatchStatistics{.rdt = rdt0}, CombineStatistics{}, f);
    return report(s, o, r);
  }  // end of apply

  /*!
   * \return the primitive of the backend, as expected by
   * `mfront::gpb::integrateBatch` and the similar functions of
   * `MFront/GenericParallelBehaviour/Integrate.hxx`, called by the generated
   * functions
   * \param[in,out] s: output status
   * \param[in] policy: execution policy
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   */
  template <typename ExecutionPolicy>
  auto makeApply(OutputStatus& s,
                 const ExecutionPolicy policy,
                 const std::size_t device_stack_size = 0) {
    return [&s, policy, device_stack_size](const size_type n,
                                           const ExecutionOptions& o,
                                           const real rdt0, const auto& f) {
      return stlpar::apply(s, n, o, rdt0, f, policy, device_stack_size);
    };
  }  // end of makeApply

}  // end of namespace mfront::gpb::stlpar

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_PARALLELSTL_HXX */
