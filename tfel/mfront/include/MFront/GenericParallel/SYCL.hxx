/*!
 * \file   mfront/include/MFront/GenericParallel/SYCL.hxx
 * \brief  This file declares the primitive of the `sycl` backend of the
 * `generic-parallel` interfaces.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_SYCL_HXX
#define LIB_MFRONT_GENERICPARALLEL_SYCL_HXX

#include <cstddef>
#include <utility>
#include <sycl/sycl.hpp>
#include "MFront/GenericParallel/Batch.hxx"
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallel/StatusBlockPool.hxx"
#include "MFront/GenericParallel/SyclBackendExecutionPolicy.hxx"

namespace mfront::generic_parallel::sycl {

  namespace internals {

    //! \brief atomic operations of the statistics
    struct Atomics {
      //! \brief a simple alias
      using AtomicCounter =
          ::sycl::atomic_ref<unsigned long long,
                             ::sycl::memory_order::relaxed,
                             ::sycl::memory_scope::device,
                             ::sycl::access::address_space::global_space>;
      static void add(unsigned long long* const p, const unsigned long long v) {
        AtomicCounter{*p}.fetch_add(v);
      }
      static void min(unsigned long long* const p, const unsigned long long v) {
        AtomicCounter{*p}.fetch_min(v);
      }
    };

    /*!
     * \brief kernel applying a function at one point and updating the
     * statistics of the batch.
     *
     * The statistics are updated by atomic operations, which are only
     * performed if an event occurs at this point (see `updateStatistics`).
     * Parallel reductions are not used since their synchronisation barriers
     * make the compilation of the large kernels of some behaviours fail with
     * AdaptiveCpp on CPUs.
     *
     * \tparam Statistics: statistics updated by the devices
     * \tparam PointFunction: function object applied at one point
     */
    template <typename Statistics, typename PointFunction>
    struct Kernel {
      void operator()(const ::sycl::id<1> i) const {
        updateStatistics<Atomics>(*(this->statistics), this->f(i[0]),
                                  this->rdt0);
      }
      //! \brief function object applied at one point
      PointFunction f;
      //! \brief statistics, accessible by the device
      Statistics* statistics;
      //! \brief ordered key of the initial time step scaling factor
      unsigned long long rdt0;
    };

    /*!
     * \return the pool of blocks of shared memory holding the statistics
     * of the synchronous computations, shared by all the functions of the
     * library. The blocks are identified by the device and the context of
     * the queue.
     * \tparam Statistics: statistics updated by the devices
     */
    template <typename Statistics>
    auto& getStatisticsBlocks() noexcept {
      using Pool = StatusBlockPool<Statistics,
                                   std::pair<::sycl::device, ::sycl::context>>;
      static auto pool = Pool{};
      return pool;
    }  // end of getStatisticsBlocks

  }  // end of namespace internals

  /*!
   * \return the policy given by the calling code, checked against the
   * policy of the backend (see
   * `mfront::generic_parallel::checkBackendExecutionPolicy`). The `sycl`
   * backend has no default policy: a null policy is refused.
   * \param[in] p: policy given by the calling code
   * \param[in] has_statistics: true if the calling code gave statistics
   */
  inline CheckedBackendExecutionPolicy<SyclBackendExecutionPolicy>
  checkExecutionPolicy(const mfront_gp_AbstractBackendExecutionPolicy* const p,
                       const bool has_statistics) noexcept {
    auto r = checkBackendExecutionPolicy<SyclBackendExecutionPolicy>(
        p, has_statistics);
    if ((r.error == nullptr) && (r.policy == nullptr)) {
      r.error = "no execution policy given: the sycl backend requires a queue";
    }
    return r;
  }  // end of checkExecutionPolicy

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * If the computations are asynchronous, this function returns once they
   * are submitted to the queue. Otherwise, it waits for their completion.
   *
   * \tparam StatisticsHandler: handler of the statistics (see
   * `MFront/GenericParallel/Batch.hxx`)
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] e: execution of the batch. The policy must be a
   * `SyclBackendExecutionPolicy`.
   * \param[in] h: handler of the statistics
   * \param[in] f: function
   * \param[in] device_stack_size: minimal stack size, in bytes, of the
   * threads of the device when the computations are performed on an NVIDIA
   * GPU through the CUDA runtime (see `setMinimalCUDAStackSize`). The stack
   * size is not modified if this value is null.
   * \return the status
   */
  template <typename StatisticsHandler, typename PointFunction>
  int apply(typename StatisticsHandler::OutputStatus& s,
            const size_type n,
            const BatchExecution<typename StatisticsHandler::Statistics>& e,
            const StatisticsHandler& h,
            const PointFunction& f,
            const std::size_t device_stack_size = 0) {
    using Statistics = typename StatisticsHandler::Statistics;
    using Kernel = internals::Kernel<Statistics, PointFunction>;
    if (n == 0) {
      return s.status;
    }
    const auto p = checkExecutionPolicy(e.policy, e.statistics != nullptr);
    if (p.error != nullptr) {
      return report(s, -5, p.error);
    }
    if (!setMinimalCUDAStackSize(device_stack_size,
                                 p.common_options.preserve_stack_size)) {
      return report(s, -7, stack_size_error_message);
    }
    auto q = p.policy->queue;
    // the kernel must be executed after the computations already submitted
    // to the queue, which may compute the inputs. This is guaranteed by an
    // in-order queue. With an out-of-order queue, the computations are
    // synchronous and the queue is synchronized before the kernel.
    if (!q.is_in_order()) {
      if (p.common_options.asynchronous) {
        return report(s, -5,
                      "asynchronous computations require an in-order SYCL "
                      "queue");
      }
      q.wait_and_throw();
    }
    const auto range = ::sycl::range<1>{n};
    const auto rdt0 = toOrderedKey(h.initial().rdt);
    if (e.statistics != nullptr) {
      auto event = q.parallel_for(range, Kernel{f, e.statistics, rdt0});
      if (!p.common_options.asynchronous) {
        event.wait_and_throw();
      }
      return s.status;
    }
    // synchronous computations: the statistics are held by a block of
    // shared memory, initialized and read by the host. The device only
    // accesses this block if an event occurs (see internals::Kernel): no
    // data is transferred otherwise.
    const auto block = internals::getStatisticsBlocks<Statistics>().acquire(
        {q.get_device(), q.get_context()},
        [&q] { return ::sycl::malloc_shared<Statistics>(1, q); });
    auto* const statistics = block.get();
    if (statistics == nullptr) {
      return report(s, -7, "unable to allocate the statistics of the batch");
    }
    h.initialize(*statistics);
    q.parallel_for(range, Kernel{f, statistics, rdt0}).wait_and_throw();
    return h.report(s, *statistics);
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
      return sycl::apply(s, n, e, h, f, device_stack_size);
    };
  }  // end of makeApply

}  // end of namespace mfront::generic_parallel::sycl

#endif /* LIB_MFRONT_GENERICPARALLEL_SYCL_HXX */
