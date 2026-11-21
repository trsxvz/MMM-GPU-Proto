/*!
 * \file   mfront/include/MFront/GenericParallel/RAJA.hxx
 * \brief  This file declares the primitive of the `raja` backend of the
 * `generic-parallel` interfaces.
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_RAJA_HXX
#define LIB_MFRONT_GENERICPARALLEL_RAJA_HXX

#include <new>
#include "RAJA/RAJA.hpp"
#include "MFront/GenericParallel/Batch.hxx"
#include "MFront/GenericParallel/RAJAExecutionPolicy.hxx"
#include "MFront/GenericParallel/RAJABackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/StatusBlockPool.hxx"

namespace mfront::generic_parallel::raja {

  namespace internals {

    //! \brief atomic operations of the statistics
    struct Atomics {
      RAJA_HOST_DEVICE static void add(unsigned long long* const p,
                                       const unsigned long long v) {
        RAJA::atomicAdd<RAJA::auto_atomic>(p, v);
      }
      RAJA_HOST_DEVICE static void min(unsigned long long* const p,
                                       const unsigned long long v) {
        RAJA::atomicMin<RAJA::auto_atomic>(p, v);
      }
    };

    /*!
     * \brief function object applying a function at one point and updating
     * the statistics of the batch by atomic operations, which are only
     * performed if an event occurs at this point (see `updateStatistics`).
     * \tparam Statistics: statistics updated by the devices
     * \tparam PointFunction: function object applied at one point
     */
    template <typename Statistics, typename PointFunction>
    struct AtomicUpdate {
      //! \brief apply the function at the `i`-th point
      RAJA_HOST_DEVICE void operator()(const size_type i) const {
        const auto r = this->f(i);
        if ((this->statistics != nullptr) &&
            (updateStatistics<Atomics>(*(this->statistics), r, this->rdt0))) {
          signalEvents(this->events);
        }
      }
      //! \brief function object applied at one point
      PointFunction f;
      /*!
       * \brief statistics, or null if the function never reports an event
       * (see the `hasEvents` method of the handlers of the statistics)
       */
      Statistics* statistics;
      //! \brief flag set when the statistics are updated (see
      //! `StatisticsBlock`), or null
      int* events;
      //! \brief ordered key of the initial time step scaling factor
      unsigned long long rdt0;
    };

    /*!
     * \return a block of statistics (see `StatisticsBlock`), allocated in
     * pinned host memory by a resource of camp, the statistics being
     * allocated in the memory of its device. The block is reset before its
     * first use. An exception is thrown if an allocation fails.
     * \tparam Statistics: statistics updated by the devices
     * \param[in] r: resource
     */
    template <typename Statistics, typename Resource>
    StatisticsBlock<Statistics>* allocateStatisticsBlock(Resource& r) {
      using camp::resources::MemoryAccess;
      auto* const statistics =
          r.template allocate<Statistics>(1, MemoryAccess::Device);
      try {
        auto* const block = r.template allocate<StatisticsBlock<Statistics>>(
            1, MemoryAccess::Pinned);
        return ::new (block)
            StatisticsBlock<Statistics>{.statistics = statistics, .events = 1};
      } catch (...) {
        r.deallocate(statistics, MemoryAccess::Device);
        throw;
      }
    }  // end of allocateStatisticsBlock

    /*!
     * \return the pool of blocks holding the statistics of the synchronous
     * computations with the resources of the given type (see
     * `StatisticsBlock`), by device, shared by all the functions of the
     * library. The blocks are kept until the end of the process.
     * \tparam Statistics: statistics updated by the devices
     * \tparam Resource: resource of camp
     */
    template <typename Statistics, typename Resource>
    auto& getStatisticsBlocks() noexcept {
      using Pool = StatusBlockPool<StatisticsBlock<Statistics>, int>;
      static auto pool = Pool{};
      return pool;
    }  // end of getStatisticsBlocks

  }  // end of namespace internals

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * The statistics are updated by atomic operations, only performed if an
   * event occurs at a point. If the calling code gives statistics, they are
   * updated, and this function returns once the computations are submitted
   * if the computations are asynchronous. Otherwise, the statistics are held
   * by a local variable with an execution policy of the host, and otherwise
   * by a block of a pool, the statistics being in the memory of the device
   * of the resource (see `StatisticsBlock`), as by the GPU backends: a call
   * without event performs neither copy nor additional synchronization. A
   * synchronous computation waits for the completion of the computations.
   *
   * \tparam ExecutionPolicy: execution policy of RAJA
   * \tparam StatisticsHandler: handler of the statistics (see
   * `MFront/GenericParallel/Batch.hxx`)
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] e: execution of the batch. The policy is a
   * `RAJABackendExecutionPolicy` holding a resource of camp matching the
   * execution policy (`camp::resources::Cuda`, `Hip` or `Host`). If null,
   * the default stream of the current device is used with CUDA and HIP, the
   * default resource otherwise (see `getDefaultResource`).
   * \param[in] h: handler of the statistics
   * \param[in] f: function
   * \return the status
   */
  template <typename ExecutionPolicy,
            typename StatisticsHandler,
            typename PointFunction>
  int apply(typename StatisticsHandler::OutputStatus& s,
            const size_type n,
            const BatchExecution<typename StatisticsHandler::Statistics>& e,
            const StatisticsHandler& h,
            const PointFunction& f) {
    using Statistics = typename StatisticsHandler::Statistics;
    using Update = internals::AtomicUpdate<Statistics, PointFunction>;
    using ResourceType = Resource<ExecutionPolicy>;
    using BackendExecutionPolicy = RAJABackendExecutionPolicy<ResourceType>;
    if (n == 0) {
      return s.status;
    }
    const auto p = checkBackendExecutionPolicy<BackendExecutionPolicy>(
        e.policy, e.statistics != nullptr);
    if (p.error != nullptr) {
      return report(s, -5, p.error);
    }
    const auto error = ensureStackSize<ExecutionPolicy>(
        f, p.common_options.preserve_stack_size);
    if (!error.empty()) {
      return report(s, -7, error);
    }
    auto resource = p.policy == nullptr ? getDefaultResource<ExecutionPolicy>()
                                        : p.policy->resource;
    const auto range = RAJA::TypedRangeSegment<size_type>(0, n);
    const auto rdt0 = toOrderedKey(h.initial().rdt);
    if ((e.statistics != nullptr) || (!h.hasEvents())) {
      // statistics given by the calling code: the result is given by the
      // report function of the library. Without event, the statistics, if
      // any, are not accessed by the kernel, as for the GPU backends
      RAJA::forall<ExecutionPolicy>(resource, range,
                                    Update{f, e.statistics, nullptr, rdt0});
      if (!p.common_options.asynchronous) {
        resource.wait();
      }
      return s.status;
    }
    // synchronous computations: the statistics are held by a local variable
    // with an execution policy of the host, and otherwise by a block of the
    // pool (see `StatisticsBlock`)
    auto statistics = Statistics{};
    h.initialize(statistics);
    if constexpr (!isDeviceExecutionPolicy<ExecutionPolicy>) {
      RAJA::forall<ExecutionPolicy>(resource, range,
                                    Update{f, &statistics, nullptr, rdt0});
      resource.wait();
    } else {
      const auto block =
          internals::getStatisticsBlocks<Statistics, ResourceType>().acquire(
              getDevice(resource), [&resource] {
                return internals::allocateStatisticsBlock<Statistics>(resource);
              });
      auto* const b = block.get();
      if (b == nullptr) {
        return report(s, -7, "unable to allocate the statistics of the batch");
      }
      if (holdsEvents(*b)) {
        // synchronous reset, after an event or at the first use
        const auto neutral = getNeutralStatistics(statistics);
        resource.memcpy(b->statistics, &neutral, sizeof(Statistics));
        resource.wait();
        clearEvents(*b);
      }
      RAJA::forall<ExecutionPolicy>(
          resource, range, Update{f, b->statistics, &(b->events), rdt0});
      resource.wait();
      if (holdsEvents(*b)) {
        auto block_statistics = Statistics{};
        resource.memcpy(&block_statistics, b->statistics, sizeof(Statistics));
        resource.wait();
        statistics = mergeBlockStatistics(statistics, block_statistics);
      }
    }
    return h.report(s, statistics);
  }  // end of apply

  /*!
   * \return the primitive of the backend, called by the functions of the
   * interfaces (see for instance `mfront::gpb::integrateBatch`)
   * \tparam ExecutionPolicy: execution policy
   */
  template <typename ExecutionPolicy>
  auto makeApply() {
    return [](auto& s, const size_type n, const auto& e, const auto& h,
              const auto& f) {
      return raja::apply<ExecutionPolicy>(s, n, e, h, f);
    };
  }  // end of makeApply

}  // end of namespace mfront::generic_parallel::raja

#endif /* LIB_MFRONT_GENERICPARALLEL_RAJA_HXX */
