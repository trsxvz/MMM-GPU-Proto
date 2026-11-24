/*!
 * \file   mfront/include/MFront/GenericParallel/Arcane.hxx
 * \brief  This file declares the primitive of the `arcane` backend of the
 * `generic-parallel` interfaces, based on the accelerator API of Arccore.
 * \date   23/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_ARCANE_HXX
#define LIB_MFRONT_GENERICPARALLEL_ARCANE_HXX

#include <new>
#include <arccore/accelerator/Atomic.h>
#include <arccore/common/accelerator/Memory.h>
#include "MFront/GenericParallel/Batch.hxx"
#include "MFront/GenericParallel/StatusBlockPool.hxx"
#include "MFront/GenericParallel/ArcaneUtilities.hxx"
#include "MFront/GenericParallel/ArcaneBackendExecutionPolicy.hxx"

namespace mfront::generic_parallel::arcane {

  namespace internals {

    //! \brief atomic operations of the statistics, provided by Arccore
    struct Atomics {
      ARCCORE_HOST_DEVICE static void add(unsigned long long* const p,
                                          const unsigned long long v) {
        using ::Arcane::Accelerator::eAtomicOperation;
        ::Arcane::Accelerator::doAtomic<eAtomicOperation::Add>(p, v);
      }
      ARCCORE_HOST_DEVICE static void min(unsigned long long* const p,
                                          const unsigned long long v) {
        using ::Arcane::Accelerator::eAtomicOperation;
        ::Arcane::Accelerator::doAtomic<eAtomicOperation::Min>(p, v);
      }
      /*!
       * \brief set a flag (see `signalEvents`), by an addition: contrary to
       * the exchanges, the additions of the threads of a warp of an NVIDIA
       * device are merged in a single atomic operation by its compiler
       * \return true if the flag was already set
       */
      ARCCORE_HOST_DEVICE static bool testAndSet(unsigned long long* const p) {
        using ::Arcane::Accelerator::eAtomicOperation;
        return ::Arcane::Accelerator::doAtomic<eAtomicOperation::Add>(
                   p, 1ull) != 0;
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
      ARCCORE_HOST_DEVICE void operator()(const size_type i) const {
        const auto r = this->f(i);
        if (this->statistics != nullptr) {
          updateStatistics<Atomics>(*(this->statistics), r, this->rdt0,
                                    this->signalled, this->events);
        }
      }
      //! \brief function object applied at one point
      PointFunction f;
      /*!
       * \brief statistics, or null if the function never reports an event
       * (see the `hasEvents` method of the handlers of the statistics)
       */
      Statistics* statistics;
      //! \brief flag set by the threads which signal an event, if the
      //! statistics are the ones of a block (see `StatisticsBlock`)
      unsigned long long* signalled;
      //! \brief flag set by the first thread which signals an event (see
      //! `signalEvents`), or null
      int* events;
      //! \brief ordered key of the initial time step scaling factor
      unsigned long long rdt0;
    };

    /*!
     * \return a block of statistics (see `StatisticsBlock`), allocated in
     * the pinned memory of the host, the statistics being allocated in the
     * memory of the current device, or a null pointer if an allocation
     * failed. The block is reset before its first use. The allocations are
     * never freed.
     * \tparam Statistics: statistics updated by the devices
     */
    template <typename Statistics>
    StatisticsBlock<Statistics>* allocateStatisticsBlock() {
      auto* const device =
          allocate<DeviceStatistics<Statistics>>(getDeviceMemoryAllocator());
      if (device == nullptr) {
        return nullptr;
      }
      auto* const block =
          allocate<StatisticsBlock<Statistics>>(getPinnedMemoryAllocator());
      if (block == nullptr) {
        return nullptr;
      }
      return ::new (block) StatisticsBlock<Statistics>{
          .device = device, .neutral = false, .events = 0, .statistics = {}};
    }  // end of allocateStatisticsBlock

    /*!
     * \return the pool of blocks holding the statistics of the synchronous
     * computations on the accelerators (see `StatisticsBlock`), shared by all
     * the functions of the library. The blocks are identified by the device
     * of the queue and are kept until the end of the process.
     * \tparam Statistics: statistics updated by the devices
     */
    template <typename Statistics>
    auto& getStatisticsBlocks() noexcept {
      using Pool =
          StatusBlockPool<StatisticsBlock<Statistics>, ::Arcane::Int32>;
      static auto pool = Pool{};
      return pool;
    }  // end of getStatisticsBlocks

  }  // end of namespace internals

  /*!
   * \return the policy given by the calling code, checked against the
   * policy of the backend (see
   * `mfront::generic_parallel::checkBackendExecutionPolicy`). The `arcane`
   * backend has no default policy: a null policy is refused, as a queue
   * which is null or whose execution policy is not supported by the
   * library (see `checkQueue`).
   * \param[in] p: policy given by the calling code
   * \param[in] has_statistics: true if the calling code gave statistics
   */
  inline CheckedBackendExecutionPolicy<ArcaneBackendExecutionPolicy>
  checkExecutionPolicy(const mfront_gp_AbstractBackendExecutionPolicy* const p,
                       const bool has_statistics) {
    auto r = checkBackendExecutionPolicy<ArcaneBackendExecutionPolicy>(
        p, has_statistics);
    if (r.error != nullptr) {
      return r;
    }
    if (r.policy == nullptr) {
      r.error =
          "no execution policy given: the arcane backend requires a queue";
      return r;
    }
    r.error = checkQueue(r.policy->queue);
    return r;
  }  // end of checkExecutionPolicy

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * The statistics are updated by atomic operations, only performed if an
   * event occurs at a point. If the calling code gives statistics, they are
   * updated, and this function returns once the commands are launched if
   * the computations are asynchronous. Otherwise, the statistics are held
   * by a local variable with the queues of the host, and otherwise by a
   * block of a pool, the statistics being in the memory of the device of the
   * queue (see `StatisticsBlock`), as by the GPU backends: a call without
   * event performs neither copy nor additional synchronization. A
   * synchronous computation waits for the completion of the
   * commands of the queue. The device of a CUDA or HIP queue is made the
   * current device of the calling thread during the call (stack size of the
   * threads, allocation of the blocks).
   *
   * \tparam StatisticsHandler: handler of the statistics (see
   * `MFront/GenericParallel/Batch.hxx`)
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] e: execution of the batch. The policy must be an
   * `ArcaneBackendExecutionPolicy`.
   * \param[in] h: handler of the statistics
   * \param[in] f: function
   * \return the status
   */
  template <typename StatisticsHandler, typename PointFunction>
  int apply(typename StatisticsHandler::OutputStatus& s,
            const size_type n,
            const BatchExecution<typename StatisticsHandler::Statistics>& e,
            const StatisticsHandler& h,
            const PointFunction& f) {
    using Statistics = typename StatisticsHandler::Statistics;
    using Update = internals::AtomicUpdate<Statistics, PointFunction>;
    if (n == 0) {
      return s.status;
    }
    const auto p = checkExecutionPolicy(e.policy, e.statistics != nullptr);
    if (p.error != nullptr) {
      return report(s, -5, p.error);
    }
    auto queue = p.policy->queue;
    const auto current_device = CurrentDeviceGuard(queue);
    if (current_device.failed) {
      return report(s, -7, "unable to select the device of the queue");
    }
    const auto error =
        ensureStackSize(queue, f, p.common_options.preserve_stack_size);
    if (!error.empty()) {
      return report(s, -7, error);
    }
    const auto rdt0 = toOrderedKey(h.initial().rdt);
    if ((e.statistics != nullptr) || (!h.hasEvents())) {
      // statistics given by the calling code: the result is given by the
      // report function of the library. Without event, the statistics, if
      // any, are not accessed by the kernel, as for the GPU backends
      launch(queue, n, Update{f, e.statistics, nullptr, nullptr, rdt0});
      if (!p.common_options.asynchronous) {
        queue.barrier();
      }
      return s.status;
    }
    // synchronous computations: the statistics are held by a local variable
    // with the queues of the host, and otherwise by a block of the pool (see
    // `StatisticsBlock`)
    auto statistics = Statistics{};
    h.initialize(statistics);
    if (!queue.isAcceleratorPolicy()) {
      launch(queue, n, Update{f, &statistics, nullptr, nullptr, rdt0});
      queue.barrier();
      return h.report(s, statistics);
    }
    const auto block = internals::getStatisticsBlocks<Statistics>().acquire(
        queue.deviceId().asInt32(),
        internals::allocateStatisticsBlock<Statistics>);
    auto* const b = block.get();
    if (b == nullptr) {
      return report(s, -7, "unable to allocate the statistics of the batch");
    }
    // the copies of Arccore are synchronous
    using ::Arcane::Accelerator::MemoryCopyArgs;
    if (!b->neutral) {
      // reset, after an event or at the first use
      const auto neutral = getNeutralStatistics(statistics);
      queue.copyMemory(MemoryCopyArgs(
          b->device, &neutral, static_cast<::Arcane::Int64>(sizeof(neutral))));
    }
    b->neutral = false;
    b->events = 0;
    launch(queue, n,
           Update{f, &(b->device->statistics), &(b->device->signalled),
                  &(b->events), rdt0});
    queue.barrier();
    if (b->events == 0) {
      b->neutral = true;
      return h.report(s, statistics);
    }
    queue.copyMemory(
        MemoryCopyArgs(&(b->statistics), &(b->device->statistics),
                       static_cast<::Arcane::Int64>(sizeof(Statistics))));
    return h.report(s, mergeBlockStatistics(statistics, b->statistics));
  }  // end of apply

  /*!
   * \return the primitive of the backend, called by the functions of the
   * interfaces (see for instance `mfront::gpb::integrateBatch`)
   */
  inline auto makeApply() {
    return [](auto& s, const size_type n, const auto& e, const auto& h,
              const auto& f) { return arcane::apply(s, n, e, h, f); };
  }  // end of makeApply

}  // end of namespace mfront::generic_parallel::arcane

#endif /* LIB_MFRONT_GENERICPARALLEL_ARCANE_HXX */
