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

#include "MFront/GenericParallel/Batch.hxx"
#include "MFront/GenericParallel/StatusBlockPool.hxx"
#include "MFront/GenericParallel/ArcaneUtilities.hxx"
#include "MFront/GenericParallel/ArcaneBackendExecutionPolicy.hxx"

namespace mfront::generic_parallel::arcane {

  namespace internals {

    //! \brief atomic operations of the statistics
    struct Atomics {
      ARCCORE_HOST_DEVICE static void add(unsigned long long* const p,
                                          const unsigned long long v) {
        ::mfront::generic_parallel::arcane::atomicAdd(p, v);
      }
      ARCCORE_HOST_DEVICE static void min(unsigned long long* const p,
                                          const unsigned long long v) {
        ::mfront::generic_parallel::arcane::atomicMin(p, v);
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
        updateStatistics<Atomics>(*(this->statistics), this->f(i), this->rdt0);
      }
      //! \brief function object applied at one point
      PointFunction f;
      //! \brief statistics
      Statistics* statistics;
      //! \brief ordered key of the initial time step scaling factor
      unsigned long long rdt0;
    };

    /*!
     * \return the pool of blocks of unified memory holding the statistics
     * of the synchronous computations, shared by all the functions of the
     * library. The blocks are identified by the allocator of Arccore and
     * are kept until the end of the process.
     * \tparam Statistics: statistics updated by the devices
     */
    template <typename Statistics>
    auto& getStatisticsBlocks() noexcept {
      using Pool = StatusBlockPool<Statistics, ::Arcane::IMemoryAllocator*>;
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
   * by a block of unified memory, initialized and read by the host. In both
   * cases, a synchronous computation waits for the completion of the
   * commands of the queue (`barrier`).
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
    const auto error =
        ensureStackSize(queue, f, p.common_options.preserve_stack_size);
    if (!error.empty()) {
      return report(s, -7, error);
    }
    const auto rdt0 = toOrderedKey(h.initial().rdt);
    if (e.statistics != nullptr) {
      launch(queue, n, Update{f, e.statistics, rdt0});
      if (!p.common_options.asynchronous) {
        queue.barrier();
      }
      return s.status;
    }
    // synchronous computations: the device only accesses the block if an
    // event occurs, as for the GPU backends
    auto* const allocator = getUnifiedMemoryAllocator();
    const auto block = internals::getStatisticsBlocks<Statistics>().acquire(
        allocator, [allocator] { return allocate<Statistics>(allocator); });
    auto* const statistics = block.get();
    if (statistics == nullptr) {
      return report(s, -7, "unable to allocate the statistics of the batch");
    }
    h.initialize(*statistics);
    launch(queue, n, Update{f, statistics, rdt0});
    queue.barrier();
    return h.report(s, *statistics);
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
