/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/Arcane.hxx
 * \brief  This file declares the primitive of the `arcane` backend of the
 * `generic-parallel` behaviour interface, based on the accelerator API of
 * Arccore.
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
#include "MFront/GenericParallel/ArcaneUtilities.hxx"
#include "MFront/GenericParallel/ArcaneBackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/StatusBlockPool.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::arcane {

  namespace internals {

    /*!
     * \brief function object applying a function at one point and
     * accumulating the statistics of the batch by atomic operations, which
     * are only performed if the point failed, gave an unreliable result or
     * decreased the time step scaling factor.
     * \tparam PointFunction: function object applied at one point
     */
    template <typename PointFunction>
    struct AtomicUpdate {
      //! \brief apply the function at the `i`-th point
      ARCCORE_HOST_DEVICE void operator()(const size_type i) const {
        namespace utilities = ::mfront::generic_parallel::arcane;
        const auto r = this->f(i);
        if (r.number_of_failures != 0) {
          utilities::atomicAdd(&(this->statistics->number_of_failures), 1ull);
          utilities::atomicMin(&(this->statistics->first_failure),
                               r.first_failure);
        }
        if (r.number_of_unreliable_results != 0) {
          utilities::atomicAdd(
              &(this->statistics->number_of_unreliable_results), 1ull);
          utilities::atomicMin(&(this->statistics->first_unreliable_result),
                               r.first_unreliable_result);
        }
        const auto rdt = toOrderedKey(r.rdt);
        if (rdt < this->initial_rdt) {
          utilities::atomicMin(&(this->statistics->rdt), rdt);
        }
      }
      //! \brief function object applied at one point
      PointFunction f;
      //! \brief statistics
      AtomicBatchStatistics* statistics;
      //! \brief ordered key of the initial time step scaling factor
      unsigned long long initial_rdt;
    };

    /*!
     * \return the pool of blocks of unified memory holding the statistics
     * of the synchronous computations, shared by all the functions of the
     * library. The blocks are identified by the allocator of Arccore and
     * are kept until the end of the process.
     */
    inline auto& getStatisticsBlocks() noexcept {
      using Pool = mfront::generic_parallel::StatusBlockPool<
          AtomicBatchStatistics, ::Arcane::IMemoryAllocator*>;
      static auto pool = Pool{};
      return pool;
    }  // end of getStatisticsBlocks

    /*!
     * \return the policy given by the calling code, checked against the
     * policy of the backend (see
     * `mfront::generic_parallel::checkBackendExecutionPolicy`). The `arcane`
     * backend has no default policy: a null policy is refused, as a queue
     * which is null or whose execution policy is not supported by the
     * library (see `mfront::generic_parallel::arcane::checkQueue`).
     * \param[in] p: policy given by the calling code
     * \param[in] has_statistics: true if the calling code gave statistics
     */
    inline ::mfront::generic_parallel::CheckedBackendExecutionPolicy<
        ::mfront::generic_parallel::ArcaneBackendExecutionPolicy>
    checkExecutionPolicy(
        const mfront_gp_AbstractBackendExecutionPolicy* const p,
        const bool has_statistics) {
      auto r = ::mfront::generic_parallel::checkBackendExecutionPolicy<
          ::mfront::generic_parallel::ArcaneBackendExecutionPolicy>(
          p, has_statistics);
      if (r.error != nullptr) {
        return r;
      }
      if (r.policy == nullptr) {
        r.error =
            "no execution policy given: the arcane backend requires a queue";
        return r;
      }
      r.error = ::mfront::generic_parallel::arcane::checkQueue(r.policy->queue);
      return r;
    }  // end of checkExecutionPolicy

  }  // end of namespace internals

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * The statistics are accumulated by atomic operations, only performed if
   * an event occurs at a point. If the calling code gives statistics, they
   * are accumulated in them, and this function returns once the commands
   * are launched if the computations are asynchronous. Otherwise, the
   * statistics are held by a block of unified memory, initialized and read
   * by the host. In both cases, a synchronous computation waits for the
   * completion of the commands of the queue (`barrier`).
   *
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] e: execution of the batch. The policy must be an
   * `ArcaneBackendExecutionPolicy`.
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \return the status
   */
  template <typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            const BatchExecution& e,
            const real rdt0,
            const PointFunction& f) {
    namespace utilities = ::mfront::generic_parallel::arcane;
    if (n == 0) {
      return s.status;
    }
    const auto p =
        internals::checkExecutionPolicy(e.policy, e.statistics != nullptr);
    if (p.error != nullptr) {
      return report(s, -5, p.error);
    }
    auto queue = p.policy->queue;
    const auto error = utilities::ensureStackSize(
        queue, f, p.common_options.preserve_stack_size);
    if (!error.empty()) {
      return report(s, -7, error);
    }
    const auto rdt_key = toOrderedKey(rdt0);
    if (e.statistics != nullptr) {
      utilities::launch(
          queue, n,
          internals::AtomicUpdate<PointFunction>{f, e.statistics, rdt_key});
      if (!p.common_options.asynchronous) {
        queue.barrier();
      }
      return s.status;
    }
    // synchronous computations: the device only accesses the block if an
    // event occurs, as for the GPU backends
    auto* const allocator = utilities::getUnifiedMemoryAllocator();
    const auto block =
        internals::getStatisticsBlocks().acquire(allocator, [allocator] {
          return utilities::allocate<AtomicBatchStatistics>(allocator);
        });
    auto* const statistics = block.get();
    if (statistics == nullptr) {
      return report(s, -7, "unable to allocate the statistics of the batch");
    }
    mfront_gpb_initializeStatistics(statistics, rdt0);
    utilities::launch(
        queue, n,
        internals::AtomicUpdate<PointFunction>{f, statistics, rdt_key});
    queue.barrier();
    return report(s, getBatchStatistics(*statistics));
  }  // end of apply

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
    const auto p = internals::checkExecutionPolicy(policy, true);
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

  /*!
   * \return the primitive of the backend, as expected by
   * `mfront::gpb::integrateBatch` and the similar functions of
   * `MFront/GenericParallelBehaviour/Integrate.hxx`, called by the generated
   * functions
   * \param[in,out] s: output status
   */
  inline auto makeApply(OutputStatus& s) {
    return [&s](const size_type n, const BatchExecution& e, const real rdt0,
                const auto& f) { return arcane::apply(s, n, e, rdt0, f); };
  }  // end of makeApply

}  // end of namespace mfront::gpb::arcane

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_ARCANE_HXX */
