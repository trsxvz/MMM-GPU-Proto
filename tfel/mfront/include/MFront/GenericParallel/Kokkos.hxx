/*!
 * \file   mfront/include/MFront/GenericParallel/Kokkos.hxx
 * \brief  This file declares the primitive of the `kokkos` backend of the
 * `generic-parallel` interfaces.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_KOKKOS_HXX
#define LIB_MFRONT_GENERICPARALLEL_KOKKOS_HXX

#include <Kokkos_Core.hpp>
#include "MFront/GenericParallel/Batch.hxx"
#include "MFront/GenericParallel/KokkosStackSize.hxx"
#include "MFront/GenericParallel/StatusBlockPool.hxx"
#include "MFront/GenericParallel/KokkosBackendExecutionPolicy.hxx"

static_assert(Kokkos::has_shared_space,
              "the Kokkos backend of the generic-parallel interface requires "
              "a shared memory space (Kokkos::SharedSpace)");

namespace mfront::generic_parallel::kokkos {

  namespace internals {

    //! \brief atomic operations of the statistics
    struct Atomics {
      KOKKOS_FUNCTION static void add(unsigned long long* const p,
                                      const unsigned long long v) {
        Kokkos::atomic_add(p, v);
      }
      KOKKOS_FUNCTION static void min(unsigned long long* const p,
                                      const unsigned long long v) {
        Kokkos::atomic_min(p, v);
      }
    };

    /*!
     * \brief functor applying a function at one point and updating the
     * statistics of the batch by atomic operations, which are only performed
     * if an event occurs at this point (see `updateStatistics`).
     * \tparam Statistics: statistics updated by the devices
     * \tparam PointFunction: function object applied at one point
     */
    template <typename Statistics, typename PointFunction>
    struct AtomicUpdate {
      //! \brief apply the function at the `i`-th point
      KOKKOS_FUNCTION void operator()(const size_type i) const {
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
     * \return the pool of blocks of shared memory holding the statistics of
     * the synchronous computations, shared by all the functions of the
     * library. Kokkos uses one device per process. The blocks are freed when
     * Kokkos is finalized.
     * \tparam Statistics: statistics updated by the devices
     */
    template <typename Statistics>
    auto& getStatisticsBlocks() noexcept {
      using Pool = StatusBlockPool<Statistics, int>;
      static auto pool = Pool{};
      return pool;
    }  // end of getStatisticsBlocks

  }  // end of namespace internals

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * Kokkos must have been initialized by the calling code. The statistics
   * are updated by atomic operations, only performed if an event occurs at
   * a point. If the calling code gives statistics, they are updated, and
   * this function returns once the computations are submitted if the
   * computations are asynchronous. Otherwise, the statistics are held by a
   * block of shared memory, initialized and read by the host. In both cases,
   * a synchronous computation waits for the completion of the computations.
   *
   * \tparam ExecutionSpace: execution space
   * \tparam StatisticsHandler: handler of the statistics (see
   * `MFront/GenericParallel/Batch.hxx`)
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] e: execution of the batch. The policy is a
   * `KokkosBackendExecutionPolicy<ExecutionSpace>`, or null for the default
   * instance of the execution space.
   * \param[in] h: handler of the statistics
   * \param[in] f: function
   * \return the status
   */
  template <typename ExecutionSpace,
            typename StatisticsHandler,
            typename PointFunction>
  int apply(typename StatisticsHandler::OutputStatus& s,
            const size_type n,
            const BatchExecution<typename StatisticsHandler::Statistics>& e,
            const StatisticsHandler& h,
            const PointFunction& f) {
    using Statistics = typename StatisticsHandler::Statistics;
    using Update = internals::AtomicUpdate<Statistics, PointFunction>;
    using ExecutionPolicy = KokkosBackendExecutionPolicy<ExecutionSpace>;
    if (n == 0) {
      return s.status;
    }
    if (!Kokkos::is_initialized()) {
      return report(s, -7, "Kokkos is not initialized");
    }
    const auto p = checkBackendExecutionPolicy<ExecutionPolicy>(
        e.policy, e.statistics != nullptr);
    if (p.error != nullptr) {
      return report(s, -5, p.error);
    }
    const auto error = ensureStackSize<ExecutionSpace>(
        f, p.common_options.preserve_stack_size);
    if (!error.empty()) {
      return report(s, -7, error);
    }
    const auto space =
        p.policy == nullptr ? ExecutionSpace{} : p.policy->execution_space;
    const auto range =
        Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<size_type>>(space,
                                                                          0, n);
    const auto rdt0 = toOrderedKey(h.initial().rdt);
    if (e.statistics != nullptr) {
      Kokkos::parallel_for("mfront::generic_parallel::kokkos::apply", range,
                           Update{f, e.statistics, rdt0});
      if (!p.common_options.asynchronous) {
        space.fence();
      }
      return s.status;
    }
    // synchronous computations: the device only accesses the block if an
    // event occurs, as for the GPU backends
    const auto block =
        internals::getStatisticsBlocks<Statistics>().acquire(0, [] {
          auto* const p = static_cast<Statistics*>(
              Kokkos::kokkos_malloc<Kokkos::SharedSpace>(
                  "mfront::generic_parallel::kokkos::statistics",
                  sizeof(Statistics)));
          Kokkos::push_finalize_hook(
              [p] { Kokkos::kokkos_free<Kokkos::SharedSpace>(p); });
          return p;
        });
    auto* const statistics = block.get();
    if (statistics == nullptr) {
      return report(s, -7, "unable to allocate the statistics of the batch");
    }
    h.initialize(*statistics);
    Kokkos::parallel_for("mfront::generic_parallel::kokkos::apply", range,
                         Update{f, statistics, rdt0});
    space.fence();
    return h.report(s, *statistics);
  }  // end of apply

  /*!
   * \return the primitive of the backend, called by the functions of the
   * interfaces (see for instance `mfront::gpb::integrateBatch`)
   * \tparam ExecutionSpace: execution space
   */
  template <typename ExecutionSpace>
  auto makeApply() {
    return [](auto& s, const size_type n, const auto& e, const auto& h,
              const auto& f) {
      return kokkos::apply<ExecutionSpace>(s, n, e, h, f);
    };
  }  // end of makeApply

}  // end of namespace mfront::generic_parallel::kokkos

#endif /* LIB_MFRONT_GENERICPARALLEL_KOKKOS_HXX */
