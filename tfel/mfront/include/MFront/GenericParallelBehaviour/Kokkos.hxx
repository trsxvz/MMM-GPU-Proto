/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/Kokkos.hxx
 * \brief  This file implements the integration of a batch of integration
 * points with Kokkos.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_KOKKOS_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_KOKKOS_HXX

#include <stdexcept>
#include <type_traits>
#include <Kokkos_Core.hpp>
#include "MFront/GenericParallel/KokkosStackSize.hxx"
#include "MFront/GenericParallel/StatusBlockPool.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

static_assert(Kokkos::has_shared_space,
              "the Kokkos backend of the generic-parallel interface requires "
              "a shared memory space (Kokkos::SharedSpace)");

namespace mfront::gpb::kokkos {

  namespace internals {

    /*!
     * \brief functor applying a function at one point and accumulating the
     * statistics of the batch by atomic operations, which are only performed
     * if the point failed, gave an unreliable result or decreased the time
     * step scaling factor.
     * \tparam PointFunction: function object applied at one point
     */
    template <typename PointFunction>
    struct AtomicUpdate {
      //! \brief apply the function at the `i`-th point
      KOKKOS_FUNCTION void operator()(const size_type i) const {
        const auto r = this->f(i);
        if (r.number_of_failures != 0) {
          Kokkos::atomic_add(&(this->statistics->number_of_failures), 1ull);
          Kokkos::atomic_min(&(this->statistics->first_failure_index),
                             static_cast<unsigned long long>(i));
        }
        if (r.number_of_unreliable_results != 0) {
          Kokkos::atomic_add(&(this->statistics->number_of_unreliable_results),
                             1ull);
        }
        const auto rdt = toOrderedKey(r.rdt);
        if (rdt < this->initial_rdt) {
          Kokkos::atomic_min(&(this->statistics->rdt), rdt);
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
     * \return the pool of blocks of shared memory holding the statistics of
     * the synchronous computations, shared by all the functions of the
     * library. Kokkos uses one device per process. The blocks are freed when
     * Kokkos is finalized.
     */
    inline auto& getStatisticsBlocks() noexcept {
      using Pool =
          mfront::generic_parallel::StatusBlockPool<AtomicBatchStatistics, int>;
      static auto pool = Pool{};
      return pool;
    }  // end of getStatisticsBlocks

  }  // end of namespace internals

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * Kokkos must have been initialized by the calling code. The statistics
   * are accumulated by atomic operations, only performed if an event occurs
   * at a point. If the execution options give statistics, this function
   * returns once the computations are submitted. Otherwise, the statistics
   * are held by a block of shared memory, initialized and read by the host,
   * and this function waits for the completion of the computations.
   *
   * \tparam ExecutionSpace: execution space
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] o: execution options. The queue is a pointer to an instance
   * of the execution space (default instance if null).
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \return the status
   */
  template <typename ExecutionSpace, typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            const ExecutionOptions& o,
            const real rdt0,
            const PointFunction& f) {
    if (!Kokkos::is_initialized()) {
      return report(s, -3, "Kokkos is not initialized");
    }
    if (n == 0) {
      return s.status;
    }
    const auto e = mfront::generic_parallel::kokkos::ensureStackSize<
        ExecutionSpace>(f, o.preserve_stack_size != 0);
    if (!e.empty()) {
      return report(s, -3, e);
    }
    const auto space = o.queue == nullptr
                           ? ExecutionSpace{}
                           : *static_cast<const ExecutionSpace*>(o.queue);
    const auto range =
        Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<size_type>>(space,
                                                                          0, n);
    const auto rdt_key = toOrderedKey(rdt0);
    if (o.statistics != nullptr) {
      Kokkos::parallel_for(
          "mfront::gpb::kokkos::apply", range,
          internals::AtomicUpdate<PointFunction>{f, o.statistics, rdt_key});
      return s.status;
    }
    // synchronous computations: the device only accesses the block if an
    // event occurs, as for the GPU backends
    const auto block = internals::getStatisticsBlocks().acquire(0, [] {
      auto* const p = static_cast<AtomicBatchStatistics*>(
          Kokkos::kokkos_malloc<Kokkos::SharedSpace>(
              "mfront::gpb::kokkos::statistics",
              sizeof(AtomicBatchStatistics)));
      Kokkos::push_finalize_hook(
          [p] { Kokkos::kokkos_free<Kokkos::SharedSpace>(p); });
      return p;
    });
    auto* const statistics = block.get();
    if (statistics == nullptr) {
      return report(s, -3, "unable to allocate the statistics of the batch");
    }
    mfront_gpb_initializeBatchStatistics(statistics, rdt0);
    Kokkos::parallel_for(
        "mfront::gpb::kokkos::apply", range,
        internals::AtomicUpdate<PointFunction>{f, statistics, rdt_key});
    space.fence();
    return report(s, getBatchStatistics(*statistics));
  }  // end of apply

  /*!
   * \return a copy of a data in the memory space of an execution space,
   * freed when Kokkos is finalized
   * \tparam ExecutionSpace: execution space
   * \param[in] v: data
   */
  template <typename ExecutionSpace, typename T>
  const T* copyToDevice(const T& v) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "the shared data must be trivially copyable");
    using MemorySpace = typename ExecutionSpace::memory_space;
    using Bytes = Kokkos::View<char*, MemorySpace, Kokkos::MemoryUnmanaged>;
    using HostBytes =
        Kokkos::View<const char*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
    if (!Kokkos::is_initialized()) {
      throw std::runtime_error("Kokkos is not initialized");
    }
    auto* const p = static_cast<T*>(Kokkos::kokkos_malloc<MemorySpace>(
        "mfront::gpb::kokkos::copyToDevice", sizeof(T)));
    Kokkos::deep_copy(Bytes(reinterpret_cast<char*>(p), sizeof(T)),
                      HostBytes(reinterpret_cast<const char*>(&v), sizeof(T)));
    Kokkos::push_finalize_hook([p] { Kokkos::kokkos_free<MemorySpace>(p); });
    return p;
  }  // end of copyToDevice

  /*!
   * \return a copy of a data shared by all the integration points (see the
   * `@SharedData` keyword) in the memory space of an execution space, made
   * at the first call (Kokkos uses one device per process) and freed when
   * Kokkos is finalized.
   * \tparam ExecutionSpace: execution space
   * \tparam get_shared_data: function returning the data
   */
  template <typename ExecutionSpace, auto get_shared_data>
  const auto* getSharedDataCopy() {
    static const auto* const p =
        copyToDevice<ExecutionSpace>(get_shared_data());
    return p;
  }  // end of getSharedDataCopy

  /*!
   * \return the primitive of the backend, as expected by
   * `mfront::gpb::integrateBatch` and the similar functions of
   * `MFront/GenericParallelBehaviour/Integrate.hxx`, called by the generated
   * functions
   * \tparam ExecutionSpace: execution space
   * \param[in,out] s: output status
   */
  template <typename ExecutionSpace>
  auto makeApply(OutputStatus& s) {
    return [&s](const size_type n, const ExecutionOptions& o, const real rdt0,
                const auto& f) {
      return kokkos::apply<ExecutionSpace>(s, n, o, rdt0, f);
    };
  }  // end of makeApply

}  // end of namespace mfront::gpb::kokkos

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_KOKKOS_HXX */
