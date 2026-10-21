/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/RAJA.hxx
 * \brief  This file implements the integration of a batch of integration
 * points with RAJA.
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_RAJA_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_RAJA_HXX

#include <stdexcept>
#include <type_traits>
#include "RAJA/RAJA.hpp"
#include "MFront/GenericParallel/RAJAExecutionPolicy.hxx"
#include "MFront/GenericParallel/StatusBlockPool.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::raja {

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
      RAJA_HOST_DEVICE void operator()(const size_type i) const {
        const auto r = this->f(i);
        if (r.number_of_failures != 0) {
          RAJA::atomicAdd<RAJA::auto_atomic>(
              &(this->statistics->number_of_failures), 1ull);
          RAJA::atomicMin<RAJA::auto_atomic>(&(this->statistics->first_failure),
                                             r.first_failure);
        }
        if (r.number_of_unreliable_results != 0) {
          RAJA::atomicAdd<RAJA::auto_atomic>(
              &(this->statistics->number_of_unreliable_results), 1ull);
          RAJA::atomicMin<RAJA::auto_atomic>(
              &(this->statistics->first_unreliable_result),
              r.first_unreliable_result);
        }
        const auto rdt = toOrderedKey(r.rdt);
        if (rdt < this->initial_rdt) {
          RAJA::atomicMin<RAJA::auto_atomic>(&(this->statistics->rdt), rdt);
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
     * \return the pool of blocks of managed memory holding the statistics of
     * the synchronous computations with the resources of the given type, by
     * device, shared by all the functions of the library. The blocks are
     * kept until the end of the process.
     * \tparam Resource: resource of camp
     */
    template <typename Resource>
    auto& getStatisticsBlocks() noexcept {
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
   * The statistics are accumulated by atomic operations, only performed if
   * an event occurs at a point. If the execution options give statistics,
   * this function returns once the computations are submitted. Otherwise,
   * the statistics are held by a block of managed memory, initialized and
   * read by the host, and this function waits for the completion of the
   * computations.
   *
   * \tparam ExecutionPolicy: execution policy of RAJA
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] o: execution options. The queue is a pointer to a resource of
   * camp matching the execution policy (`camp::resources::Cuda`, `Hip` or
   * `Host`). If null, the default stream of the current device is used with
   * CUDA and HIP, the default resource otherwise (see
   * `mfront::generic_parallel::raja::getDefaultResource`).
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \return the status
   */
  template <typename ExecutionPolicy, typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            const ExecutionOptions& o,
            const real rdt0,
            const PointFunction& f) {
    using Resource = mfront::generic_parallel::raja::Resource<ExecutionPolicy>;
    if (n == 0) {
      return s.status;
    }
    const auto e =
        mfront::generic_parallel::raja::ensureStackSize<ExecutionPolicy>(
            f, o.preserve_stack_size != 0);
    if (!e.empty()) {
      return report(s, -7, e);
    }
    auto resource =
        mfront::generic_parallel::raja::getResource<ExecutionPolicy>(o.queue);
    const auto range = RAJA::TypedRangeSegment<size_type>(0, n);
    const auto rdt_key = toOrderedKey(rdt0);
    if (o.statistics != nullptr) {
      RAJA::forall<ExecutionPolicy>(
          resource, range,
          internals::AtomicUpdate<PointFunction>{f, o.statistics, rdt_key});
      return s.status;
    }
    // synchronous computations: the device only accesses the block if an
    // event occurs, as for the GPU backends
    const auto block = internals::getStatisticsBlocks<Resource>().acquire(
        mfront::generic_parallel::raja::getDevice(resource), [&resource] {
          return resource.template allocate<AtomicBatchStatistics>(
              1, camp::resources::MemoryAccess::Managed);
        });
    auto* const statistics = block.get();
    if (statistics == nullptr) {
      return report(s, -7, "unable to allocate the statistics of the batch");
    }
    mfront_gpb_initializeStatistics(statistics, rdt0);
    RAJA::forall<ExecutionPolicy>(
        resource, range,
        internals::AtomicUpdate<PointFunction>{f, statistics, rdt_key});
    resource.wait();
    return report(s, getBatchStatistics(*statistics));
  }  // end of apply

  /*!
   * \return a copy of a data in the memory of a resource (memory of the
   * device of the resource for the policies of CUDA and HIP), kept until the
   * end of the process
   * \param[in] resource: resource
   * \param[in] v: data
   */
  template <typename Resource, typename T>
  const T* copyToDevice(Resource& resource, const T& v) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "the shared data must be trivially copyable");
    auto* const p =
        resource.template allocate<T>(1, camp::resources::MemoryAccess::Device);
    if (p == nullptr) {
      throw std::runtime_error(
          "mfront::gpb::raja::copyToDevice: allocation failed");
    }
    resource.memcpy(p, &v, sizeof(T));
    resource.wait();
    return p;
  }  // end of copyToDevice

  /*!
   * \return a copy of a data shared by all the integration points (see the
   * `@SharedData` keyword) in the memory of the default resource of an
   * execution policy (memory of the current device for the policies of CUDA
   * and HIP). A copy is made at the first call on each device and kept until
   * the end of the process, as with the `cuda` and `hip` backends.
   * \tparam ExecutionPolicy: execution policy
   * \tparam get_shared_data: function returning the data
   */
  template <typename ExecutionPolicy, auto get_shared_data>
  const auto* getSharedDataCopy() {
    using T = std::remove_cvref_t<decltype(get_shared_data())>;
    static auto copies = SharedDataCopies<T, int>{};
    auto resource =
        mfront::generic_parallel::raja::getDefaultResource<ExecutionPolicy>();
    return copies.get(
        mfront::generic_parallel::raja::getDevice(resource),
        [&resource] { return copyToDevice(resource, get_shared_data()); });
  }  // end of getSharedDataCopy

  /*!
   * \return the primitive of the backend, as expected by
   * `mfront::gpb::integrateBatch` and the similar functions of
   * `MFront/GenericParallelBehaviour/Integrate.hxx`, called by the generated
   * functions
   * \tparam ExecutionPolicy: execution policy
   * \param[in,out] s: output status
   */
  template <typename ExecutionPolicy>
  auto makeApply(OutputStatus& s) {
    return [&s](const size_type n, const ExecutionOptions& o, const real rdt0,
                const auto& f) {
      return raja::apply<ExecutionPolicy>(s, n, o, rdt0, f);
    };
  }  // end of makeApply

}  // end of namespace mfront::gpb::raja

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_RAJA_HXX */
