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
#include "MFront/GenericParallelBehaviour/Integrate.hxx"
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#include "MFront/GenericParallelBehaviour/GPU.hxx"
#endif

namespace mfront::gpb::kokkos {

  namespace internals {

    /*!
     * \brief functor computing the statistics of a batch with
     * `Kokkos::parallel_reduce`.
     * \tparam PointFunction: function object applied at one point
     */
    template <typename PointFunction>
    struct Reduction {
      //! \brief type of the reduced value
      using value_type = BatchStatistics;
      //! \brief apply the function at the `i`-th point
      KOKKOS_FUNCTION void operator()(const size_type i, value_type& s) const {
        s = combine(s, this->f(i));
      }
      //! \brief combine two partial results
      KOKKOS_FUNCTION void join(value_type& s, const value_type& s2) const {
        s = combine(s, s2);
      }
      //! \brief initialize a partial result
      KOKKOS_FUNCTION void init(value_type& s) const {
        s = value_type{.rdt = this->initial_rdt};
      }
      //! \brief function object applied at one point
      PointFunction f;
      //! \brief initial value of the time step scaling factor
      real initial_rdt;
    };

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

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)

    /*!
     * \brief kernel which is never launched: its attributes give the stack
     * size required to apply a function at one point. Kokkos does not
     * increase the stack size of the threads, which is too small by default
     * for some behaviours.
     */
    template <typename PointFunction>
    __global__ void stackSizeProbe(const PointFunction f) {
      f(0);
    }  // end of stackSizeProbe

#endif /* defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) */

    /*!
     * \brief prepare the execution of the given function in the given
     * execution space.
     * \return true on success
     * \param[out] s: output status
     * \param[in] preserve_stack_size: if true, the stack size of the
     * threads of the device is not modified
     */
    template <typename ExecutionSpace, typename PointFunction>
    bool prepare(OutputStatus& s, const bool preserve_stack_size) {
#if defined(KOKKOS_ENABLE_CUDA)
      constexpr auto device = std::is_same_v<ExecutionSpace, Kokkos::Cuda>;
#elif defined(KOKKOS_ENABLE_HIP)
      constexpr auto device = std::is_same_v<ExecutionSpace, Kokkos::HIP>;
#else
      constexpr auto device = false;
#endif
      if constexpr (device) {
        return gpu::internals::ensureStackSize(s, stackSizeProbe<PointFunction>,
                                               preserve_stack_size);
      }
      static_cast<void>(s);
      static_cast<void>(preserve_stack_size);
      return true;
    }  // end of prepare

  }  // end of namespace internals

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * Kokkos must have been initialized by the calling code. If the execution
   * options give statistics, the statistics are accumulated by atomic
   * operations and this function returns once the computations are
   * submitted. Otherwise, the statistics are computed by a parallel
   * reduction and this function waits for the completion of the
   * computations.
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
    if (!internals::prepare<ExecutionSpace, PointFunction>(
            s, o.preserve_stack_size != 0)) {
      return s.status;
    }
    const auto space = o.queue == nullptr
                           ? ExecutionSpace{}
                           : *static_cast<const ExecutionSpace*>(o.queue);
    const auto range =
        Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<size_type>>(space,
                                                                          0, n);
    if (o.statistics != nullptr) {
      Kokkos::parallel_for("mfront::gpb::kokkos::apply", range,
                           internals::AtomicUpdate<PointFunction>{
                               f, o.statistics, toOrderedKey(rdt0)});
      return s.status;
    }
    // the reduction in a scalar is synchronous
    auto r = BatchStatistics{};
    Kokkos::parallel_reduce("mfront::gpb::kokkos::apply", range,
                            internals::Reduction<PointFunction>{f, rdt0}, r);
    return report(s, r);
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
