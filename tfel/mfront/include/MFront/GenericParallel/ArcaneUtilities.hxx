/*!
 * \file   mfront/include/MFront/GenericParallel/ArcaneUtilities.hxx
 * \brief  This file declares the functions shared by the `arcane` backends
 * of the `generic-parallel` interfaces: check of the queue, launch of the
 * loops, atomic operations, unified memory and stack size of the threads of
 * the device.
 * \date   23/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_ARCANEUTILITIES_HXX
#define LIB_MFRONT_GENERICPARALLEL_ARCANEUTILITIES_HXX

#include <atomic>
#include <string>
#include <algorithm>
#include <type_traits>
#include <arccore/base/ArccoreGlobal.h>
#include <arccore/common/MemoryUtils.h>
#include <arccore/common/IMemoryAllocator.h>
#include <arccore/common/AllocatedMemoryInfo.h>
#include <arccore/common/MemoryAllocationArgs.h>
#include <arccore/common/accelerator/RunQueue.h>
#include <arccore/accelerator/RunCommandLoop.h>
#if defined(ARCCORE_COMPILING_CUDA_OR_HIP)
#include "MFront/GenericParallel/GPUStackSize.hxx"
#endif

namespace mfront::generic_parallel::arcane {

  /*!
   * \return a message explaining why the commands can't be launched on a
   * queue by the code of the generated library, or a null pointer. The
   * `Thread` and `Sequential` queues are always supported, the CUDA (HIP)
   * queues only if the library is compiled by a CUDA (HIP) compiler.
   * \param[in] q: queue
   */
  inline const char* checkQueue(const ::Arcane::Accelerator::RunQueue& q) {
    using ::Arcane::Accelerator::eExecutionPolicy;
    if (q.isNull()) {
      return "the queue of the execution policy is null";
    }
    switch (q.executionPolicy()) {
      case eExecutionPolicy::Sequential:
      case eExecutionPolicy::Thread:
        return nullptr;
#if defined(ARCCORE_COMPILING_CUDA)
      case eExecutionPolicy::CUDA:
        return nullptr;
#endif
#if defined(ARCCORE_COMPILING_HIP)
      case eExecutionPolicy::HIP:
        return nullptr;
#endif
      default:
        break;
    }
    return "the execution policy of the queue is not supported by the "
           "library (the CUDA and HIP queues require a library compiled by "
           "a CUDA or HIP compiler)";
  }  // end of checkQueue

  /*!
   * \brief launch the commands applying a function at the points `[0, n[`
   * on a queue.
   *
   * The loops of Arccore have 32 bits indices, and the kernels of Arccore
   * compute the index of a thread on 32 bits: the ranges are split in
   * commands of at most `2^30` points, so that the index of a thread of the
   * grid of a command never overflows. The commands are launched
   * asynchronously, whatever the asynchronous state of the queue, which is
   * restored afterwards: the calling code waits for their completion
   * (`barrier`) if needed. With the `Thread` and `Sequential` queues, the
   * commands are completed on return. If the launch of a command fails, the
   * commands already launched are completed before the exception is
   * propagated, so that the memory they access may be released.
   *
   * \tparam Index: type of the indices taken by the function
   * \param[in] q: queue
   * \param[in] n: number of points
   * \param[in] f: function, copied by the commands
   */
  template <typename Index, typename Function>
  void launch(::Arcane::Accelerator::RunQueue& q,
              const Index n,
              const Function& f) {
    using ::Arcane::Int32;
    using ::Arcane::Int64;
    constexpr auto chunk = Int64{1} << 30;
    const auto size = static_cast<Int64>(n);
    ::Arcane::Accelerator::RunQueue::ScopedAsync asynchronous_launch(&q);
    try {
      for (auto offset = Int64{}; offset < size; offset += chunk) {
        const auto m = static_cast<Int32>(std::min(chunk, size - offset));
        auto command = ::Arcane::Accelerator::makeCommand(q);
        command << RUNCOMMAND_LOOP1(iter, m) {
          auto [i] = iter();
          f(static_cast<Index>(offset + i));
        };
      }
    } catch (...) {
      try {
        q.barrier();
      } catch (...) {
      }
      throw;
    }
  }  // end of launch

  /*!
   * \brief atomically add a value to a counter of the statistics. The
   * atomic operations of Arccore do not support the unsigned 64 bits
   * integers used by the statistics.
   * \param[in] p: counter
   * \param[in] v: value
   */
  ARCCORE_HOST_DEVICE inline void atomicAdd(unsigned long long* const p,
                                            const unsigned long long v) {
#if defined(ARCCORE_DEVICE_TARGET_CUDA) || defined(ARCCORE_DEVICE_TARGET_HIP)
    ::atomicAdd(p, v);
#else
    std::atomic_ref<unsigned long long>(*p).fetch_add(
        v, std::memory_order_relaxed);
#endif
  }  // end of atomicAdd

  /*!
   * \brief atomically replace a counter of the statistics by the minimum of
   * its value and of a given value (see `atomicAdd`)
   * \param[in] p: counter
   * \param[in] v: value
   */
  ARCCORE_HOST_DEVICE inline void atomicMin(unsigned long long* const p,
                                            const unsigned long long v) {
#if defined(ARCCORE_DEVICE_TARGET_CUDA) || defined(ARCCORE_DEVICE_TARGET_HIP)
    ::atomicMin(p, v);
#else
    auto c = std::atomic_ref<unsigned long long>(*p);
    auto current = c.load(std::memory_order_relaxed);
    while ((v < current) &&
           (!c.compare_exchange_weak(current, v, std::memory_order_relaxed))) {
    }
#endif
  }  // end of atomicMin

  /*!
   * \return the allocator of the unified memory of Arccore, accessible by
   * the host and by the devices, or a null pointer if it is not available.
   * If no accelerator runtime is initialized, Arccore returns the allocator
   * of the host.
   */
  inline ::Arcane::IMemoryAllocator* getUnifiedMemoryAllocator() noexcept {
    try {
      return ::Arcane::MemoryUtils::getAllocator(
          ::Arcane::eMemoryResource::UnifiedMemory);
    } catch (...) {
    }
    return nullptr;
  }  // end of getUnifiedMemoryAllocator

  /*!
   * \return an object allocated by an allocator of Arccore, which is never
   * freed, or a null pointer if the allocator is null
   * \tparam T: type of the object, which must be trivially copyable
   * \param[in] a: allocator
   */
  template <typename T>
  T* allocate(::Arcane::IMemoryAllocator* const a) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "the allocated objects must be trivially copyable");
    if (a == nullptr) {
      return nullptr;
    }
    const auto m = a->allocate(::Arcane::MemoryAllocationArgs{},
                               static_cast<::Arcane::Int64>(sizeof(T)));
    return static_cast<T*>(m.baseAddress());
  }  // end of allocate

#if defined(ARCCORE_COMPILING_CUDA_OR_HIP)

  /*!
   * \brief kernel which is never launched: its attributes give the stack
   * size required to apply a function at one point. Arccore does not
   * increase the stack size of the threads, which is too small by default
   * for some behaviours.
   */
  template <typename Function>
  __global__ void stackSizeProbe(const Function f) {
    f(0);
  }  // end of stackSizeProbe

#endif /* defined(ARCCORE_COMPILING_CUDA_OR_HIP) */

  /*!
   * \brief increase the stack size of the threads of the current device of
   * the calling thread (which must be the device of the queue) if the given
   * function, applied at one point on a CUDA or HIP queue, requires it.
   * Nothing is done for the other queues.
   * \return an error message, empty on success
   * \param[in] q: queue
   * \param[in] f: function applied at one point
   * \param[in] preserve: if true, the stack size is not modified, and an
   * error is returned if it is too small
   */
  template <typename Function>
  std::string ensureStackSize(const ::Arcane::Accelerator::RunQueue& q,
                              const Function& f,
                              const bool preserve) {
    static_cast<void>(q);
    static_cast<void>(f);
    static_cast<void>(preserve);
#if defined(ARCCORE_COMPILING_CUDA_OR_HIP)
    if (q.isAcceleratorPolicy()) {
      const auto r = gpu::ensureStackSize<stackSizeProbe<Function>>(preserve);
      if (!r.succeeded()) {
        return r.getErrorMessage();
      }
    }
#endif
    return {};
  }  // end of ensureStackSize

}  // end of namespace mfront::generic_parallel::arcane

#endif /* LIB_MFRONT_GENERICPARALLEL_ARCANEUTILITIES_HXX */
