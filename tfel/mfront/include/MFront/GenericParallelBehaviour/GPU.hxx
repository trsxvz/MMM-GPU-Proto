/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/GPU.hxx
 * \brief  This file implements the integration of a batch of integration
 * points with CUDA or HIP. This file must be compiled by a CUDA or a HIP
 * compiler. The HIP runtime is used if the `hip/hip_runtime.h` header has
 * been included before this file or if the HIP language is enabled, the
 * CUDA runtime otherwise. Including the HIP runtime header first is required
 * to use HIP on NVIDIA GPUs, since `hipcc` then calls `nvcc`.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_GPU_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_GPU_HXX

#include <mutex>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include "MFront/GenericParallel/GPUStackSize.hxx"
#define MFRONT_GPB_GPU_API(f) MFRONT_GENERIC_PARALLEL_GPU_API(f)
#include "MFront/GenericParallel/StatusBlockPool.hxx"
#if defined(__HIPCC__) || defined(__HIP__)
#include "MFront/GenericParallel/HIPBackendExecutionPolicy.hxx"
#else
#include "MFront/GenericParallel/CUDABackendExecutionPolicy.hxx"
#endif
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::gpu {

  namespace internals {

#if defined(__HIPCC__) || defined(__HIP__)
    //! \brief device attribute giving the warp size
    inline constexpr auto warp_size_attribute = hipDeviceAttributeWarpSize;
    //! \brief execution policy of the backend
    using ExecutionPolicy =
        ::mfront::generic_parallel::HIPBackendExecutionPolicy;
#else
    //! \brief device attribute giving the warp size
    inline constexpr auto warp_size_attribute = cudaDevAttrWarpSize;
    //! \brief execution policy of the backend
    using ExecutionPolicy =
        ::mfront::generic_parallel::CUDABackendExecutionPolicy;
#endif

    /*!
     * \brief report an error of the runtime. Only the first error is
     * reported.
     * \return true if the given code denotes an error
     * \param[in,out] s: output status
     * \param[in] e: error code
     */
    inline bool failed(OutputStatus& s,
                       const MFRONT_GPB_GPU_API(Error_t) e) noexcept {
      if (e == MFRONT_GPB_GPU_API(Success)) {
        return false;
      }
      if (s.status != -7) {
        report(s, -7, MFRONT_GPB_GPU_API(GetErrorString)(e));
      }
      return true;
    }  // end of failed

    /*!
     * \brief increase the stack size of the threads if the given kernel
     * requires it (see `mfront::generic_parallel::gpu::ensureStackSize`).
     * \return true on success
     * \tparam kernel: kernel
     * \param[in,out] s: output status
     * \param[in] preserve: if true, the stack size is not modified and an
     * error is reported if it is too small
     */
    template <auto kernel>
    bool ensureStackSize(OutputStatus& s, const bool preserve) {
      const auto r =
          mfront::generic_parallel::gpu::ensureStackSize<kernel>(preserve);
      if (failed(s, r.error)) {
        return false;
      }
      if (r.required_stack_size != 0) {
        report(s, -7, r.getErrorMessage());
        return false;
      }
      return true;
    }  // end of ensureStackSize

    //! \brief sum of the values of a warp, in the first thread of the warp
    template <typename T>
    __device__ inline T warpSum(T v) {
      for (int o = warpSize / 2; o > 0; o /= 2) {
#if defined(__HIP_PLATFORM_AMD__)
        v += __shfl_down(v, o);
#else
        v += __shfl_down_sync(0xffffffffu, v, o);
#endif
      }
      return v;
    }  // end of warpSum

    //! \return true if the predicate is true for a thread of the warp
    __device__ inline bool warpAny(const bool p) {
#if defined(__HIP_PLATFORM_AMD__)
      return __any(p);
#else
      return __any_sync(0xffffffffu, p);
#endif
    }  // end of warpAny

    //! \brief minimum of the values of a warp, in the first thread of the warp
    __device__ inline unsigned long long warpMin(unsigned long long v) {
      for (int o = warpSize / 2; o > 0; o /= 2) {
#if defined(__HIP_PLATFORM_AMD__)
        const auto w = __shfl_down(v, o);
#else
        const auto w = __shfl_down_sync(0xffffffffu, v, o);
#endif
        v = w < v ? w : v;
      }
      return v;
    }  // end of warpMin

    /*!
     * \brief apply a function at the points of a batch, one thread per
     * point, and update the statistics of the batch. The number of threads
     * per block must be a multiple of the warp size.
     *
     * \param[in] f: function
     * \param[in] n: number of points
     * \param[out] s: statistics
     * \param[in] rdt0: ordered key of the initial time step scaling factor
     */
    template <typename PointFunction>
    __global__ void apply(const PointFunction f,
                          const size_type n,
                          AtomicBatchStatistics* const s,
                          const unsigned long long rdt0) {
      const auto i = static_cast<size_type>(blockIdx.x) * blockDim.x +
                     static_cast<size_type>(threadIdx.x);
      auto failures = 0u;
      auto first_failure = ::mfront::generic_parallel::no_event;
      auto unreliable_results = 0u;
      auto first_unreliable_result = ::mfront::generic_parallel::no_event;
      auto rdt = rdt0;
      if (i < n) {
        const auto r = f(i);
        failures = static_cast<unsigned int>(r.number_of_failures);
        first_failure = r.first_failure;
        unreliable_results =
            static_cast<unsigned int>(r.number_of_unreliable_results);
        first_unreliable_result = r.first_unreliable_result;
        rdt = toOrderedKey(r.rdt);
      }
      // the statistics are only reduced and updated if a point of the warp
      // failed, gave an unreliable result or decreased the time step
      // scaling factor: otherwise, the statistics are not accessed
      if (!warpAny((failures != 0) || (unreliable_results != 0) ||
                   (rdt < rdt0))) {
        return;
      }
      failures = warpSum(failures);
      first_failure = warpMin(first_failure);
      unreliable_results = warpSum(unreliable_results);
      first_unreliable_result = warpMin(first_unreliable_result);
      rdt = warpMin(rdt);
      if (threadIdx.x % warpSize == 0) {
        if (failures != 0) {
          atomicAdd(&(s->number_of_failures),
                    static_cast<unsigned long long>(failures));
          atomicMin(&(s->first_failure), first_failure);
        }
        if (unreliable_results != 0) {
          atomicAdd(&(s->number_of_unreliable_results),
                    static_cast<unsigned long long>(unreliable_results));
          atomicMin(&(s->first_unreliable_result), first_unreliable_result);
        }
        if (rdt < rdt0) {
          atomicMin(&(s->rdt), rdt);
        }
      }
    }  // end of apply

    /*!
     * \return the pool of blocks of managed memory holding the statistics
     * of the synchronous computations, shared by all the functions of the
     * library
     */
    inline auto& getStatisticsBlocks() noexcept {
      using Pool =
          mfront::generic_parallel::StatusBlockPool<AtomicBatchStatistics, int>;
      static auto pool = Pool{};
      return pool;
    }  // end of getStatisticsBlocks

  }  // end of namespace internals

  /*!
   * \brief apply a function at the points of a batch, one thread per point,
   * and report the statistics of the batch.
   *
   * The kernel is launched on the stream given by the execution policy. If
   * the computations are asynchronous, this function returns once the
   * kernel is launched. Otherwise, it waits for its completion.
   *
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] e: execution of the batch. The policy is a
   * `CUDABackendExecutionPolicy` (`HIPBackendExecutionPolicy` with HIP), or
   * null for the default stream.
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \param[in] threads_per_block: number of threads per block, rounded up to
   * a multiple of the warp size
   * \return the status
   */
  template <typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            const BatchExecution& e,
            const real rdt0,
            const PointFunction& f,
            const unsigned int threads_per_block) {
    using Stream = MFRONT_GPB_GPU_API(Stream_t);
    using internals::failed;
    if (n == 0) {
      return s.status;
    }
    const auto p = ::mfront::generic_parallel::checkBackendExecutionPolicy<
        internals::ExecutionPolicy>(e.policy, e.statistics != nullptr);
    if (p.error != nullptr) {
      return report(s, -5, p.error);
    }
    constexpr auto k = internals::apply<PointFunction>;
    if (!internals::ensureStackSize<k>(s,
                                       p.common_options.preserve_stack_size)) {
      return s.status;
    }
    auto device = 0;
    auto warp_size = 0;
    if (failed(s, MFRONT_GPB_GPU_API(GetDevice)(&device)) ||
        failed(s, MFRONT_GPB_GPU_API(DeviceGetAttribute)(
                      &warp_size, internals::warp_size_attribute, device))) {
      return s.status;
    }
    const auto ws = static_cast<unsigned int>(warp_size);
    const auto nthreads = ((threads_per_block + ws - 1) / ws) * ws;
    const auto nblocks = (n + nthreads - 1) / nthreads;
    // as for the material properties, the number of blocks must fit in the
    // grid of the device
    if (nblocks >
        static_cast<decltype(nblocks)>(std::numeric_limits<int>::max())) {
      return report(s, -5, "too many points for the grid of the device");
    }
    const auto stream = p.policy == nullptr ? Stream{} : p.policy->stream;
    const auto rdt_key = toOrderedKey(rdt0);
    if (e.statistics != nullptr) {
      // statistics given by the calling code: the result is given by the
      // report function of the library
      k<<<nblocks, nthreads, 0, stream>>>(f, n, e.statistics, rdt_key);
      if (failed(s, MFRONT_GPB_GPU_API(GetLastError)()) ||
          p.common_options.asynchronous) {
        return s.status;
      }
      failed(s, MFRONT_GPB_GPU_API(StreamSynchronize)(stream));
      return s.status;
    }
    // synchronous computations: the statistics are held by a block of
    // managed memory, initialized and read by the host. The device only
    // accesses this block if an event occurs (see internals::apply): no
    // data is transferred otherwise.
    const auto block = internals::getStatisticsBlocks().acquire(device, [] {
      AtomicBatchStatistics* p = nullptr;
      const auto e = MFRONT_GPB_GPU_API(MallocManaged)(
          reinterpret_cast<void**>(&p), sizeof(AtomicBatchStatistics));
      return e == MFRONT_GPB_GPU_API(Success) ? p : nullptr;
    });
    auto* const statistics = block.get();
    if (statistics == nullptr) {
      return report(s, -7, "unable to allocate the statistics of the batch");
    }
    mfront_gpb_initializeStatistics(statistics, rdt0);
    k<<<nblocks, nthreads, 0, stream>>>(f, n, statistics, rdt_key);
    if (failed(s, MFRONT_GPB_GPU_API(GetLastError)()) ||
        failed(s, MFRONT_GPB_GPU_API(StreamSynchronize)(stream))) {
      return s.status;
    }
    return report(s, getBatchStatistics(*statistics));
  }  // end of apply

  /*!
   * \return a copy of a data in the memory of the current device
   * \param[in] v: data
   */
  template <typename T>
  const T* copyToDevice(const T& v) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "the shared data must be trivially copyable");
    constexpr auto success = MFRONT_GPB_GPU_API(Success);
    T* p = nullptr;
    if (MFRONT_GPB_GPU_API(Malloc)(reinterpret_cast<void**>(&p), sizeof(T)) !=
        success) {
      throw std::runtime_error("memory allocation failed on the device");
    }
    if (MFRONT_GPB_GPU_API(Memcpy)(p, &v, sizeof(T),
                                   MFRONT_GPB_GPU_API(MemcpyHostToDevice)) !=
        success) {
      static_cast<void>(MFRONT_GPB_GPU_API(Free)(p));
      throw std::runtime_error("copy to the device failed");
    }
    return p;
  }  // end of copyToDevice

  /*!
   * \return a copy of a data shared by all the integration points (see the
   * `@SharedData` keyword) in the memory of the current device. A copy is
   * made at the first call on each device and is never freed.
   * \tparam get_shared_data: function returning the data
   */
  template <auto get_shared_data>
  const auto* getSharedDataCopy() {
    using T = std::remove_cvref_t<decltype(get_shared_data())>;
    static auto copies = SharedDataCopies<T, int>{};
    auto device = int{};
    if (MFRONT_GPB_GPU_API(GetDevice)(&device) != MFRONT_GPB_GPU_API(Success)) {
      throw std::runtime_error("unable to get the current device");
    }
    return copies.get(device, [] { return copyToDevice(get_shared_data()); });
  }  // end of getSharedDataCopy

  /*!
   * \return the primitive of the backend, as expected by
   * `mfront::gpb::integrateBatch` and the similar functions of
   * `MFront/GenericParallelBehaviour/Integrate.hxx`, called by the generated
   * functions
   * \param[in,out] s: output status
   * \param[in] threads_per_block: number of threads per block
   */
  inline auto makeApply(OutputStatus& s, const unsigned int threads_per_block) {
    return [&s, threads_per_block](const size_type n, const BatchExecution& e,
                                   const real rdt0, const auto& f) {
      return gpu::apply(s, n, e, rdt0, f, threads_per_block);
    };
  }  // end of makeApply

}  // end of namespace mfront::gpb::gpu

#undef MFRONT_GPB_GPU_API

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_GPU_HXX */
