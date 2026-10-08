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
#include <stdexcept>
#include <string>
#include <type_traits>
#if defined(__HIPCC__) || defined(__HIP__)
#include <hip/hip_runtime.h>
#define MFRONT_GPB_GPU_API(f) hip##f
#else
#include <cuda_runtime.h>
#define MFRONT_GPB_GPU_API(f) cuda##f
#endif
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::gpu {

  namespace internals {

    //! \brief device attribute giving the warp size
#if defined(__HIPCC__) || defined(__HIP__)
    inline constexpr auto warp_size_attribute = hipDeviceAttributeWarpSize;
#else
    inline constexpr auto warp_size_attribute = cudaDevAttrWarpSize;
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
      if (s.status != -3) {
        report(s, -3, MFRONT_GPB_GPU_API(GetErrorString)(e));
      }
      return true;
    }  // end of failed

    /*!
     * \brief increase the stack size of the threads if the given kernel
     * requires it. The default stack size (1024 bytes with CUDA) is too
     * small for some behaviours and is not increased automatically.
     * \return true on success
     * \param[in,out] s: output status
     * \param[in] k: kernel
     * \param[in] preserve: if true, the stack size is not modified and an
     * error is reported if it is too small
     */
    template <typename Kernel>
    bool ensureStackSize(OutputStatus& s, Kernel k, const bool preserve) {
      auto a = MFRONT_GPB_GPU_API(FuncAttributes){};
      // the HIP runtime only provides the untyped overload
      if (failed(s, MFRONT_GPB_GPU_API(FuncGetAttributes)(
                        &a, reinterpret_cast<const void*>(k)))) {
        return false;
      }
      // the static estimate of the local memory of the kernel may be
      // slightly lower than the actual stack usage (by a few tens of bytes
      // for a kernel calling the pow function of libdevice), so a margin is
      // added
      constexpr auto margin = std::size_t{512};
      constexpr auto granularity = std::size_t{1024};
      const auto required =
          ((a.localSizeBytes + margin + granularity - 1) / granularity) *
          granularity;
      auto lock = std::lock_guard<std::mutex>{
          mfront::generic_parallel::getStackSizeMutex()};
      auto l = std::size_t{};
      if (failed(s, MFRONT_GPB_GPU_API(DeviceGetLimit)(
                        &l, MFRONT_GPB_GPU_API(LimitStackSize)))) {
        return false;
      }
      if (required <= l) {
        return true;
      }
      if (preserve) {
        report(s, -3,
               "the stack size of the threads of the device must be at "
               "least " +
                   std::to_string(required) +
                   " bytes, but it is preserved by the execution options");
        return false;
      }
      return !failed(s, MFRONT_GPB_GPU_API(DeviceSetLimit)(
                            MFRONT_GPB_GPU_API(LimitStackSize), required));
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
      auto unreliable_results = 0u;
      auto rdt = rdt0;
      if (i < n) {
        const auto r = f(i);
        failures = static_cast<unsigned int>(r.number_of_failures);
        unreliable_results =
            static_cast<unsigned int>(r.number_of_unreliable_results);
        rdt = toOrderedKey(r.rdt);
      }
      // the atomic operations are only performed if a point of the warp
      // failed or decreased the time step scaling factor
      failures = warpSum(failures);
      unreliable_results = warpSum(unreliable_results);
      rdt = warpMin(rdt);
      if (threadIdx.x % warpSize == 0) {
        if (failures != 0) {
          atomicAdd(&(s->number_of_failures),
                    static_cast<unsigned long long>(failures));
        }
        if (unreliable_results != 0) {
          atomicAdd(&(s->number_of_unreliable_results),
                    static_cast<unsigned long long>(unreliable_results));
        }
        if (rdt < rdt0) {
          atomicMin(&(s->rdt), rdt);
        }
      }
    }  // end of apply

  }  // end of namespace internals

  /*!
   * \brief apply a function at the points of a batch, one thread per point,
   * and report the statistics of the batch.
   *
   * The kernel is launched on the queue given by the execution options. If
   * the execution options give statistics, this function returns once the
   * kernel is launched. Otherwise, it waits for its completion.
   *
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] o: execution options. The queue is a stream (default stream
   * if null).
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \param[in] threads_per_block: number of threads per block, rounded up to
   * a multiple of the warp size
   * \return the status
   */
  template <typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            const ExecutionOptions& o,
            const real rdt0,
            const PointFunction& f,
            const unsigned int threads_per_block) {
    using Stream = MFRONT_GPB_GPU_API(Stream_t);
    using internals::failed;
    if (n == 0) {
      return s.status;
    }
    constexpr auto k = internals::apply<PointFunction>;
    if (!internals::ensureStackSize(s, k, o.preserve_stack_size != 0)) {
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
    const auto stream = static_cast<Stream>(o.queue);
    const auto rdt_key = toOrderedKey(rdt0);
    if (o.statistics != nullptr) {
      k<<<nblocks, nthreads, 0, stream>>>(f, n, o.statistics, rdt_key);
      failed(s, MFRONT_GPB_GPU_API(GetLastError)());
      return s.status;
    }
    auto statistics = makeAtomicBatchStatistics(rdt0);
    AtomicBatchStatistics* device_statistics = nullptr;
    if (failed(s, MFRONT_GPB_GPU_API(MallocAsync)(
                      reinterpret_cast<void**>(&device_statistics),
                      sizeof(AtomicBatchStatistics), stream))) {
      return s.status;
    }
    const auto launched = [&] {
      if (failed(s, MFRONT_GPB_GPU_API(MemcpyAsync)(
                        device_statistics, &statistics, sizeof(statistics),
                        MFRONT_GPB_GPU_API(MemcpyHostToDevice), stream))) {
        return false;
      }
      k<<<nblocks, nthreads, 0, stream>>>(f, n, device_statistics, rdt_key);
      if (failed(s, MFRONT_GPB_GPU_API(GetLastError)())) {
        return false;
      }
      return !failed(s, MFRONT_GPB_GPU_API(MemcpyAsync)(
                            &statistics, device_statistics, sizeof(statistics),
                            MFRONT_GPB_GPU_API(MemcpyDeviceToHost), stream));
    }();
    const auto freed =
        !failed(s, MFRONT_GPB_GPU_API(FreeAsync)(device_statistics, stream));
    const auto synchronized =
        !failed(s, MFRONT_GPB_GPU_API(StreamSynchronize)(stream));
    if (!(launched && freed && synchronized)) {
      return s.status;
    }
    return report(s, getBatchStatistics(statistics));
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
    return [&s, threads_per_block](const size_type n, const ExecutionOptions& o,
                                   const real rdt0, const auto& f) {
      return gpu::apply(s, n, o, rdt0, f, threads_per_block);
    };
  }  // end of makeApply

}  // end of namespace mfront::gpb::gpu

#undef MFRONT_GPB_GPU_API

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_GPU_HXX */
