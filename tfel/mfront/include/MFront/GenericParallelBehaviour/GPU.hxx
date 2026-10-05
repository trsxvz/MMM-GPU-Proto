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
#if defined(__HIPCC__) || defined(__HIP__)
#include <hip/hip_runtime.h>
#define MFRONT_GPB_GPU_API(f) hip##f
#else
#include <cuda_runtime.h>
#define MFRONT_GPB_GPU_API(f) cuda##f
#endif
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
     * \brief increase the stack size of the threads if the given kernel
     * requires it. The default stack size (1024 bytes with CUDA) is too
     * small for some behaviours and is not increased automatically.
     * \return the error code of the runtime
     * \param[in] k: kernel
     */
    template <typename Kernel>
    MFRONT_GPB_GPU_API(Error_t)
    ensureStackSize(Kernel k) {
      using Error = MFRONT_GPB_GPU_API(Error_t);
      // the stack size is a property of the context, shared by all threads
      static std::mutex m;
      auto a = MFRONT_GPB_GPU_API(FuncAttributes){};
      // the HIP runtime only provides the untyped overload
      if (const auto e = MFRONT_GPB_GPU_API(FuncGetAttributes)(
              &a, reinterpret_cast<const void*>(k));
          e != MFRONT_GPB_GPU_API(Success)) {
        return e;
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
      auto lock = std::lock_guard<std::mutex>{m};
      auto l = std::size_t{};
      if (const auto e = MFRONT_GPB_GPU_API(DeviceGetLimit)(
              &l, MFRONT_GPB_GPU_API(LimitStackSize));
          e != MFRONT_GPB_GPU_API(Success)) {
        return e;
      }
      if (required <= l) {
        return Error{MFRONT_GPB_GPU_API(Success)};
      }
      return MFRONT_GPB_GPU_API(DeviceSetLimit)(
          MFRONT_GPB_GPU_API(LimitStackSize), required);
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
   * The kernel is launched on the given queue. This function waits for its
   * completion.
   *
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] queue: stream (default stream if null)
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \param[in] threads_per_block: number of threads per block, rounded up to
   * a multiple of the warp size
   * \return the status
   */
  template <typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            void* const queue,
            const real rdt0,
            const PointFunction& f,
            const unsigned int threads_per_block) {
    using Stream = MFRONT_GPB_GPU_API(Stream_t);
    using Error = MFRONT_GPB_GPU_API(Error_t);
    if (n == 0) {
      return s.status;
    }
    // report a runtime error, only the first error is reported
    auto failed = [&s](const Error e) {
      if (e == MFRONT_GPB_GPU_API(Success)) {
        return false;
      }
      if (s.status != -3) {
        report(s, -3, MFRONT_GPB_GPU_API(GetErrorString)(e));
      }
      return true;
    };
    auto device = 0;
    auto warp_size = 0;
    if (failed(MFRONT_GPB_GPU_API(GetDevice)(&device)) ||
        failed(MFRONT_GPB_GPU_API(DeviceGetAttribute)(
            &warp_size, internals::warp_size_attribute, device))) {
      return s.status;
    }
    const auto ws = static_cast<unsigned int>(warp_size);
    const auto nthreads = ((threads_per_block + ws - 1) / ws) * ws;
    const auto nblocks = (n + nthreads - 1) / nthreads;
    const auto stream = static_cast<Stream>(queue);
    auto statistics = AtomicBatchStatistics::initial(rdt0);
    AtomicBatchStatistics* device_statistics = nullptr;
    if (failed(MFRONT_GPB_GPU_API(MallocAsync)(
            reinterpret_cast<void**>(&device_statistics),
            sizeof(AtomicBatchStatistics), stream))) {
      return s.status;
    }
    const auto launched = [&] {
      if (failed(MFRONT_GPB_GPU_API(MemcpyAsync)(
              device_statistics, &statistics, sizeof(statistics),
              MFRONT_GPB_GPU_API(MemcpyHostToDevice), stream))) {
        return false;
      }
      constexpr auto k = internals::apply<PointFunction>;
      if (failed(internals::ensureStackSize(k))) {
        return false;
      }
      k<<<nblocks, nthreads, 0, stream>>>(f, n, device_statistics,
                                          statistics.rdt);
      if (failed(MFRONT_GPB_GPU_API(GetLastError)())) {
        return false;
      }
      return !failed(MFRONT_GPB_GPU_API(MemcpyAsync)(
          &statistics, device_statistics, sizeof(statistics),
          MFRONT_GPB_GPU_API(MemcpyDeviceToHost), stream));
    }();
    const auto freed =
        !failed(MFRONT_GPB_GPU_API(FreeAsync)(device_statistics, stream));
    const auto synchronized =
        !failed(MFRONT_GPB_GPU_API(StreamSynchronize)(stream));
    if (!(launched && freed && synchronized)) {
      return s.status;
    }
    return report(s, statistics.get());
  }  // end of apply

  namespace internals {

    /*!
     * \return the primitive of the backend, as expected by
     * `mfront::gpb::integrateBatch` and the similar functions
     * \param[in,out] s: output status
     * \param[in] threads_per_block: number of threads per block
     */
    inline auto makeApply(OutputStatus& s,
                          const unsigned int threads_per_block) {
      return [&s, threads_per_block](const size_type n, void* const queue,
                                     const real rdt0, const auto& f) {
        return gpu::apply(s, n, queue, rdt0, f, threads_per_block);
      };
    }  // end of makeApply

  }  // end of namespace internals

  /*!
   * \brief integrate the behaviour at the points of a batch (see
   * `mfront::gpb::integrateBatch`).
   *
   * The kernel is launched on the queue given by the batch data. This
   * function waits for its completion.
   *
   * \param[out] s: output status
   * \param[in] d: batch data
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] threads_per_block: number of threads per block, rounded up to
   * a multiple of the warp size
   * \return the status
   */
  template <typename Behaviour, typename PointTraits, typename ParametersType>
  int integrate(OutputStatus& s,
                const BatchData& d,
                const ParametersType& parameters,
                const tfel::material::OutOfBoundsPolicy p,
                const unsigned int threads_per_block) {
    return integrateBatch<Behaviour, PointTraits>(
        s, d, parameters, p, internals::makeApply(s, threads_per_block));
  }  // end of integrate

  /*!
   * \brief execute an initialize function at the points of a batch (see
   * `mfront::gpb::executeInitializeFunctionOnBatch`).
   *
   * \param[out] s: output status
   * \param[in] d: batch data
   * \param[in] inputs: values of the inputs of the initialize function
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] threads_per_block: number of threads per block
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename InitializeFunction,
            size_type number_of_inputs,
            typename ParametersType>
  int executeInitializeFunction(OutputStatus& s,
                                const BatchData& d,
                                const ConstField& inputs,
                                const ParametersType& parameters,
                                const tfel::material::OutOfBoundsPolicy p,
                                const unsigned int threads_per_block) {
    return executeInitializeFunctionOnBatch<
        Behaviour, PointTraits, InitializeFunction, number_of_inputs>(
        s, d, inputs, parameters, p,
        internals::makeApply(s, threads_per_block));
  }  // end of executeInitializeFunction

  /*!
   * \brief execute a post-processing at the points of a batch (see
   * `mfront::gpb::executePostProcessingOnBatch`).
   *
   * \param[out] s: output status
   * \param[out] outputs: values of the outputs of the post-processing
   * \param[in] d: batch data
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] threads_per_block: number of threads per block
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename PostProcessing,
            bool use_initial_state,
            size_type number_of_outputs,
            typename ParametersType>
  int executePostProcessing(OutputStatus& s,
                            const Field& outputs,
                            const BatchData& d,
                            const ParametersType& parameters,
                            const tfel::material::OutOfBoundsPolicy p,
                            const unsigned int threads_per_block) {
    return executePostProcessingOnBatch<Behaviour, PointTraits, PostProcessing,
                                        use_initial_state, number_of_outputs>(
        s, outputs, d, parameters, p,
        internals::makeApply(s, threads_per_block));
  }  // end of executePostProcessing

  /*!
   * \brief rotate the values of the points of a batch of an orthotropic
   * behaviour (see `mfront::gpb::rotateBatch`).
   *
   * \param[out] s: output status
   * \param[in] d: rotation data
   * \param[in] threads_per_block: number of threads per block
   * \return the status
   */
  template <typename PointRotation>
  int rotate(OutputStatus& s,
             const RotationData& d,
             const unsigned int threads_per_block) {
    return rotateBatch<PointRotation>(
        s, d, internals::makeApply(s, threads_per_block));
  }  // end of rotate

}  // end of namespace mfront::gpb::gpu

#undef MFRONT_GPB_GPU_API

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_GPU_HXX */
