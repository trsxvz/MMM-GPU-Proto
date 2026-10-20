/*!
 * \file   mfront/include/MFront/GenericParallel/RAJAExecutionPolicy.hxx
 * \brief  This file declares the default execution policies of the raja
 * backends of the `generic-parallel` interfaces, and the functions common to
 * these backends.
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_RAJAEXECUTIONPOLICY_HXX
#define LIB_MFRONT_GENERICPARALLEL_RAJAEXECUTIONPOLICY_HXX

#include <string>
#include <type_traits>
#include "RAJA/RAJA.hpp"
#if defined(RAJA_CUDA_ACTIVE) || defined(RAJA_HIP_ACTIVE)
#include "MFront/GenericParallel/GPUStackSize.hxx"
#endif

namespace mfront::generic_parallel::raja {

  /*!
   * \brief default execution policy: the asynchronous loop of the GPU if the
   * code is compiled for CUDA or HIP (128 threads per block, as the default
   * of the `cuda` and `hip` backends), the parallel loop of OpenMP if RAJA
   * is built with OpenMP, the sequential loop otherwise.
   */
#if defined(RAJA_CUDA_ACTIVE)
  using DefaultExecutionPolicy = RAJA::cuda_exec_async<128>;
#elif defined(RAJA_HIP_ACTIVE)
  using DefaultExecutionPolicy = RAJA::hip_exec_async<128>;
#elif defined(RAJA_ENABLE_OPENMP)
  using DefaultExecutionPolicy = RAJA::omp_parallel_for_exec;
#else
  using DefaultExecutionPolicy = RAJA::seq_exec;
#endif

  /*!
   * \brief default execution policy on the host: the parallel loop of
   * OpenMP if RAJA is built with OpenMP, the sequential loop otherwise.
   */
#if defined(RAJA_ENABLE_OPENMP)
  using DefaultHostExecutionPolicy = RAJA::omp_parallel_for_exec;
#else
  using DefaultHostExecutionPolicy = RAJA::seq_exec;
#endif

  /*!
   * \brief true if the execution policy is supported by the raja backends:
   * the policies of OpenMP target and SYCL are refused, since the resources
   * of camp associated with them can't allocate the managed memory holding
   * the statistics, and since the `raja` configuration of TFEL only
   * compiles for CUDA and for the host
   */
  template <typename ExecutionPolicy>
  inline constexpr bool isSupportedExecutionPolicy =
      !(RAJA::type_traits::is_target_openmp_policy<ExecutionPolicy>::value ||
        RAJA::type_traits::is_sycl_policy<ExecutionPolicy>::value);

  /*!
   * \brief true if the loops of the execution policy run on a GPU (CUDA,
   * HIP), which may not access the memory of the host
   */
  template <typename ExecutionPolicy>
  inline constexpr bool isDeviceExecutionPolicy =
      RAJA::type_traits::is_device_exec_policy<ExecutionPolicy>::value;

  /*!
   * \brief resource of camp associated with an execution policy
   * (`camp::resources::Cuda`, `camp::resources::Hip` or
   * `camp::resources::Host`, for instance). The queue given by the
   * execution options is a pointer to such a resource.
   */
  template <typename ExecutionPolicy>
  using Resource = RAJA::resources::resource_from_pol_t<ExecutionPolicy>;

  /*!
   * \return the resource used when the queue of the execution options is
   * null: the default stream of the current device with CUDA and HIP, as
   * for the `cuda` and `hip` backends (the default resource of camp is a
   * stream created by camp on the device 0), the default resource of camp
   * otherwise. The unsupported execution policies are refused at compile
   * time (see `isSupportedExecutionPolicy`).
   * \tparam ExecutionPolicy: execution policy
   */
  template <typename ExecutionPolicy>
  Resource<ExecutionPolicy> getDefaultResource() {
    static_assert(isSupportedExecutionPolicy<ExecutionPolicy>,
                  "the execution policies of OpenMP target and SYCL are not "
                  "supported by the raja backend");
    using ResourceType = Resource<ExecutionPolicy>;
#if defined(RAJA_CUDA_ACTIVE)
    if constexpr (std::is_same_v<ResourceType, camp::resources::Cuda>) {
      return camp::resources::Cuda::CudaFromStream(0);
    }
#endif
#if defined(RAJA_HIP_ACTIVE)
    if constexpr (std::is_same_v<ResourceType, camp::resources::Hip>) {
      return camp::resources::Hip::HipFromStream(0);
    }
#endif
    return ResourceType::get_default();
  }  // end of getDefaultResource

  /*!
   * \return the resource given by the execution options, or the default
   * resource (see `getDefaultResource`) if the queue is null
   * \tparam ExecutionPolicy: execution policy
   * \param[in] queue: queue of the execution options
   */
  template <typename ExecutionPolicy>
  Resource<ExecutionPolicy> getResource(const void* const queue) {
    if (queue == nullptr) {
      return getDefaultResource<ExecutionPolicy>();
    }
    return *static_cast<const Resource<ExecutionPolicy>*>(queue);
  }  // end of getResource

  /*!
   * \return the identifier of the device of a resource (CUDA and HIP), 0
   * otherwise. The memory allocated by the resources (statistics, copies of
   * the shared data) is distinguished by device.
   * \param[in] r: resource
   */
  template <typename ResourceType>
  int getDevice(const ResourceType& r) {
    static_cast<void>(r);
#if defined(RAJA_CUDA_ACTIVE)
    if constexpr (std::is_same_v<ResourceType, camp::resources::Cuda>) {
      return r.get_device();
    }
#endif
#if defined(RAJA_HIP_ACTIVE)
    if constexpr (std::is_same_v<ResourceType, camp::resources::Hip>) {
      return r.get_device();
    }
#endif
    return 0;
  }  // end of getDevice

#if defined(RAJA_CUDA_ACTIVE) || defined(RAJA_HIP_ACTIVE)

  /*!
   * \brief kernel which is never launched: its attributes give the stack
   * size required to apply a function at one point. RAJA does not increase
   * the stack size of the threads, which is too small by default for some
   * behaviours.
   */
  template <typename Function>
  __global__ void stackSizeProbe(const Function f) {
    f(0);
  }  // end of stackSizeProbe

#endif /* defined(RAJA_CUDA_ACTIVE) || defined(RAJA_HIP_ACTIVE) */

  /*!
   * \brief increase the stack size of the threads of the device if the
   * given function, applied at one point with the given execution policy,
   * requires it. Nothing is done if the execution policy is not one of CUDA
   * or HIP.
   * \return an error message, empty on success
   * \tparam ExecutionPolicy: execution policy
   * \param[in] f: function applied at one point
   * \param[in] preserve: if true, the stack size is not modified, and an
   * error is returned if it is too small
   */
  template <typename ExecutionPolicy, typename Function>
  std::string ensureStackSize(const Function& f, const bool preserve) {
    static_cast<void>(f);
    static_cast<void>(preserve);
#if defined(RAJA_CUDA_ACTIVE) || defined(RAJA_HIP_ACTIVE)
    if constexpr (RAJA::type_traits::is_device_exec_policy<
                      ExecutionPolicy>::value) {
      const auto r = gpu::ensureStackSize<stackSizeProbe<Function>>(preserve);
      if (!r.succeeded()) {
        return r.getErrorMessage();
      }
    }
#endif
    return {};
  }  // end of ensureStackSize

}  // end of namespace mfront::generic_parallel::raja

#endif /* LIB_MFRONT_GENERICPARALLEL_RAJAEXECUTIONPOLICY_HXX */
