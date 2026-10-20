/*!
 * \file   mfront/include/MFront/GenericParallel/OpenACCStackSize.hxx
 * \brief  This file declares a function increasing the stack size of the
 * threads of the NVIDIA GPU used by OpenACC, shared by the openacc backends
 * of the `generic-parallel` interfaces.
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_OPENACCSTACKSIZE_HXX
#define LIB_MFRONT_GENERICPARALLEL_OPENACCSTACKSIZE_HXX

#include <cstddef>
#if defined(__NVCOMPILER_OPENACC_GPU)
#include <mutex>
#include <dlfcn.h>
#include <openacc.h>
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#endif /* defined(__NVCOMPILER_OPENACC_GPU) */

namespace mfront::generic_parallel::openacc {

  /*!
   * \brief increase the stack size of the threads of the current NVIDIA GPU
   * of OpenACC, if the code is compiled by `nvc++` for a GPU and if the
   * current device type of OpenACC is an NVIDIA GPU. Nothing is done
   * otherwise.
   *
   * The OpenACC runtime of `nvc++` uses the CUDA driver without the CUDA
   * runtime: the stack size is set through the driver, on the CUDA context
   * of OpenACC (`acc_get_current_cuda_context`), the device being
   * initialized if needed. The driver is looked up among the loaded
   * libraries, so that this function does not require the CUDA headers.
   * The `NV_ACC_CUDA_STACKSIZE` environment variable of `nvc++` only sets
   * the stack size when the device is initialized.
   *
   * \return true on success, false if the stack size could not be set or
   * is too small and must be preserved
   * \param[in] stack_size: minimal stack size, in bytes. Nothing is done if
   * null.
   * \param[in] preserve: if true, the stack size is only checked
   */
  inline bool setMinimalStackSize(const std::size_t stack_size,
                                  const bool preserve = false) noexcept {
#if defined(__NVCOMPILER_OPENACC_GPU)
    // `CU_LIMIT_STACK_SIZE` in `CUlimit` and `CUDA_SUCCESS` in `CUresult`
    constexpr int stack_size_limit = 0;
    constexpr int success = 0;
    using PushCurrent = int (*)(void*);
    using PopCurrent = int (*)(void**);
    using GetLimit = int (*)(std::size_t*, int);
    using SetLimit = int (*)(int, std::size_t);
    if ((stack_size == 0) || (acc_get_device_type() != acc_device_nvidia)) {
      return true;
    }
    auto* context = acc_get_current_cuda_context();
    if (context == nullptr) {
      // initialization of the current device of OpenACC
      acc_set_device_num(acc_get_device_num(acc_device_nvidia),
                         acc_device_nvidia);
      context = acc_get_current_cuda_context();
    }
    auto* const h = ::dlopen("libcuda.so.1", RTLD_LAZY | RTLD_NOLOAD);
    if ((context == nullptr) || (h == nullptr)) {
      if (h != nullptr) {
        ::dlclose(h);
      }
      return false;
    }
    const auto push =
        reinterpret_cast<PushCurrent>(::dlsym(h, "cuCtxPushCurrent_v2"));
    const auto pop =
        reinterpret_cast<PopCurrent>(::dlsym(h, "cuCtxPopCurrent_v2"));
    const auto get_limit =
        reinterpret_cast<GetLimit>(::dlsym(h, "cuCtxGetLimit"));
    const auto set_limit =
        reinterpret_cast<SetLimit>(::dlsym(h, "cuCtxSetLimit"));
    auto r = (push != nullptr) && (pop != nullptr) &&  //
             (get_limit != nullptr) && (set_limit != nullptr);
    if (r) {
      auto lock = std::lock_guard<std::mutex>{getStackSizeMutex()};
      r = push(context) == success;
      if (r) {
        auto current = std::size_t{};
        r = get_limit(&current, stack_size_limit) == success;
        if ((r) && (current < stack_size)) {
          r = (!preserve) &&
              (set_limit(stack_size_limit, stack_size) == success);
        }
        void* previous = nullptr;
        r = (pop(&previous) == success) && r;
      }
    }
    ::dlclose(h);
    return r;
#else
    static_cast<void>(stack_size);
    static_cast<void>(preserve);
    return true;
#endif
  }  // end of setMinimalStackSize

}  // end of namespace mfront::generic_parallel::openacc

#endif /* LIB_MFRONT_GENERICPARALLEL_OPENACCSTACKSIZE_HXX */
