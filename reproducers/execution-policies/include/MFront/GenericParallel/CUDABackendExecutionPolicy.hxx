/*!
 * \file   MFront/GenericParallel/CUDABackendExecutionPolicy.hxx
 * \brief  Policy of the `cuda` backend. Defined in this header only.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_CUDABACKENDEXECUTIONPOLICY_HXX
#define LIB_MFRONT_GENERICPARALLEL_CUDABACKENDEXECUTIONPOLICY_HXX

#include "StandIns.hxx"
#include "MFront/GenericParallel/BackendExecutionPolicyBase.hxx"

namespace mfront::generic_parallel {

  //! \brief options of the `cuda` backend
  struct CUDABackendExecutionOptions {
    BackendExecutionCommonPolicy common_options;
    //! \brief stream, the default stream by default
    cudaStream_t stream = nullptr;
  };  // end of struct CUDABackendExecutionOptions

  //! \brief policy of the `cuda` backend
  struct CUDABackendExecutionPolicy final : BackendExecutionPolicyBase {
    static const char* identifier() noexcept { return "cuda"; }
    explicit CUDABackendExecutionPolicy(
        const CUDABackendExecutionOptions& opts) noexcept
        : BackendExecutionPolicyBase(opts.common_options),
          stream(opts.stream) {}
    const char* getIdentifier() const noexcept override { return identifier(); }
    //! \brief stream
    cudaStream_t stream = nullptr;
  };  // end of struct CUDABackendExecutionPolicy

  inline CUDABackendExecutionPolicy makeCUDABackendExecutionPolicy(
      const CUDABackendExecutionOptions& opts = {}) noexcept {
    return CUDABackendExecutionPolicy{opts};
  }

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_CUDABACKENDEXECUTIONPOLICY_HXX */
