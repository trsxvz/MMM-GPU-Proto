/*!
 * \file   MFront/GenericParallel/ParallelSTLBackendExecutionPolicy.hxx
 * \brief  Policy of the `stlpar` backend (common options only). Defined in
 * this header only.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_PARALLELSTLBACKENDEXECUTIONPOLICY_HXX
#define LIB_MFRONT_GENERICPARALLEL_PARALLELSTLBACKENDEXECUTIONPOLICY_HXX

#include "MFront/GenericParallel/BackendExecutionPolicyBase.hxx"

namespace mfront::generic_parallel {

  //! \brief policy of the `stlpar` backend
  struct ParallelSTLBackendExecutionPolicy final : BackendExecutionPolicyBase {
    static const char* identifier() noexcept { return "stlpar"; }
    using BackendExecutionPolicyBase::BackendExecutionPolicyBase;
    const char* getIdentifier() const noexcept override { return identifier(); }
  };  // end of struct ParallelSTLBackendExecutionPolicy

  inline ParallelSTLBackendExecutionPolicy
  makeParallelSTLBackendExecutionPolicy(
      const BackendExecutionCommonPolicy& opts = {}) noexcept {
    return ParallelSTLBackendExecutionPolicy{opts};
  }

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_PARALLELSTLBACKENDEXECUTIONPOLICY_HXX */
