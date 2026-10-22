/*!
 * \file   MFront/GenericParallel/SyclBackendExecutionPolicy.hxx
 * \brief  Policy of the `sycl` backend, which holds a `sycl::queue`. Its
 * identifier includes the implementation of SYCL. Defined in this header
 * only.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_SYCLBACKENDEXECUTIONPOLICY_HXX
#define LIB_MFRONT_GENERICPARALLEL_SYCLBACKENDEXECUTIONPOLICY_HXX

#include "StandIns.hxx"
#include "MFront/GenericParallel/BackendExecutionPolicyBase.hxx"

namespace mfront::generic_parallel {

  //! \brief policy of the `sycl` backend
  struct SyclBackendExecutionPolicy final : BackendExecutionPolicyBase {
    static const char* identifier() noexcept {
      return "sycl:" MFRONT_GENERICPARALLEL_SYCL_IMPLEMENTATION;
    }
    SyclBackendExecutionPolicy(const sycl::queue& q,
                               const BackendExecutionCommonPolicy& opts)
        : BackendExecutionPolicyBase(opts), queue(q) {}
    const char* getIdentifier() const noexcept override { return identifier(); }
    //! \brief queue
    sycl::queue queue;
  };  // end of struct SyclBackendExecutionPolicy

  //! \brief the queue is required: the arrays belong to its context
  inline SyclBackendExecutionPolicy makeSyclBackendExecutionPolicy(
      const sycl::queue& q, const BackendExecutionCommonPolicy& opts = {}) {
    return SyclBackendExecutionPolicy{q, opts};
  }

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_SYCLBACKENDEXECUTIONPOLICY_HXX */
