/*!
 * \file   MFront/GenericParallel/BackendExecutionPolicyBase.hxx
 * \brief  Options common to most backends. Defined in this header only.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_BACKENDEXECUTIONPOLICYBASE_HXX
#define LIB_MFRONT_GENERICPARALLEL_BACKENDEXECUTIONPOLICYBASE_HXX

#include <cstring>
#include "MFront/GenericParallel/AbstractBackendExecutionPolicy.hxx"

namespace mfront::generic_parallel {

  //! \brief options common to most backends
  struct BackendExecutionCommonPolicy {
    //! \brief do not increase the stack size of the device
    bool preserve_stack_size = false;
    //! \brief return before the end of the computations
    bool asynchronous = false;
  };  // end of struct BackendExecutionCommonPolicy

  //! \brief base class of the policies of the backends
  struct BackendExecutionPolicyBase : AbstractBackendExecutionPolicy,
                                      BackendExecutionCommonPolicy {
    explicit BackendExecutionPolicyBase(
        const BackendExecutionCommonPolicy& opts) noexcept
        : BackendExecutionCommonPolicy(opts) {}
  };  // end of struct BackendExecutionPolicyBase

  /*!
   * \return the policy given by a generated function, or a null pointer if
   * its identifier is not the one of `Policy`
   * \param[in] p: policy given by the caller
   */
  template <typename Policy>
  const Policy* getBackendExecutionPolicy(
      const AbstractBackendExecutionPolicy& p) noexcept {
    if (std::strcmp(p.getIdentifier(), Policy::identifier()) != 0) {
      return nullptr;
    }
    // valid: the bases are not virtual
    return static_cast<const Policy*>(&p);
  }  // end of getBackendExecutionPolicy

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_BACKENDEXECUTIONPOLICYBASE_HXX */
