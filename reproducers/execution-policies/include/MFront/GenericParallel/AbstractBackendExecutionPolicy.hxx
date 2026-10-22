/*!
 * \file   MFront/GenericParallel/AbstractBackendExecutionPolicy.hxx
 * \brief  Base class of the execution policies, defined in the runtime
 * library. A generated function compares the identifier of the policy it
 * receives with the one of its backend, then uses `static_cast`.
 *
 * The macro `MFRONT_GENERICPARALLEL_IDENTIFIER_FIRST` declares
 * `getIdentifier` before the destructor. The position of a virtual function
 * in the table of virtual functions then does not depend on the C++ ABI:
 * the ABI of Visual Studio gives one entry to a virtual destructor, the
 * Itanium ABI (gcc, clang, MinGW) two.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_ABSTRACTBACKENDEXECUTIONPOLICY_HXX
#define LIB_MFRONT_GENERICPARALLEL_ABSTRACTBACKENDEXECUTIONPOLICY_HXX

#include "MFront/GenericParallel/Config.hxx"

namespace mfront::generic_parallel {

  //! \brief base class of the execution policies
  struct MFRONT_GENERICPARALLEL_RUNTIME_VISIBILITY
      AbstractBackendExecutionPolicy {
#ifdef MFRONT_GENERICPARALLEL_IDENTIFIER_FIRST
    //! \return the identifier of the policy, for instance "cuda"
    virtual const char* getIdentifier() const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractBackendExecutionPolicy();
#else
    //! \brief destructor
    virtual ~AbstractBackendExecutionPolicy();
    //! \return the identifier of the policy, for instance "cuda"
    virtual const char* getIdentifier() const noexcept = 0;
#endif
  };  // end of struct AbstractBackendExecutionPolicy

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_ABSTRACTBACKENDEXECUTIONPOLICY_HXX */
