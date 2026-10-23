/*!
 * \file   mfront/include/MFront/GenericParallel/ArcaneBackendExecutionPolicy.hxx
 * \brief  This file declares the execution policy of the `arcane` backend
 * of the `generic-parallel` interfaces, which holds a queue of the
 * accelerator API of Arccore.
 * \date   23/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_ARCANEBACKENDEXECUTIONPOLICY_HXX
#define LIB_MFRONT_GENERICPARALLEL_ARCANEBACKENDEXECUTIONPOLICY_HXX

#include <arccore/base/ArccoreGlobal.h>
#include <arccore/common/accelerator/RunQueue.h>
#include "MFront/GenericParallel/BackendExecutionPolicyBase.hxx"

namespace mfront::generic_parallel {

  /*!
   * \brief execution policy of the `arcane` backend. The `arcane` backend
   * has no default policy: the queue gives the execution policy of Arccore
   * (`CUDA`, `HIP`, `Thread` or `Sequential`), which is only known at
   * runtime, so the queue is required.
   *
   * The identifier only depends on the accelerator mode of Arccore, so that
   * a calling code compiled by a C++ compiler and a library compiled by a
   * CUDA compiler agree on it. A library compiled by a CUDA (HIP) compiler
   * accepts the CUDA (HIP), `Thread` and `Sequential` queues.
   */
  struct ArcaneBackendExecutionPolicy final : BackendExecutionPolicyBase {
    //! \return the identifier of the policy
    static constexpr const char* identifier() noexcept {
#if defined(ARCCORE_HAS_CUDA)
      return "arcane:cuda";
#elif defined(ARCCORE_HAS_HIP)
      return "arcane:hip";
#elif defined(ARCCORE_HAS_SYCL)
#error "the SYCL accelerator mode of Arccore is not supported"
#else
      return "arcane:host";
#endif
    }
    /*!
     * \brief constructor
     * \param[in] q: queue on which the commands are launched. The queue has
     * a reference semantics: the policy refers to the queue of the calling
     * code.
     * \param[in] opts: common options
     */
    ArcaneBackendExecutionPolicy(const ::Arcane::Accelerator::RunQueue& q,
                                 const BackendExecutionCommonPolicy& opts = {})
        : BackendExecutionPolicyBase(opts), queue(q) {}
    const char* getIdentifier() const noexcept override { return identifier(); }
    //! \brief queue on which the commands are launched
    ::Arcane::Accelerator::RunQueue queue;
  };  // end of struct ArcaneBackendExecutionPolicy

  /*!
   * \return an execution policy of the `arcane` backend
   * \param[in] q: queue on which the commands are launched (required)
   * \param[in] opts: common options
   */
  inline ArcaneBackendExecutionPolicy makeArcaneBackendExecutionPolicy(
      const ::Arcane::Accelerator::RunQueue& q,
      const BackendExecutionCommonPolicy& opts = {}) {
    return ArcaneBackendExecutionPolicy{q, opts};
  }  // end of makeArcaneBackendExecutionPolicy

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_ARCANEBACKENDEXECUTIONPOLICY_HXX */
