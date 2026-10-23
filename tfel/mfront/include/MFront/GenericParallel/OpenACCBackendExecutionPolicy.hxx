/*!
 * \file mfront/include/MFront/GenericParallel/OpenACCBackendExecutionPolicy.hxx
 * \brief  This file declares the execution policy of the `openacc` backend
 * of the `generic-parallel` interfaces.
 * \date   22/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_OPENACCBACKENDEXECUTIONPOLICY_HXX
#define LIB_MFRONT_GENERICPARALLEL_OPENACCBACKENDEXECUTIONPOLICY_HXX

#include "MFront/GenericParallel/BackendExecutionPolicyBase.hxx"

namespace mfront::generic_parallel {

  /*!
   * \brief execution policy of the `openacc` backend: the computations are
   * performed synchronously, and only the common options are used
   */
  struct OpenACCBackendExecutionPolicy final : BackendExecutionPolicyBase {
    //! \return the identifier of the policy
    static constexpr const char* identifier() noexcept { return "openacc"; }
    /*!
     * \brief constructor
     * \param[in] opts: common options
     */
    explicit OpenACCBackendExecutionPolicy(
        const BackendExecutionCommonPolicy& opts = {}) noexcept
        : BackendExecutionPolicyBase(opts) {}
    const char* getIdentifier() const noexcept override { return identifier(); }
  };  // end of struct OpenACCBackendExecutionPolicy

  /*!
   * \return an execution policy of the `openacc` backend
   * \param[in] opts: common options
   */
  inline OpenACCBackendExecutionPolicy makeOpenACCBackendExecutionPolicy(
      const BackendExecutionCommonPolicy& opts = {}) noexcept {
    return OpenACCBackendExecutionPolicy{opts};
  }  // end of makeOpenACCBackendExecutionPolicy

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_OPENACCBACKENDEXECUTIONPOLICY_HXX */
