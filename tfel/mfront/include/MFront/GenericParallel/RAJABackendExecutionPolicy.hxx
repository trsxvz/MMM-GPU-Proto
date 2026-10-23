/*!
 * \file   mfront/include/MFront/GenericParallel/RAJABackendExecutionPolicy.hxx
 * \brief  This file declares the execution policy of the `raja` backend of
 * the `generic-parallel` interfaces: a class template on the resource of
 * camp associated with the execution policy of RAJA selected by the
 * backend, whose identifier is made of the name of the backend and of the
 * platform of the resource.
 * \date   22/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_RAJABACKENDEXECUTIONPOLICY_HXX
#define LIB_MFRONT_GENERICPARALLEL_RAJABACKENDEXECUTIONPOLICY_HXX

#include "MFront/GenericParallel/RAJAExecutionPolicy.hxx"
#include "MFront/GenericParallel/BackendExecutionPolicyBase.hxx"

namespace mfront::generic_parallel {

  namespace raja {

    /*!
     * \brief identifier of the policies using a resource of camp, null for
     * the resources which are not supported (see
     * `isSupportedExecutionPolicy`)
     */
    template <typename ResourceType>
    inline constexpr const char* policy_identifier = nullptr;
    //! \brief specialisation for the host
    template <>
    inline constexpr const char* policy_identifier<camp::resources::Host> =
        "raja:host";
#if defined(RAJA_CUDA_ACTIVE)
    //! \brief specialisation for CUDA
    template <>
    inline constexpr const char* policy_identifier<camp::resources::Cuda> =
        "raja:cuda";
#endif
#if defined(RAJA_HIP_ACTIVE)
    //! \brief specialisation for HIP
    template <>
    inline constexpr const char* policy_identifier<camp::resources::Hip> =
        "raja:hip";
#endif

  }  // end of namespace raja

  //! \brief options of the execution policy of the `raja` backend
  template <typename ResourceType>
  struct RAJABackendExecutionOptions {
    //! \brief common options
    BackendExecutionCommonPolicy common_options;
    /*!
     * \brief resource on which the loops are executed, the default resource
     * by default (see `raja::getDefaultResourceOfType`)
     */
    ResourceType resource = raja::getDefaultResourceOfType<ResourceType>();
  };  // end of struct RAJABackendExecutionOptions

  /*!
   * \brief execution policy of the `raja` backend
   * \tparam ResourceType: resource of camp associated with the execution
   * policy of RAJA selected by the `execution_policy` option of the backend
   * (see `raja::Resource`)
   */
  template <typename ResourceType>
  struct RAJABackendExecutionPolicy final : BackendExecutionPolicyBase {
    static_assert(raja::policy_identifier<ResourceType> != nullptr,
                  "unsupported resource");
    //! \return the identifier of the policy, for instance "raja:cuda"
    static constexpr const char* identifier() noexcept {
      return raja::policy_identifier<ResourceType>;
    }
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    explicit RAJABackendExecutionPolicy(
        const RAJABackendExecutionOptions<ResourceType>& opts = {})
        : BackendExecutionPolicyBase(opts.common_options),
          resource(opts.resource) {}
    const char* getIdentifier() const noexcept override { return identifier(); }
    //! \brief resource on which the loops are executed
    ResourceType resource;
  };  // end of struct RAJABackendExecutionPolicy

  /*!
   * \return an execution policy of the `raja` backend
   * \param[in] opts: options
   */
  template <typename ResourceType>
  RAJABackendExecutionPolicy<ResourceType> makeRAJABackendExecutionPolicy(
      const RAJABackendExecutionOptions<ResourceType>& opts = {}) {
    return RAJABackendExecutionPolicy<ResourceType>{opts};
  }  // end of makeRAJABackendExecutionPolicy

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_RAJABACKENDEXECUTIONPOLICY_HXX */
