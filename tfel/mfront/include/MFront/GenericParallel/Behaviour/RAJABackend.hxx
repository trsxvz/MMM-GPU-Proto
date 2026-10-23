/*!
 * \file   mfront/include/MFront/GenericParallel/Behaviour/RAJABackend.hxx
 * \brief  This file declares the RAJA backend of the `generic-parallel`
 * behaviour interface
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_RAJABACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_RAJABACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/Behaviour/AbstractBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  /*!
   * \brief backend based on RAJA.
   *
   * The generated library must be compiled with a compiler supporting the
   * execution policy (for instance `nvcc` with extended lambdas for the
   * policies of CUDA) and linked with RAJA.
   *
   * Options:
   *
   * - `execution_policy`: execution policy of RAJA
   *   (`mfront::generic_parallel::raja::DefaultExecutionPolicy` by default,
   *   see `MFront/GenericParallel/RAJAExecutionPolicy.hxx`), for instance
   *   `RAJA::cuda_exec<128>` or `RAJA::omp_parallel_for_exec`. The
   *   execution policies of OpenMP target and SYCL are not supported.
   *
   * The execution policy of the backend (`RAJABackendExecutionPolicy`)
   * holds a resource of camp matching the execution policy of RAJA. With
   * the default policy, the default stream of the current device is used
   * with CUDA and HIP, the default resource of camp otherwise (see
   * `mfront::generic_parallel::raja::getDefaultResource`).
   */
  struct MFRONT_VISIBILITY_EXPORT RAJABackend final : AbstractBackend {
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    RAJABackend(const tfel::utilities::DataMap&);
    //
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::string getSourceFileExtension() const override;
    void writeSpecificIncludes(std::ostream&) const override;
    [[nodiscard]] std::string getApply(std::string_view) const override;
    [[nodiscard]] std::string getMemorySpace() const override;
    [[nodiscard]] std::string getExecutionPolicyType() const override;
    [[nodiscard]] std::string getSharedDataCopy(
        std::string_view, std::string_view) const override;
    //! \brief destructor
    ~RAJABackend() noexcept override;

   private:
    //! \brief execution policy
    std::string execution_policy =
        "mfront::generic_parallel::raja::DefaultExecutionPolicy";
  };  // end of RAJABackend

}  // end of namespace mfront::generic_parallel::behaviour

#endif /* LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_RAJABACKEND_HXX */
