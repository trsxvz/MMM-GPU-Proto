/*!
 * \file mfront/include/MFront/GenericParallel/MaterialProperty/RAJABackend.hxx
 * \brief  This file declares the RAJA backend of the `generic-parallel`
 * material property interface
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_RAJABACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_RAJABACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/MaterialProperty/CxxProgrammingModelBackendBase.hxx"

namespace mfront::generic_parallel::material_property {

  /*!
   * \brief backend based on RAJA.
   *
   * Options:
   *
   * - `execution_policy`: execution policy of RAJA
   *   (`mfront::generic_parallel::raja::DefaultExecutionPolicy` by default,
   *   see `MFront/GenericParallel/RAJAExecutionPolicy.hxx`). The arguments
   *   and the output must be accessible from this execution policy.
   *
   * The execution policy of the backend (`RAJABackendExecutionPolicy`)
   * holds a resource of camp matching the execution policy of RAJA. With
   * the default policy, the default stream of the current device is used
   * with CUDA and HIP, the default resource of camp otherwise (see
   * `mfront::generic_parallel::raja::getDefaultResource`).
   * The execution policies of OpenMP target and SYCL are not supported.
   */
  struct MFRONT_VISIBILITY_EXPORT RAJABackend : CxxProgrammingModelBackendBase {
    //! \brief default constructor
    RAJABackend();
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    RAJABackend(const tfel::utilities::DataMap&);
    //
    [[nodiscard]] std::string getName() const override;
    void writeSpecificIncludesInSourceFile(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&) const override;
    [[nodiscard]] std::string getMemorySpace() const override;
    [[nodiscard]] std::string getExecutionPolicyType() const override;
    [[nodiscard]] std::string getExecutionPolicyHeader() const override;
    //! \brief destructor
    ~RAJABackend() noexcept override;

   protected:
    //! \brief the statuses of the points are merged with atomic operations
    [[nodiscard]] bool mergesPointStatusesInKernels() const override;
    [[nodiscard]] std::string getAtomicIncrement(
        const std::string_view) const override;
    [[nodiscard]] std::string getAtomicMinimum(
        const std::string_view, const std::string_view) const override;
    [[nodiscard]] std::string getKernelSpecifiers() const override;
    void writeKernelCall(std::ostream&,
                         const GenericParallelMaterialPropertyInterface&,
                         const MaterialPropertyDescription&,
                         const FileDescription&) const override;

   private:
    //! \brief execution policy
    std::string execution_policy =
        "mfront::generic_parallel::raja::DefaultExecutionPolicy";
  };  // end of RAJABackend

}  // end of namespace mfront::generic_parallel::material_property

#endif /* LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_RAJABACKEND_HXX */
