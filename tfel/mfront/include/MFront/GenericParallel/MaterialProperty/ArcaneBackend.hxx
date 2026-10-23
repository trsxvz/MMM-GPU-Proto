/*!
 * \file   mfront/include/MFront/GenericParallel/MaterialProperty/ArcaneBackend.hxx
 * \brief  This file declares the arcane backend of the `generic-parallel`
 * material property interface
 * \date   23/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_ARCANEBACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_ARCANEBACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/MaterialProperty/CxxProgrammingModelBackendBase.hxx"

namespace mfront::generic_parallel::material_property {

  /*!
   * \brief backend based on the accelerator API of Arccore, the core
   * component of the Arcane framework.
   *
   * The execution policy of the backend (`ArcaneBackendExecutionPolicy`)
   * holds a queue of Arccore, which is required. The arguments and the
   * output must be accessible by the device of the queue (`CUDA`, `HIP`,
   * `Thread` or `Sequential`), which is selected at runtime.
   *
   * This backend has no option.
   */
  struct MFRONT_VISIBILITY_EXPORT ArcaneBackend
      : CxxProgrammingModelBackendBase {
    //! \brief default constructor
    ArcaneBackend();
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    ArcaneBackend(const tfel::utilities::DataMap&);
    //
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::string getSourceFileExtension() const override;
    void writeSpecificIncludesInSourceFile(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&) const override;
    [[nodiscard]] std::string getMemorySpace() const override;
    [[nodiscard]] std::string getExecutionPolicyType() const override;
    [[nodiscard]] std::string getExecutionPolicyHeader() const override;
    //! \brief destructor
    ~ArcaneBackend() noexcept override;

   protected:
    //! \brief the statuses of the points are merged with atomic operations
    [[nodiscard]] bool mergesPointStatusesInKernels() const override;
    [[nodiscard]] std::string getAtomicIncrement(
        const std::string_view) const override;
    [[nodiscard]] std::string getAtomicMinimum(
        const std::string_view, const std::string_view) const override;
    [[nodiscard]] std::string getKernelSpecifiers() const override;
    void writeExecutionPreparation(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&) const override;
    void writeKernelCall(std::ostream&,
                         const GenericParallelMaterialPropertyInterface&,
                         const MaterialPropertyDescription&,
                         const FileDescription&) const override;
  };  // end of ArcaneBackend

}  // end of namespace mfront::generic_parallel::material_property

#endif /* LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_ARCANEBACKEND_HXX */
