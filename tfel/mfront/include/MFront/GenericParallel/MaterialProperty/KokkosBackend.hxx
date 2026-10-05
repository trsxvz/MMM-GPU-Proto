/*!
 * \file
 * mfront/include/MFront/GenericParallel/MaterialProperty/KokkosBackend.hxx
 * \brief  This file declares the backend of the `generic-parallel` material
 * property interface based on Kokkos
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_KOKKOSBACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_KOKKOSBACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/MaterialProperty/CxxProgrammingModelBackendBase.hxx"

namespace mfront::generic_parallel::material_property {

  /*!
   * \brief backend based on Kokkos.
   *
   * Options:
   *
   * - `execution_space`: execution space (`Kokkos::DefaultExecutionSpace`
   *   by default). The arguments and the output must be accessible from
   *   this execution space.
   *
   * Kokkos must be initialized by the calling code.
   */
  struct MFRONT_VISIBILITY_EXPORT KokkosBackend
      : CxxProgrammingModelBackendBase {
    //! \brief default constructor
    KokkosBackend();
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    KokkosBackend(const tfel::utilities::DataMap&);
    //
    [[nodiscard]] std::string getName() const override;
    void writeSpecificIncludesInSourceFile(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&) const override;
    //! \brief destructor
    ~KokkosBackend() noexcept override;

   protected:
    //
    [[nodiscard]] std::string getKernelSpecifiers() const override;
    void writeKernelCall(std::ostream&,
                         const GenericParallelMaterialPropertyInterface&,
                         const MaterialPropertyDescription&,
                         const FileDescription&,
                         const bool) const override;
    //! \brief execution space
    std::string execution_space = "Kokkos::DefaultExecutionSpace";
  };  // end of KokkosBackend

}  // end of namespace mfront::generic_parallel::material_property

#endif /* LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_KOKKOSBACKEND_HXX */
