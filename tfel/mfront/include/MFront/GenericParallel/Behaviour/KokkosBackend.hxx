/*!
 * \file   mfront/include/MFront/GenericParallel/Behaviour/KokkosBackend.hxx
 * \brief  This file declares the Kokkos backend of the `generic-parallel`
 * behaviour interface
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_KOKKOSBACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_KOKKOSBACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/Behaviour/AbstractBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  /*!
   * \brief backend based on Kokkos.
   *
   * The generated library must be compiled with the compiler used to build
   * Kokkos (for instance `nvcc_wrapper`) and linked with a shared build of
   * Kokkos, which is shared with the calling code.
   *
   * Options:
   *
   * - `execution_space`: execution space (`Kokkos::DefaultExecutionSpace` by
   *   default), for instance `Kokkos::DefaultHostExecutionSpace`.
   */
  struct MFRONT_VISIBILITY_EXPORT KokkosBackend final : AbstractBackend {
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    KokkosBackend(const tfel::utilities::DataMap&);
    //
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::string getSourceFileExtension() const override;
    void writeSpecificIncludes(std::ostream&) const override;
    [[nodiscard]] std::string getApply(std::string_view) const override;
    [[nodiscard]] std::string getSharedDataCopy(
        std::string_view, std::string_view) const override;
    //! \brief destructor
    ~KokkosBackend() noexcept override;

   private:
    //! \brief execution space
    std::string execution_space = "Kokkos::DefaultExecutionSpace";
  };  // end of KokkosBackend

}  // end of namespace mfront::generic_parallel::behaviour

#endif /* LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_KOKKOSBACKEND_HXX */
