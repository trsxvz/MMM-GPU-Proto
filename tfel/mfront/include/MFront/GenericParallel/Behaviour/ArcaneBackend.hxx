/*!
 * \file   mfront/include/MFront/GenericParallel/Behaviour/ArcaneBackend.hxx
 * \brief  This file declares the arcane backend of the `generic-parallel`
 * behaviour interface
 * \date   23/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_ARCANEBACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_ARCANEBACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/Behaviour/AbstractBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  /*!
   * \brief backend based on the accelerator API of Arccore, the core
   * component of the Arcane framework.
   *
   * The generated source (`.cu` file) must be compiled with the compiler of
   * the accelerator mode of Arccore (for instance `clang` with `-x cuda` if
   * Arccore is built with CUDA, see the `cuda` compilation options of the
   * `arcane` configuration) and the library linked with the libraries of
   * Arccore.
   *
   * The execution policy of the backend (`ArcaneBackendExecutionPolicy`)
   * holds a queue of Arccore, which is required. The execution policy of
   * the queue (`CUDA`, `HIP`, `Thread` or `Sequential`) is selected at
   * runtime: a library compiled for CUDA computes on the GPU or on the
   * cores of the host, depending on the queue.
   *
   * This backend has no option.
   */
  struct MFRONT_VISIBILITY_EXPORT ArcaneBackend final : AbstractBackend {
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    ArcaneBackend(const tfel::utilities::DataMap&);
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
    ~ArcaneBackend() noexcept override;
  };  // end of ArcaneBackend

}  // end of namespace mfront::generic_parallel::behaviour

#endif /* LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_ARCANEBACKEND_HXX */
