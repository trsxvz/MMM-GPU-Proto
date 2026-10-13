/*!
 * \file
 * mfront/include/MFront/GenericParallel/MaterialProperty/OpenMPBackend.hxx
 * \brief  This file declares the backend of the `generic-parallel` material
 * property interface based on OpenMP
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_OPENMPBACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_OPENMPBACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/MaterialProperty/CxxProgrammingModelBackendBase.hxx"

namespace mfront::generic_parallel::material_property {

  /*!
   * \brief backend based on OpenMP.
   *
   * Options:
   *
   * - `offload`: if true (default), the computations are offloaded to the
   *   default device (`target` construct). Otherwise, they are performed by
   *   the threads of the host (`parallel` construct). If the computations
   *   are offloaded, the arguments and the output must be accessible from
   *   the device (unified or managed memory).
   * - `device_stack_size`: minimal stack size, in bytes, of the threads of the
   *   device when the computations are offloaded to an NVIDIA GPU through the
   *   CUDA runtime (ignored otherwise).
   */
  struct MFRONT_VISIBILITY_EXPORT OpenMPBackend
      : CxxProgrammingModelBackendBase {
    //! \brief default constructor
    OpenMPBackend();
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    OpenMPBackend(const tfel::utilities::DataMap&);
    //
    [[nodiscard]] std::string getName() const override;
    void writeSpecificIncludesInSourceFile(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&) const override;
    [[nodiscard]] std::string getMemorySpace() const override;
    //! \brief destructor
    ~OpenMPBackend() noexcept override;

   protected:
    //
    void writeKernelCall(std::ostream&,
                         const GenericParallelMaterialPropertyInterface&,
                         const MaterialPropertyDescription&,
                         const FileDescription&,
                         const bool) const override;
    //! \brief offload the computations to the default device
    bool offload = true;
    //! \brief minimal stack size of the threads of the device (0 if not set)
    int device_stack_size = 0;
  };  // end of OpenMPBackend

}  // end of namespace mfront::generic_parallel::material_property

#endif /* LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_OPENMPBACKEND_HXX */
