/*!
 * \file MFront/GenericParallel/MaterialProperty/SYCLBackend.hxx
 *
 * \brief  This file declares the `SYCLBackend` class
 * \author Thomas Helfer
 * \date   18/05/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_SYCLBACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_SYCLBACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/MaterialProperty/CxxProgrammingModelBackendBase.hxx"

namespace mfront::generic_parallel::material_property {

  //! \brief This class implements a backend for SYCL programming model
  struct MFRONT_VISIBILITY_EXPORT SYCLBackend : CxxProgrammingModelBackendBase {
    //! \brief default constructor
    SYCLBackend();
    /*!
     * \brief constructor
     *
     * \param[in] opts: options
     */
    SYCLBackend(const tfel::utilities::DataMap&);
    //
    [[nodiscard]] std::string getName() const override;
    void writeSpecificIncludesInSourceFile(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&) const override;
    //! \brief destructor
    [[nodiscard]] std::string getMemorySpace() const override;
    ~SYCLBackend() noexcept override;

   protected:
    //
    void writeExecutionPreparation(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&) const override;
    void writeKernelCall(std::ostream&,
                         const GenericParallelMaterialPropertyInterface&,
                         const MaterialPropertyDescription&,
                         const FileDescription&,
                         const bool) const override;
    bool handlesDataTransfer() const override;
    void writeDataTransfersToDevice(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&,
        const bool) const override;
    void writeDataTransfersToHost(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&,
        const bool) const override;
    /*!
     * \brief enumeration giving if data transfer are handled
     */
    enum struct DataTransferPolicy { AUTOMATIC, NONE };
    //
    DataTransferPolicy data_transfer_policy = DataTransferPolicy::AUTOMATIC;
    /*!
     * \brief minimal stack size of the threads of the device when the
     * computations are performed on an NVIDIA GPU through the CUDA runtime
     * (0 if not specified)
     */
    int device_stack_size = 0;
  };  // end of SYCLBackend

}  // end of namespace mfront::generic_parallel::material_property

#endif /* LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_SYCLBACKEND_HXX */
