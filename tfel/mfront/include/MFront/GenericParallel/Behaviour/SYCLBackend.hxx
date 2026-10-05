/*!
 * \file   mfront/include/MFront/GenericParallel/Behaviour/SYCLBackend.hxx
 * \brief  This file declares the SYCL backend of the `generic-parallel`
 * behaviour interface
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_SYCLBACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_SYCLBACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/Behaviour/AbstractBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  /*!
   * \brief backend based on SYCL. The points are integrated on the device
   * associated with the queue given by the batch data.
   *
   * Options:
   *
   * - `data_transfer_policy`: `automatic` (or `auto`, default), the fields
   *   which are not accessible by the device are copied in device memory and
   *   the outputs are copied back, or `none`, all the fields are assumed to
   *   be accessible by the device. This option has the same meaning as for
   *   material properties.
   * - `device_stack_size`: minimal stack size, in bytes, of the threads of
   *   the device when the computations are performed on an NVIDIA GPU
   *   through the CUDA runtime (AdaptiveCpp).
   */
  struct MFRONT_VISIBILITY_EXPORT SYCLBackend final : AbstractBackend {
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    SYCLBackend(const tfel::utilities::DataMap&);
    //
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::string getSourceFileExtension() const override;
    void writeSpecificIncludes(std::ostream&) const override;
    void writeCall(std::ostream&,
                   std::string_view,
                   std::string_view,
                   std::string_view) const override;
    //! \brief destructor
    ~SYCLBackend() noexcept override;

   private:
    //! \brief if true, the fields not accessible by the device are copied
    bool automatic_data_transfers = true;
    //! \brief stack size of the threads of the device (0 if not specified)
    int device_stack_size = 0;
  };  // end of SYCLBackend

}  // end of namespace mfront::generic_parallel::behaviour

#endif /* LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_SYCLBACKEND_HXX */
