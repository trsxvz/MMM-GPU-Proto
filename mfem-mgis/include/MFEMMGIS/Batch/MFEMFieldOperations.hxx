/*!
 * \file   include/MFEMMGIS/Batch/MFEMFieldOperations.hxx
 * \brief  This file declares the `MFEMFieldOperations` class
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_BATCH_MFEMFIELDOPERATIONS_HXX
#define LIB_MFEM_MGIS_BATCH_MFEMFIELDOPERATIONS_HXX

#include "MGIS/FieldOperations.hxx"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::batch {

  /*!
   * \brief copies and fills of views of MGIS (`mgis::FieldView`) by kernels
   * of MFEM (`mfem::forall`), for the update and the revert of the states of
   * a material (`mgis::behaviour::update`, `mgis::behaviour::revert`) whose
   * arrays are located on the device.
   *
   * The operations run on the device if the object was constructed with
   * `true`: the views must then be located in the memory of the device
   * (pointers given by `Read`, `Write`, `ReadWrite` of a `mfem::Memory`).
   * Otherwise, the views must be located on the host, and the operations
   * are the ones of `mgis::HostFieldOperations`. Nothing is detected. On a
   * GPU, the operations are asynchronous, on the stream of MFEM: an access
   * on the host must go through the copies of MFEM (`HostRead`, ...), which
   * wait for them.
   */
  struct MFEM_MGIS_EXPORT MFEMFieldOperations final : mgis::FieldOperations {
    /*!
     * \brief constructor
     * \param[in] d: if true, the views are located in the memory of the
     * device of MFEM (if a device is enabled), and the operations run on the
     * device. Otherwise, the views are located in the memory of the host,
     * and the operations run on the host, as `mgis::HostFieldOperations`.
     */
    explicit MFEMFieldOperations(const bool d) noexcept;
    //! \brief destructor
    ~MFEMFieldOperations() override;

   private:
    void do_copy(const mgis::FieldView<real>&,
                 const mgis::FieldView<const real>&) override;
    void do_fill(const mgis::FieldView<real>&, const real) override;
    //! \brief if true, the operations run on the device
    const bool use_device;
  };  // end of struct MFEMFieldOperations

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_MFEMFIELDOPERATIONS_HXX */
