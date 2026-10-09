/*!
 * \file   include/MFEMMGIS/Batch/InnerForces.hxx
 * \brief  This file declares the computation of the inner forces of the
 * standard mechanical behaviour integrators by kernels of MFEM, for all the
 * elements of a material at once
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_BATCH_INNERFORCES_HXX
#define LIB_MFEM_MGIS_BATCH_INNERFORCES_HXX

#include "mfem/linalg/vector.hpp"
#include "MGIS/FieldView.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"

namespace mfem_mgis {
  // forward declaration
  struct PartialQuadratureSpace;
}  // end of namespace mfem_mgis

namespace mfem_mgis::batch {

  /*!
   * \brief adds the inner forces of the elements of a material, `Bᵀ σ w det
   * J` summed over the integration points, to a vector of values element by
   * element, in the layout of the element restriction of MFEM (`nd x vdim x
   * ne`, native ordering of the degrees of freedom), on the device if a
   * device is enabled
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] f: values element by element
   * \param[in] qs: quadrature space of the material. All its elements must
   * use the integration rule of the geometric data.
   * \param[in] s: thermodynamic forces in the global frame at all the
   * points of the material, located on the device if a device is enabled
   * \param[in] g: geometric data
   * \param[in] k: kind of the gradients, which defines the operator `B`
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addInnerForces(
      Context&,
      mfem::Vector&,
      const PartialQuadratureSpace&,
      const mgis::FieldView<const real>&,
      const ElementsGeometry&,
      const GradientsKind) noexcept;

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_INNERFORCES_HXX */
