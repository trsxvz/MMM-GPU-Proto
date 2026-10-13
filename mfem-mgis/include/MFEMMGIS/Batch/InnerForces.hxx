/*!
 * \file   include/MFEMMGIS/Batch/InnerForces.hxx
 * \brief  This file declares the computation of the inner forces of the
 * standard behaviour integrators (mechanics, heat transfer, micromorphic
 * damage) by kernels of MFEM, for all the elements of a material at once
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
      Context& ctx,
      mfem::Vector& f,
      const PartialQuadratureSpace& qs,
      const mgis::FieldView<const real>& s,
      const ElementsGeometry& g,
      const GradientsKind k) noexcept;
  /*!
   * \brief adds the inner forces of the elements of a material, as
   * `addInnerForces`, to values stored by element of the material
   * (`node + nd * (component + vdim * k)` for the `k`-th element of the
   * material) instead of by element of the mesh
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] f: values by element of the material
   * \param[in] qs: quadrature space of the material
   * \param[in] s: thermodynamic forces in the global frame
   * \param[in] g: geometric data
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addMaterialInnerForces(
      Context& ctx,
      mfem::Vector& f,
      const PartialQuadratureSpace& qs,
      const mgis::FieldView<const real>& s,
      const ElementsGeometry& g,
      const GradientsKind k) noexcept;

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_INNERFORCES_HXX */
