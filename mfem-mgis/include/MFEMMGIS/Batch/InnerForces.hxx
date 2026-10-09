/*!
 * \file   include/MFEMMGIS/Batch/InnerForces.hxx
 * \brief  This file declares the computation of the inner forces of the
 * standard mechanical behaviour integrators by kernels of MFEM, for all the
 * elements of a material at once
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_BATCH_INNERFORCES_HXX
#define LIB_MFEM_MGIS_BATCH_INNERFORCES_HXX

#include <map>
#include "mfem/linalg/vector.hpp"
#include "MGIS/FieldView.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"

namespace mfem {
  // forward declarations
  class DofToQuad;
  class GeometricFactors;
}  // end of namespace mfem

namespace mfem_mgis {
  // forward declaration
  struct PartialQuadratureSpace;
}  // end of namespace mfem_mgis

namespace mfem_mgis::batch {

  /*!
   * \brief geometric data of all the elements of the mesh for one
   * integration rule: gradients of the shape functions on the reference
   * element (`mfem::DofToQuad`, owned by the finite element) and jacobians
   * of the element transformations (`mfem::GeometricFactors`, owned by the
   * mesh)
   */
  struct MFEM_MGIS_EXPORT ElementsGeometry {
    //! \brief integration rule
    const mfem::IntegrationRule* rule = nullptr;
    //! \brief values and gradients of the shape functions at the points
    const mfem::DofToQuad* maps = nullptr;
    //! \brief jacobians and their determinants at the points
    const mfem::GeometricFactors* factors = nullptr;
    //! \brief weights of the integration points (in the memory of MFEM)
    mfem::Vector weights;
    //! \brief number of degrees of freedom per element and per component
    int nd = 0;
    //! \brief number of integration points per element
    int nq = 0;
    //! \brief dimension of the space
    int dim = 0;
    //! \brief number of elements
    int ne = 0;
  };

  //! \brief geometric data for the integration rules of several materials
  using ElementsGeometrySet =
      std::map<const mfem::IntegrationRule*, ElementsGeometry>;

  /*!
   * \brief computes the geometric data of all the elements of the mesh for
   * an integration rule
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] g: geometric data
   * \param[in] fes: finite element space. Its elements must have the same
   * geometry.
   * \param[in] ir: integration rule
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool computeElementsGeometry(
      Context&,
      ElementsGeometry&,
      const mfem::FiniteElementSpace&,
      const mfem::IntegrationRule&) noexcept;
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
