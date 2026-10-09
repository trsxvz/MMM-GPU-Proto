/*!
 * \file   include/MFEMMGIS/Batch/ElementsGeometry.hxx
 * \brief  This file declares the geometric data of the elements used by the
 * kernels of the batch path
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_BATCH_ELEMENTSGEOMETRY_HXX
#define LIB_MFEM_MGIS_BATCH_ELEMENTSGEOMETRY_HXX

#include <map>
#include <memory>
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem {
  // forward declarations
  class DofToQuad;
  class GeometricFactors;
  class IntegrationRule;
  class FiniteElementSpace;
}  // end of namespace mfem

namespace mfem_mgis::batch {

  /*!
   * \brief geometric data of all the elements of the mesh for one
   * integration rule: gradients of the shape functions on the reference
   * element (`mfem::DofToQuad`, owned by the finite element) and jacobians
   * of the element transformations (`mfem::GeometricFactors`)
   */
  struct MFEM_MGIS_EXPORT ElementsGeometry {
    //! \brief integration rule
    const mfem::IntegrationRule* rule = nullptr;
    //! \brief values and gradients of the shape functions at the points
    const mfem::DofToQuad* maps = nullptr;
    /*!
     * \brief jacobians and their determinants at the points: owned by the
     * mesh if it has nodes, by `owned_factors` otherwise
     */
    const mfem::GeometricFactors* factors = nullptr;
    /*!
     * \brief jacobians and their determinants computed from a copy of the
     * vertices of a mesh without nodes (see `computeElementsGeometry`)
     */
    std::shared_ptr<const mfem::GeometricFactors> owned_factors;
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
   *
   * For a mesh without nodes, the jacobians are computed from a copy of the
   * vertices: `mfem::Mesh::GetGeometricFactors` would give nodes to the
   * mesh (`mfem::Mesh::EnsureNodes`), and the element transformations of
   * all the computations which follow would use them, which makes the
   * computations element by element about 50% slower.
   *
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

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_ELEMENTSGEOMETRY_HXX */
