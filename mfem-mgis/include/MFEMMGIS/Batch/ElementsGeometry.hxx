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
  class IntegrationRule;
  class FiniteElementSpace;
}  // end of namespace mfem

namespace mfem_mgis::batch {

  /*!
   * \brief inverses of the jacobians of the element transformations and
   * determinants of the jacobians at the integration points (see
   * `computeElementsGeometry`). The kernels only use the jacobians through
   * their inverses, which are computed once.
   */
  struct MFEM_MGIS_EXPORT ElementsJacobians {
    /*!
     * \brief inverses of the jacobians at the points (`nq x dim x dim x
     * ne`, layout of `mfem::GeometricFactors::J`: the component `(i, j)` of
     * the inverse at the point `q` of the element `e` is `invJ[q + nq * (i +
     * dim * (j + dim * e))]`)
     */
    mfem::Vector invJ;
    //! \brief determinants of the jacobians (`nq x ne`)
    mfem::Vector detJ;
  };

  /*!
   * \brief geometric data of all the elements of the mesh for one
   * integration rule: gradients of the shape functions on the reference
   * element (`mfem::DofToQuad`, owned by the finite element), inverses and
   * determinants of the jacobians of the element transformations
   *
   * The data are computed once: they are not updated if the mesh is
   * modified (refinement, displacement of the nodes or of the vertices).
   */
  struct MFEM_MGIS_EXPORT ElementsGeometry {
    //! \brief integration rule
    const mfem::IntegrationRule* rule = nullptr;
    /*!
     * \brief values and gradients of the shape functions of the unknowns at
     * the points, cached by their finite element for the address of the
     * integration rule (see `computeElementsGeometry`)
     */
    const mfem::DofToQuad* maps = nullptr;
    /*!
     * \brief inverses and determinants of the jacobians of the element
     * transformations
     */
    std::shared_ptr<const ElementsJacobians> jacobians;
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
   * The jacobians are computed by a kernel from the coordinates of the
   * nodes of the transformations of the elements, as by
   * `mfem::Mesh::GetElementTransformation`: the nodes of the mesh if it has
   * some, its vertices otherwise (with the finite element of the
   * transformations of MFEM, `mfem::Mesh::GetTransformationFEforElementType`).
   * The gradients of the shape functions of the transformations are
   * computed for the integration rule, without the caches of the finite
   * elements of MFEM, indexed by the address of the rule. The geometric
   * factors of MFEM are not used:
   * - for a mesh without nodes, `mfem::Mesh::GetGeometricFactors` would give
   *   nodes to the mesh (`mfem::Mesh::EnsureNodes`), and the element
   *   transformations of all the computations which follow would use them,
   *   which makes the computations element by element about 50% slower;
   *   computing them from a copy of the vertices in a finite element space
   *   of order one costs 0.7 s on the host and 2.1 s with the GPU enabled
   *   for 1.6 million tetrahedra (construction of the finite element space
   *   and projection of the coordinates on the host), against about 10 ms
   *   for the kernel on the GPU and 0.3 s on the host;
   * - for the quadrangles and hexahedra, they assume a tensor product
   *   integration rule.
   *
   * The same kernel computes the determinants and the inverses of the
   * jacobians, which are the only ones stored (see `ElementsJacobians`).
   *
   * The values and gradients of the shape functions of the unknowns
   * (`maps`) are the ones cached by the finite element of the unknowns: the
   * integration rule must not be destroyed while the finite element space
   * exists (a rule replaced by another one with a different number of
   * points at the same address is detected). NURBS meshes are refused.
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
  /*!
   * \brief computes the geometric data of a list of elements of the mesh
   * for an integration rule, as `computeElementsGeometry`, the data of the
   * `k`-th element of the list being stored as the ones of the element `k`
   * (`ne` is the number of elements of the list)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] g: geometric data
   * \param[in] fes: finite element space
   * \param[in] ir: integration rule
   * \param[in] elements: elements (local numbering of the finite element
   * space)
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool computeElementsGeometry(
      Context&,
      ElementsGeometry&,
      const mfem::FiniteElementSpace&,
      const mfem::IntegrationRule&,
      const mfem::Array<int>&) noexcept;

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_ELEMENTSGEOMETRY_HXX */
