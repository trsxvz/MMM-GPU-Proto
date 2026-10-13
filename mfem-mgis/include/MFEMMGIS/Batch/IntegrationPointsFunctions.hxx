/*!
 * \file   include/MFEMMGIS/Batch/IntegrationPointsFunctions.hxx
 * \brief  This file declares the computations by kernels on functions
 * defined at the integration points of a material, used by the
 * post-processings: integral, right-hand side of the L2 projection, local
 * projections at the nodes of the elements, mass matrices of the elements
 * \date   11/10/2026
 */

#ifndef LIB_MFEM_MGIS_BATCH_INTEGRATIONPOINTSFUNCTIONS_HXX
#define LIB_MFEM_MGIS_BATCH_INTEGRATIONPOINTSFUNCTIONS_HXX

#include <optional>
#include "mfem/linalg/vector.hpp"
#include "mfem/fem/fespace.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"

namespace mfem_mgis {
  // forward declarations
  struct PartialQuadratureSpace;
  struct ImmutablePartialQuadratureFunctionView;
}  // end of namespace mfem_mgis

namespace mfem_mgis::batch {

  /*!
   * \brief gathers the values of a function at the integration points of the
   * elements of its partial quadrature space in a vector of MFEM, which may
   * then be used on the device: the component `c` at the point `q` of the
   * `k`-th element of the space is stored at `c + nc * (q + nq * k)`, `nc`
   * being the number of components of the function and `nq` the number of
   * points per element. The values are read on the host.
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] v: values
   * \param[in] f: function. All the elements of its space must use the same
   * integration rule.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool gatherValues(
      Context& ctx,
      mfem::Vector& v,
      const ImmutablePartialQuadratureFunctionView& f) noexcept;

  /*!
   * \return the integral, over the elements of a partial quadrature space
   * owned by this process, of a scalar function given by its values at the
   * integration points (`w det J` being the weight of a point), or the
   * measure of these elements if no values are given, computed on the device
   * if a device is enabled, or an empty value on failure
   * \param[in, out] ctx: execution context
   * \param[in] qs: partial quadrature space. All its elements must use the
   * integration rule of the geometric data.
   * \param[in] v: values of the function (see `gatherValues`), or null
   * \param[in] g: geometric data of all the elements of the mesh, or of the
   * elements of the space only (`computeElementsGeometry` with the list of
   * the elements of the space), recognized by their number of elements
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<real> computeIntegral(
      Context& ctx,
      const PartialQuadratureSpace& qs,
      const mfem::Vector* v,
      const ElementsGeometry& g) noexcept;

  /*!
   * \brief adds the right-hand side of the L2 projection of a component of a
   * function given at the integration points of the elements of a partial
   * quadrature space, `N f w det J` summed over the points of each element,
   * to a vector of values element by element of a scalar finite element
   * space, in the layout of the element restriction of MFEM (`nd x ne`,
   * native ordering of the degrees of freedom), on the device if a device is
   * enabled
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] r: values element by element
   * \param[in] qs: partial quadrature space. All its elements must use the
   * integration rule of the geometric data.
   * \param[in] v: values of the function (see `gatherValues`)
   * \param[in] nc: number of components of the function
   * \param[in] c: component
   * \param[in] g: geometric data of the finite element space for the
   * integration rule of the space
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addL2ProjectionRightHandSide(
      Context& ctx,
      mfem::Vector& r,
      const PartialQuadratureSpace& qs,
      const mfem::Vector& v,
      const size_type nc,
      const size_type c,
      const ElementsGeometry& g) noexcept;

  /*!
   * \brief maximal number of degrees of freedom per element and per
   * component of the local projections computed by kernels
   */
  inline constexpr int maximumNumberOfNodesOfLocalProjections = 32;

  /*!
   * \brief computes the values at the nodes of the local L2 projections, on
   * each element of a partial quadrature space, of a function given at its
   * integration points, in a vector of values element by element of a
   * finite element space of Lagrange elements (layout of the element
   * restriction of MFEM, `nd x nc x ne`, native ordering of the degrees of
   * freedom), on the device if a device is enabled. The mass matrix of each
   * element is integrated by the integration rule of `m`, the right-hand
   * side by the one of the quadrature space (`g`), and the local system is
   * solved by a Cholesky factorization. The absolute values of the
   * determinants of the jacobians are used, so that the local mass matrices
   * are positive definite for elements of negative orientation. The values
   * of the other elements are not modified.
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] r: values element by element
   * \param[in] qs: partial quadrature space. All its elements must use the
   * integration rule of `g`.
   * \param[in] v: values of the function (see `gatherValues`)
   * \param[in] nc: number of components of the function
   * \param[in] g: geometric data of the finite element space for the
   * integration rule of the quadrature space
   * \param[in] m: geometric data of the finite element space for the
   * integration rule of the mass matrices
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool computeLocalProjections(
      Context& ctx,
      mfem::Vector& r,
      const PartialQuadratureSpace& qs,
      const mfem::Vector& v,
      const size_type nc,
      const ElementsGeometry& g,
      const ElementsGeometry& m) noexcept;

  /*!
   * \brief computes, or adds, the mass matrices of all the elements of a
   * scalar finite element space, `N_i N_j w det J` summed over the points of
   * the integration rule of the geometric data, in the layout of the
   * element matrices of MFEM (`nd x nd x ne`), on the device if a device is
   * enabled
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] m: mass matrices
   * \param[in] g: geometric data of the finite element space
   * \param[in] o: ordering of the degrees of freedom (the lexicographic
   * ordering requires elements with a tensor basis)
   * \param[in] add: if true, the mass matrices are added to `m`
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool computeMassMatrices(
      Context& ctx,
      mfem::Vector& m,
      const ElementsGeometry& g,
      const mfem::ElementDofOrdering o,
      const bool add) noexcept;

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_INTEGRATIONPOINTSFUNCTIONS_HXX */
