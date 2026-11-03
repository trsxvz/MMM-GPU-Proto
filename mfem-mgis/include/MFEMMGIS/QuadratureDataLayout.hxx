/*!
 * \file   include/MFEMMGIS/QuadratureDataLayout.hxx
 * \brief  This file declares the parameter selecting the layout of the data
 * at the integration points of the materials
 * \date   03/11/2026
 */

#ifndef LIB_MFEM_MGIS_QUADRATUREDATALAYOUT_HXX
#define LIB_MFEM_MGIS_QUADRATUREDATALAYOUT_HXX

#include <optional>
#include "MGIS/FieldView.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  /*!
   * \brief name of the parameter giving the layout of the arrays of the
   * data at the integration points of a material (states and tangent
   * operator blocks, see `Material`):
   *
   * - `PointByPoint` (default): the components of a point are contiguous
   *   (`mgis::FieldLayout::POINT_BY_POINT`);
   * - `ComponentByComponent`: the values of a component are contiguous
   *   (`mgis::FieldLayout::COMPONENT_BY_COMPONENT`), which gives coalesced
   *   accesses on a GPU to the kernels using one thread per point.
   *
   * It is a parameter of the non linear evolution problems (default layout
   * of their materials), of their `addBehaviourIntegrator` method (layout
   * of one material) and of the point-wise models (`PointWiseModel`). The
   * layout does not change the results.
   */
  inline constexpr const char* quadrature_data_layout_parameter =
      "QuadratureDataLayout";

  /*!
   * \return the layout given by the environment variable
   * `MFEM_MGIS_QUADRATURE_DATA_LAYOUT` (`PointByPoint` or
   * `ComponentByComponent`), for instance to run a whole test suite or an
   * application with the component by component layout, or the point by
   * point layout if this variable is not defined
   * \param[in, out] ctx: execution context
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<mgis::FieldLayout>
  getDefaultQuadratureDataLayout(Context&) noexcept;
  /*!
   * \return the layout given by the `QuadratureDataLayout` parameter, or
   * the given default layout if this parameter is not given
   * \param[in, out] ctx: execution context
   * \param[in] params: parameters
   * \param[in] d: default layout
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<mgis::FieldLayout>
  getQuadratureDataLayout(Context&,
                          const Parameters&,
                          const mgis::FieldLayout) noexcept;
  /*!
   * \return the layout given by the `QuadratureDataLayout` parameter, or
   * else by the `MFEM_MGIS_QUADRATURE_DATA_LAYOUT` environment variable
   * (see `getDefaultQuadratureDataLayout`)
   * \param[in] params: parameters
   * \throws std::runtime_error if the value of the parameter or of the
   * environment variable is invalid
   */
  MFEM_MGIS_EXPORT [[nodiscard]] mgis::FieldLayout getQuadratureDataLayout(
      attributes::MayThrow, const Parameters&);
  /*!
   * \return the value of the `QuadratureDataLayout` parameter selecting a
   * layout
   * \param[in] l: layout
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const char* toString(
      const mgis::FieldLayout) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_QUADRATUREDATALAYOUT_HXX */
