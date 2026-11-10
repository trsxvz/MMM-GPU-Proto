/*!
 * \file   include/MFEMMGIS/QuadratureDataLayout.hxx
 * \brief  This file declares the parameter selecting the layouts of the data
 * at the integration points of the materials
 * \date   03/11/2026
 */

#ifndef LIB_MFEM_MGIS_QUADRATUREDATALAYOUT_HXX
#define LIB_MFEM_MGIS_QUADRATUREDATALAYOUT_HXX

#include <optional>
#include "MGIS/FieldView.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Behaviour.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  /*!
   * \brief name of the parameter giving the layouts of the arrays of the
   * data at the integration points of a material (states and tangent
   * operator blocks, see `Material`). Its value is either a string, giving
   * the layout of all the arrays, or a set of parameters giving the layout
   * of each array (`Gradients`, `ThermodynamicForces`,
   * `InternalStateVariables`, `ExternalStateVariables`, `TangentOperator`)
   * and the padding of the arrays (`Padding`, see
   * `mgis::behaviour::MaterialDataManagerLayout`). The layouts are:
   *
   * - `PointByPoint`: the components of a point are contiguous
   *   (`mgis::FieldLayout::POINT_BY_POINT`);
   * - `ComponentByComponent`: the values of a component are contiguous
   *   (`mgis::FieldLayout::COMPONENT_BY_COMPONENT`), which gives coalesced
   *   accesses on a GPU to the kernels using one thread per point;
   * - `Automatic` (default): `ComponentByComponent` if the behaviour of the
   *   material is integrated by batches on a GPU, `PointByPoint` otherwise.
   *
   * It is a parameter of the non linear evolution problems (default layouts
   * of their materials), of their `addBehaviourIntegrator` method (layouts
   * of one material) and of the point-wise models (`PointWiseModel`). The
   * layouts do not change the results.
   */
  inline constexpr const char* quadrature_data_layout_parameter =
      "QuadratureDataLayout";

  /*!
   * \brief layouts of the arrays of a material, as given by the
   * `QuadratureDataLayout` parameter: an empty value selects the layout from
   * the placement of the integration of the behaviour (`Automatic`, the
   * default)
   */
  struct QuadratureDataLayoutSpecification {
    //! \brief layout of the gradients
    std::optional<mgis::FieldLayout> gradients;
    //! \brief layout of the thermodynamic forces
    std::optional<mgis::FieldLayout> thermodynamic_forces;
    //! \brief layout of the internal state variables
    std::optional<mgis::FieldLayout> internal_state_variables;
    //! \brief layout of the external state variables with several components
    std::optional<mgis::FieldLayout> external_state_variables;
    //! \brief layout of the tangent operator blocks
    std::optional<mgis::FieldLayout> tangent_operator;
    //! \brief padding of the arrays
    size_type padding = 0;
  };

  /*!
   * \return the specification of the layouts given by the environment
   * variables `MFEM_MGIS_QUADRATURE_DATA_LAYOUT` (`Automatic` by default,
   * `PointByPoint` or `ComponentByComponent`) and
   * `MFEM_MGIS_QUADRATURE_DATA_PADDING` (a non negative integer, zero by
   * default), for instance to run a whole test suite or an application with
   * the component by component layout, or with padded arrays (whose unused
   * values are NaN: a code reading the arrays without their views fails)
   * \param[in, out] ctx: execution context
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]] std::optional<QuadratureDataLayoutSpecification>
  getDefaultQuadratureDataLayoutSpecification(Context& ctx) noexcept;
  /*!
   * \return the specification of the layouts given by the
   * `QuadratureDataLayout` parameter, the layouts and the padding which are
   * not given being the ones of the given default specification
   * \param[in, out] ctx: execution context
   * \param[in] params: parameters
   * \param[in] d: default specification
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]] std::optional<QuadratureDataLayoutSpecification>
  getQuadratureDataLayoutSpecification(
      Context& ctx,
      const Parameters& params,
      const QuadratureDataLayoutSpecification& d) noexcept;
  /*!
   * \return the specification of the layouts given by the
   * `QuadratureDataLayout` parameter (the layouts and the padding which are
   * not given being the default ones, automatic without padding), or, if
   * this parameter is not given, by the environment variables (see
   * `getDefaultQuadratureDataLayoutSpecification`)
   * \param[in] params: parameters
   * \throws std::runtime_error if the value of the parameter or of the
   * environment variables is invalid
   */
  MFEM_MGIS_EXPORT [[nodiscard]] QuadratureDataLayoutSpecification
  getQuadratureDataLayoutSpecification(attributes::MayThrow,
                                       const Parameters& params);
  /*!
   * \return the layouts of the arrays of a material given by a
   * specification, the automatic layouts being chosen from the placement of
   * the integration of its behaviour: `ComponentByComponent` if it is
   * integrated by batches on a GPU, `PointByPoint` otherwise
   * \param[in] s: specification
   * \param[in] b: behaviour
   */
  MFEM_MGIS_EXPORT [[nodiscard]] mgis::behaviour::MaterialDataManagerLayout
  resolveQuadratureDataLayout(const QuadratureDataLayoutSpecification& s,
                              const Behaviour& b) noexcept;
  /*!
   * \return the layouts of the arrays of a material given by the
   * `QuadratureDataLayout` parameter (see
   * `resolveQuadratureDataLayout`), the arrays being stored point by point
   * without padding if this parameter is not given. This function is used
   * by the generators of the behaviour integrators, which receive the
   * layouts resolved by the non linear evolution problem, given only if
   * they are not point by point without padding.
   * \param[in, out] ctx: execution context
   * \param[in] params: parameters
   * \param[in] b: behaviour
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]] std::optional<mgis::behaviour::MaterialDataManagerLayout>
  getQuadratureDataLayout(Context& ctx,
                          const Parameters& params,
                          const Behaviour& b) noexcept;
  /*!
   * \return the value of the `QuadratureDataLayout` parameter selecting the
   * given layouts: a string if all the arrays have the same layout without
   * padding, a set of parameters otherwise
   * \param[in] l: layouts
   */
  MFEM_MGIS_EXPORT [[nodiscard]] Parameter toParameter(
      const mgis::behaviour::MaterialDataManagerLayout& l);
  /*!
   * \return true if all the arrays of a material are stored point by point,
   * with or without padding: the values of a point of each array are then
   * contiguous, as required by the integration point by point of the
   * behaviour
   * \param[in] l: layouts
   */
  [[nodiscard]] inline bool isStoredPointByPoint(
      const mgis::behaviour::MaterialDataManagerLayout& l) noexcept {
    const auto ol = mgis::behaviour::getUniformLayout(l);
    return ol.has_value() && (*ol == mgis::FieldLayout::POINT_BY_POINT);
  }  // end of isStoredPointByPoint
  /*!
   * \return the value of the `QuadratureDataLayout` parameter selecting a
   * layout
   * \param[in] l: layout
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const char* toString(
      const mgis::FieldLayout l) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_QUADRATUREDATALAYOUT_HXX */
