/*!
 * \file   src/QuadratureDataLayout.cxx
 * \brief  This file implements the functions related to the parameter
 * selecting the layout of the data at the integration points of the
 * materials
 * \date   03/11/2026
 */

#include <cstdlib>
#include <string>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/QuadratureDataLayout.hxx"

namespace mfem_mgis {

  /*!
   * \return the layout selected by a value of the `QuadratureDataLayout`
   * parameter or of the `MFEM_MGIS_QUADRATURE_DATA_LAYOUT` environment
   * variable
   * \param[in, out] ctx: execution context
   * \param[in] v: value
   * \param[in] origin: origin of the value, for the error message
   */
  [[nodiscard]] static std::optional<mgis::FieldLayout>
  convertToQuadratureDataLayout(Context& ctx,
                                const std::string& v,
                                const std::string& origin) noexcept {
    if (v == "PointByPoint") {
      return mgis::FieldLayout::POINT_BY_POINT;
    }
    if (v == "ComponentByComponent") {
      return mgis::FieldLayout::COMPONENT_BY_COMPONENT;
    }
    return ctx.registerErrorMessage(
        "invalid value '" + v + "' for " + origin +
        " (expected 'PointByPoint' or 'ComponentByComponent')");
  }  // end of convertToQuadratureDataLayout

  std::optional<mgis::FieldLayout> getDefaultQuadratureDataLayout(
      Context& ctx) noexcept {
    const auto* const v = "MFEM_MGIS_QUADRATURE_DATA_LAYOUT";
    const auto* const e = std::getenv(v);
    if (e == nullptr) {
      return mgis::FieldLayout::POINT_BY_POINT;
    }
    return convertToQuadratureDataLayout(
        ctx, e, "environment variable '" + std::string{v} + "'");
  }  // end of getDefaultQuadratureDataLayout

  std::optional<mgis::FieldLayout> getQuadratureDataLayout(
      Context& ctx,
      const Parameters& params,
      const mgis::FieldLayout d) noexcept {
    const auto* const n = quadrature_data_layout_parameter;
    if (!contains(params, n)) {
      return d;
    }
    const auto ov = get<std::string>(ctx, params, n);
    if (isInvalid(ov)) {
      return {};
    }
    return convertToQuadratureDataLayout(ctx, *ov,
                                         "parameter '" + std::string{n} + "'");
  }  // end of getQuadratureDataLayout

  mgis::FieldLayout getQuadratureDataLayout(attributes::MayThrow,
                                            const Parameters& params) {
    // the environment variable is only read if the parameter is not given
    auto ctx = Context{};
    const auto ol = contains(params, quadrature_data_layout_parameter)
                        ? getQuadratureDataLayout(
                              ctx, params, mgis::FieldLayout::POINT_BY_POINT)
                        : getDefaultQuadratureDataLayout(ctx);
    if (isInvalid(ol)) {
      raise(ctx.getErrorMessage());
    }
    return *ol;
  }  // end of getQuadratureDataLayout

  const char* toString(const mgis::FieldLayout l) noexcept {
    return l == mgis::FieldLayout::POINT_BY_POINT ? "PointByPoint"
                                                  : "ComponentByComponent";
  }  // end of toString

}  // end of namespace mfem_mgis
