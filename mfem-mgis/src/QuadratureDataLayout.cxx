/*!
 * \file   src/QuadratureDataLayout.cxx
 * \brief  This file implements the functions related to the parameter
 * selecting the layouts of the data at the integration points of the
 * materials
 * \date   03/11/2026
 */

#include <array>
#include <cstdlib>
#include <string>
#include <tuple>
#include <utility>
#include "mfem/general/device.hpp"
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/QuadratureDataLayout.hxx"

namespace mfem_mgis {

  //! \brief names of the arrays in the `QuadratureDataLayout` parameter
  static constexpr auto quadrature_data_arrays = std::array<const char*, 5>{
      "Gradients", "ThermodynamicForces", "InternalStateVariables",
      "ExternalStateVariables", "TangentOperator"};

  //! \brief name of the padding in the `QuadratureDataLayout` parameter
  static constexpr const char* quadrature_data_padding = "Padding";

  /*!
   * \return the members of a specification, in the order of
   * `quadrature_data_arrays`
   * \param[in] s: specification
   */
  static std::array<std::optional<mgis::FieldLayout>*, 5> getArraysLayouts(
      QuadratureDataLayoutSpecification& s) noexcept {
    return {&s.gradients, &s.thermodynamic_forces, &s.internal_state_variables,
            &s.external_state_variables, &s.tangent_operator};
  }  // end of getArraysLayouts

  /*!
   * \return the layout selected by a value of the `QuadratureDataLayout`
   * parameter or of the `MFEM_MGIS_QUADRATURE_DATA_LAYOUT` environment
   * variable, an empty value meaning `Automatic`
   * \param[in, out] ctx: execution context
   * \param[in] v: value
   * \param[in] origin: origin of the value, for the error message
   */
  [[nodiscard]] static std::optional<std::optional<mgis::FieldLayout>>
  convertToQuadratureDataLayout(Context& ctx,
                                const std::string& v,
                                const std::string& origin) noexcept {
    using Result = std::optional<std::optional<mgis::FieldLayout>>;
    if (v == "PointByPoint") {
      return Result{mgis::FieldLayout::POINT_BY_POINT};
    }
    if (v == "ComponentByComponent") {
      return Result{mgis::FieldLayout::COMPONENT_BY_COMPONENT};
    }
    if (v == "Automatic") {
      return Result{std::optional<mgis::FieldLayout>{}};
    }
    std::ignore = ctx.registerErrorMessage(
        "invalid value '" + v + "' for " + origin +
        " (expected 'PointByPoint', 'ComponentByComponent' or 'Automatic')");
    return Result{};
  }  // end of convertToQuadratureDataLayout

  std::optional<QuadratureDataLayoutSpecification>
  getDefaultQuadratureDataLayoutSpecification(Context& ctx) noexcept {
    auto s = QuadratureDataLayoutSpecification{};
    const auto* const vl = "MFEM_MGIS_QUADRATURE_DATA_LAYOUT";
    if (const auto* const e = std::getenv(vl); e != nullptr) {
      const auto ol = convertToQuadratureDataLayout(
          ctx, e, "environment variable '" + std::string{vl} + "'");
      if (isInvalid(ol)) {
        return {};
      }
      for (auto* const l : getArraysLayouts(s)) {
        *l = *ol;
      }
    }
    const auto* const vp = "MFEM_MGIS_QUADRATURE_DATA_PADDING";
    if (const auto* const e = std::getenv(vp); e != nullptr) {
      try {
        auto pos = std::size_t{};
        const auto p = std::stoi(e, &pos);
        if ((pos != std::string_view{e}.size()) || (p < 0)) {
          throw std::invalid_argument("invalid padding");
        }
        s.padding = static_cast<size_type>(p);
      } catch (...) {
        return ctx.registerErrorMessage("invalid value '" + std::string{e} +
                                        "' for environment variable '" +
                                        std::string{vp} +
                                        "' (expected a non negative integer)");
      }
    }
    return s;
  }  // end of getDefaultQuadratureDataLayoutSpecification

  std::optional<QuadratureDataLayoutSpecification>
  getQuadratureDataLayoutSpecification(
      Context& ctx,
      const Parameters& params,
      const QuadratureDataLayoutSpecification& d) noexcept {
    const auto* const n = quadrature_data_layout_parameter;
    if (!contains(params, n)) {
      return d;
    }
    const auto origin = "parameter '" + std::string{n} + "'";
    auto s = d;
    const auto oparameter = get(ctx, params, n);
    if (isInvalid(oparameter)) {
      return {};
    }
    if (is<std::string>(*oparameter)) {
      const auto ov = get<std::string>(ctx, *oparameter);
      if (isInvalid(ov)) {
        return {};
      }
      const auto ol = convertToQuadratureDataLayout(ctx, *ov, origin);
      if (isInvalid(ol)) {
        return {};
      }
      for (auto* const l : getArraysLayouts(s)) {
        *l = *ol;
      }
      return s;
    }
    const auto op = get<Parameters>(ctx, *oparameter);
    if (isInvalid(op)) {
      return {};
    }
    auto names = std::vector<std::string>(quadrature_data_arrays.begin(),
                                          quadrature_data_arrays.end());
    names.push_back(quadrature_data_padding);
    if (!checkParameters(ctx, *op, names)) {
      return {};
    }
    const auto layouts = getArraysLayouts(s);
    for (std::size_t i = 0; i != quadrature_data_arrays.size(); ++i) {
      const auto* const a = quadrature_data_arrays[i];
      if (!contains(*op, a)) {
        continue;
      }
      const auto ov = get<std::string>(ctx, *op, a);
      if (isInvalid(ov)) {
        return {};
      }
      const auto ol = convertToQuadratureDataLayout(
          ctx, *ov, origin + " (array '" + std::string{a} + "')");
      if (isInvalid(ol)) {
        return {};
      }
      *(layouts[i]) = *ol;
    }
    if (contains(*op, quadrature_data_padding)) {
      const auto op2 = get<int>(ctx, *op, quadrature_data_padding);
      if (isInvalid(op2)) {
        return {};
      }
      if (*op2 < 0) {
        return ctx.registerErrorMessage("invalid padding for " + origin +
                                        " (expected a non negative integer)");
      }
      s.padding = static_cast<size_type>(*op2);
    }
    return s;
  }  // end of getQuadratureDataLayoutSpecification

  //! \return the specification of the point by point layout without padding
  static QuadratureDataLayoutSpecification
  getPointByPointLayoutSpecification() noexcept {
    constexpr auto pbp = mgis::FieldLayout::POINT_BY_POINT;
    return {.gradients = pbp,
            .thermodynamic_forces = pbp,
            .internal_state_variables = pbp,
            .external_state_variables = pbp,
            .tangent_operator = pbp,
            .padding = 0};
  }  // end of getPointByPointLayoutSpecification

  QuadratureDataLayoutSpecification getQuadratureDataLayoutSpecification(
      attributes::MayThrow, const Parameters& params) {
    // the environment variables are only read if the parameter is not given
    auto ctx = Context{};
    const auto os = contains(params, quadrature_data_layout_parameter)
                        ? getQuadratureDataLayoutSpecification(
                              ctx, params, QuadratureDataLayoutSpecification{})
                        : getDefaultQuadratureDataLayoutSpecification(ctx);
    if (isInvalid(os)) {
      raise(ctx.getErrorMessage());
    }
    return *os;
  }  // end of getQuadratureDataLayoutSpecification

  mgis::behaviour::MaterialDataManagerLayout resolveQuadratureDataLayout(
      const QuadratureDataLayoutSpecification& s, const Behaviour& b) noexcept {
    // the behaviours integrated on a device which is not a GPU (debug device
    // of MFEM) are integrated on the CPU (see
    // `checkBehaviourExecutionLocation`)
    const auto on_gpu = isIntegratedByBatchesOnDevice(b) &&
                        mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                                             mfem::Backend::HIP_MASK);
    const auto automatic = on_gpu ? mgis::FieldLayout::COMPONENT_BY_COMPONENT
                                  : mgis::FieldLayout::POINT_BY_POINT;
    auto l = mgis::behaviour::MaterialDataManagerLayout{};
    l.gradients = s.gradients.value_or(automatic);
    l.thermodynamic_forces = s.thermodynamic_forces.value_or(automatic);
    l.internal_state_variables = s.internal_state_variables.value_or(automatic);
    l.external_state_variables = s.external_state_variables.value_or(automatic);
    l.tangent_operator = s.tangent_operator.value_or(automatic);
    l.padding = s.padding;
    return l;
  }  // end of resolveQuadratureDataLayout

  std::optional<mgis::behaviour::MaterialDataManagerLayout>
  getQuadratureDataLayout(Context& ctx,
                          const Parameters& params,
                          const Behaviour& b) noexcept {
    const auto os = getQuadratureDataLayoutSpecification(
        ctx, params, getPointByPointLayoutSpecification());
    if (isInvalid(os)) {
      return {};
    }
    return resolveQuadratureDataLayout(*os, b);
  }  // end of getQuadratureDataLayout

  Parameter toParameter(const mgis::behaviour::MaterialDataManagerLayout& l) {
    if (l.padding == 0) {
      if (const auto ol = mgis::behaviour::getUniformLayout(l);
          ol.has_value()) {
        return std::string{toString(*ol)};
      }
    }
    return Parameters{
        {quadrature_data_arrays[0], std::string{toString(l.gradients)}},
        {quadrature_data_arrays[1],
         std::string{toString(l.thermodynamic_forces)}},
        {quadrature_data_arrays[2],
         std::string{toString(l.internal_state_variables)}},
        {quadrature_data_arrays[3],
         std::string{toString(l.external_state_variables)}},
        {quadrature_data_arrays[4], std::string{toString(l.tangent_operator)}},
        {quadrature_data_padding, static_cast<int>(l.padding)}};
  }  // end of toParameter

  const char* toString(const mgis::FieldLayout l) noexcept {
    return l == mgis::FieldLayout::POINT_BY_POINT ? "PointByPoint"
                                                  : "ComponentByComponent";
  }  // end of toString

}  // end of namespace mfem_mgis
