/*!
 * \file   src/Batch/Gradients.cxx
 * \brief  This file implements the computation of the gradients of the
 * behaviours by kernels (host part: checks and arguments of the kernels)
 * \date   09/10/2026
 */

#include <limits>
#include <string>
#include <cstdint>
#include <exception>
#include "mfem/general/device.hpp"
#include "mfem/fem/fespace.hpp"
#include "mfem/fem/restriction.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES
#include "GradientsKernels.hxx"
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

namespace mfem_mgis::batch {

  /*!
   * \return true if the behaviour is a stationary non linear heat transfer
   * behaviour, recognised by its structure, whatever the names of its
   * variables: one gradient (the temperature gradient) and one
   * thermodynamic force (the heat flux), both vectors, derivatives of the
   * heat flux with respect to the temperature gradient and to the
   * temperature. The heat transfer behaviour integrators only check the
   * sizes of the variables: the behaviours that they accept and that are
   * not recognised here are refused by the kernels.
   * \param[in] b: behaviour
   */
  [[nodiscard]] static bool isHeatTransferBehaviour(
      const Behaviour& b) noexcept {
    using Variable = mgis::behaviour::Variable;
    if ((b.btype != Behaviour::GENERALBEHAVIOUR) ||
        (b.gradients.size() != 1u) || (b.thermodynamic_forces.size() != 1u) ||
        (b.to_blocks.size() != 2u)) {
      return false;
    }
    const auto& g = b.gradients[0];
    const auto& f = b.thermodynamic_forces[0];
    return (g.type == Variable::VECTOR) && (f.type == Variable::VECTOR) &&
           (b.to_blocks[0].first.name == f.name) &&
           (b.to_blocks[0].second.name == g.name) &&
           (b.to_blocks[1].first.name == f.name) &&
           (b.to_blocks[1].second.name == "Temperature");
  }  // end of isHeatTransferBehaviour

  /*!
   * \return true if the behaviour describes the transient term of the heat
   * equation, recognised by its structure: one scalar gradient (the
   * temperature at the integration point) and one scalar thermodynamic force
   * (the opposite of the rate of the enthalpy), and the derivative of the
   * force with respect to the gradient (the behaviour integrator
   * `TransientHeatTransfer` only checks the sizes and the symmetry)
   * \param[in] b: behaviour
   */
  [[nodiscard]] static bool isHeatStorageBehaviour(
      const Behaviour& b) noexcept {
    using Variable = mgis::behaviour::Variable;
    if ((b.btype != Behaviour::GENERALBEHAVIOUR) ||
        (b.gradients.size() != 1u) || (b.thermodynamic_forces.size() != 1u) ||
        (b.to_blocks.size() != 1u)) {
      return false;
    }
    const auto& g = b.gradients[0];
    const auto& f = b.thermodynamic_forces[0];
    return (g.type == Variable::SCALAR) && (f.type == Variable::SCALAR) &&
           (b.to_blocks[0].first.name == f.name) &&
           (b.to_blocks[0].second.name == g.name);
  }  // end of isHeatStorageBehaviour

  /*!
   * \return true if the behaviour is a micromorphic damage behaviour,
   * recognised by its structure (the behaviour integrators
   * `MicromorphicDamage` only check the sizes and the symmetry): a scalar
   * gradient (the micromorphic damage) and a
   * vector gradient (its gradient), their scalar and vector dual forces, and
   * the derivatives of each force with respect to the associated gradient
   * \param[in] b: behaviour
   */
  [[nodiscard]] static bool isMicromorphicDamageBehaviour(
      const Behaviour& b) noexcept {
    using Variable = mgis::behaviour::Variable;
    if ((b.btype != Behaviour::GENERALBEHAVIOUR) ||
        (b.gradients.size() != 2u) || (b.thermodynamic_forces.size() != 2u) ||
        (b.to_blocks.size() != 2u)) {
      return false;
    }
    const auto& g = b.gradients;
    const auto& f = b.thermodynamic_forces;
    return (g[0].type == Variable::SCALAR) && (g[1].type == Variable::VECTOR) &&
           (f[0].type == Variable::SCALAR) && (f[1].type == Variable::VECTOR) &&
           (b.to_blocks[0].first.name == f[0].name) &&
           (b.to_blocks[0].second.name == g[0].name) &&
           (b.to_blocks[1].first.name == f[1].name) &&
           (b.to_blocks[1].second.name == g[1].name);
  }  // end of isMicromorphicDamageBehaviour

  std::optional<GradientsKind> getGradientsKind(const Behaviour& b) noexcept {
    const auto h = b.hypothesis;
    const auto tridimensional = h == Hypothesis::TRIDIMENSIONAL;
    const auto bidimensional =
        (h == Hypothesis::PLANESTRAIN) || (h == Hypothesis::PLANESTRESS);
    if ((b.btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) &&
        (b.kinematic == Behaviour::SMALLSTRAINKINEMATIC)) {
      if (tridimensional) {
        return GradientsKind::SMALL_STRAIN_3D;
      }
      if (bidimensional) {
        return GradientsKind::SMALL_STRAIN_2D;
      }
    }
    if ((b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) &&
        (b.kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY)) {
      if (tridimensional) {
        return GradientsKind::DEFORMATION_GRADIENT_3D;
      }
      if (bidimensional) {
        return GradientsKind::DEFORMATION_GRADIENT_2D;
      }
    }
    if (isHeatTransferBehaviour(b)) {
      if (tridimensional) {
        return GradientsKind::TEMPERATURE_GRADIENT_3D;
      }
      if (bidimensional) {
        return GradientsKind::TEMPERATURE_GRADIENT_2D;
      }
    }
    if (isHeatStorageBehaviour(b)) {
      if (tridimensional) {
        return GradientsKind::TEMPERATURE_3D;
      }
      if (bidimensional) {
        return GradientsKind::TEMPERATURE_2D;
      }
    }
    if (isMicromorphicDamageBehaviour(b)) {
      if (tridimensional) {
        return GradientsKind::MICROMORPHIC_DAMAGE_3D;
      }
      if (bidimensional) {
        return GradientsKind::MICROMORPHIC_DAMAGE_2D;
      }
    }
    return {};
  }  // end of getGradientsKind

  bool areGradientsKernelsAvailable() noexcept {
#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES
    return true;
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    return false;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of areGradientsKernelsAvailable

  /*!
   * \brief resizes an array of the memory of MFEM, used on the device if a
   * device is enabled
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] v: array
   * \param[in] s: size
   */
  static bool resizeArray(Context& ctx, mfem::Vector& v, const std::int64_t s) {
    if (s > std::numeric_limits<int>::max()) {
      return ctx.registerErrorMessage("array too large for MFEM");
    }
    v.UseDevice(true);
    v.SetSize(static_cast<int>(s));
    return true;
  }  // end of resizeArray

  bool computeUnknownsDerivatives(Context& ctx,
                                  UnknownsDerivatives& d,
                                  const mfem::FiniteElementSpace& fes,
                                  const mfem::Vector& u,
                                  const mfem::IntegrationRule& ir,
                                  const bool values) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(d);
    static_cast<void>(fes);
    static_cast<void>(u);
    static_cast<void>(ir);
    static_cast<void>(values);
    return ctx.registerErrorMessage(
        "computeUnknownsDerivatives: the kernels computing the gradients are "
        "not compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      const auto& mesh = *(fes.GetMesh());
      if (mesh.GetNumGeometries(mesh.Dimension()) > 1) {
        return ctx.registerErrorMessage(
            "computeUnknownsDerivatives: meshes mixing several geometries "
            "are not supported");
      }
      d.rule = &ir;
      d.vdim = fes.GetVDim();
      d.dim = mesh.Dimension();
      d.nq = ir.GetNPoints();
      d.ne = fes.GetNE();
      const auto size = std::int64_t{d.vdim} * std::int64_t{d.dim} *
                        std::int64_t{d.nq} * std::int64_t{d.ne};
      if (!resizeArray(ctx, d.values, size)) {
        return false;
      }
      d.has_point_values = values;
      if ((values) && (!resizeArray(ctx, d.point_values,
                                    std::int64_t{d.vdim} * std::int64_t{d.nq} *
                                        std::int64_t{d.ne}))) {
        return false;
      }
      if (d.ne == 0) {
        return true;
      }
      const auto* const r =
          fes.GetElementRestriction(mfem::ElementDofOrdering::NATIVE);
      if (!resizeArray(ctx, d.element_values, r->Height())) {
        return false;
      }
      r->Mult(u, d.element_values);
      const auto* const fe = fes.GetTypicalFE();
      const auto& maps = fe->GetDofToQuad(ir, mfem::DofToQuad::FULL);
      // the values and gradients of the shape functions are cached by the
      // finite element for the address of the integration rule: a different
      // number of points reveals a rule destroyed and replaced by another at
      // the same address
      // (the number of degrees of freedom is the one of the element, which
      // creates the cache)
      if (maps.nqpt != d.nq) {
        return ctx.registerErrorMessage(
            "computeUnknownsDerivatives: the values of the shape functions "
            "cached by the finite element do not match the integration rule");
      }
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::UnknownsDerivativesKernelArguments{};
      a.number_of_elements = d.ne;
      a.nd = fe->GetDof();
      a.nq = d.nq;
      a.vdim = d.vdim;
      a.dim = d.dim;
      if (std::int64_t{a.nd} * a.vdim * a.number_of_elements !=
          std::int64_t{d.element_values.Size()}) {
        return ctx.registerErrorMessage(
            "computeUnknownsDerivatives: unexpected number of values element "
            "by element");
      }
      if (std::int64_t{a.number_of_elements} * a.nq >
          std::numeric_limits<int>::max()) {
        return ctx.registerErrorMessage(
            "computeUnknownsDerivatives: too many points for the kernel");
      }
      a.element_values = d.element_values.Read(use_device);
      a.shape_functions_gradients = maps.G.Read(use_device);
      a.derivatives = d.values.Write(use_device);
      if (values) {
        a.shape_functions = maps.B.Read(use_device);
        a.point_values = d.point_values.Write(use_device);
      }
      internals::launchUnknownsDerivativesKernel(a, use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"computeUnknownsDerivatives: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "computeUnknownsDerivatives: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of computeUnknownsDerivatives

  size_type getGradientsSize(const GradientsKind k) noexcept {
    switch (k) {
      case GradientsKind::SMALL_STRAIN_3D:
        return 6;
      case GradientsKind::SMALL_STRAIN_2D:
        return 4;
      case GradientsKind::DEFORMATION_GRADIENT_3D:
        return 9;
      case GradientsKind::DEFORMATION_GRADIENT_2D:
        return 5;
      case GradientsKind::TEMPERATURE_GRADIENT_3D:
        return 3;
      case GradientsKind::TEMPERATURE_GRADIENT_2D:
        return 2;
      case GradientsKind::MICROMORPHIC_DAMAGE_3D:
        return 4;
      case GradientsKind::MICROMORPHIC_DAMAGE_2D:
        return 3;
      case GradientsKind::TEMPERATURE_3D:
      case GradientsKind::TEMPERATURE_2D:
        break;
    }
    return 1;
  }  // end of getGradientsSize

  bool gradientsRequireUnknownsValues(const GradientsKind k) noexcept {
    return (k == GradientsKind::TEMPERATURE_3D) ||
           (k == GradientsKind::TEMPERATURE_2D) ||
           (k == GradientsKind::MICROMORPHIC_DAMAGE_3D) ||
           (k == GradientsKind::MICROMORPHIC_DAMAGE_2D);
  }  // end of gradientsRequireUnknownsValues

  bool gradientsRequireUnknownsDerivatives(const GradientsKind k) noexcept {
    return (k != GradientsKind::TEMPERATURE_3D) &&
           (k != GradientsKind::TEMPERATURE_2D);
  }  // end of gradientsRequireUnknownsDerivatives

  bool areUnknownsValuesExternalStateVariables(const GradientsKind k) noexcept {
    return (k == GradientsKind::TEMPERATURE_GRADIENT_3D) ||
           (k == GradientsKind::TEMPERATURE_GRADIENT_2D);
  }  // end of areUnknownsValuesExternalStateVariables

  bool requiresUnknownsValues(const GradientsKind k) noexcept {
    return gradientsRequireUnknownsValues(k) ||
           areUnknownsValuesExternalStateVariables(k);
  }  // end of requiresUnknownsValues

  size_type getTangentOperatorBlocksSize(const GradientsKind k) noexcept {
    if ((k == GradientsKind::MICROMORPHIC_DAMAGE_3D) ||
        (k == GradientsKind::MICROMORPHIC_DAMAGE_2D)) {
      // block diagonal tangent operator
      const auto d = getSpaceDimension(k);
      return 1 + d * d;
    }
    const auto n = getGradientsSize(k);
    // derivatives of the heat flux with respect to the temperature
    return areUnknownsValuesExternalStateVariables(k) ? n * n + n : n * n;
  }  // end of getTangentOperatorBlocksSize

  int getSpaceDimension(const GradientsKind k) noexcept {
    return ((k == GradientsKind::SMALL_STRAIN_3D) ||
            (k == GradientsKind::DEFORMATION_GRADIENT_3D) ||
            (k == GradientsKind::TEMPERATURE_GRADIENT_3D) ||
            (k == GradientsKind::TEMPERATURE_3D) ||
            (k == GradientsKind::MICROMORPHIC_DAMAGE_3D))
               ? 3
               : 2;
  }  // end of getSpaceDimension

  int getUnknownsSize(const GradientsKind k) noexcept {
    return requiresUnknownsValues(k) ? 1 : getSpaceDimension(k);
  }  // end of getUnknownsSize

#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES

  /*!
   * \brief checks the arguments of the computation of the gradients of a
   * material and fills the arguments of the kernel which do not depend on
   * the destination
   * \return true on success, false if there is nothing to do or on error
   * (see `error`)
   * \param[in, out] ctx: execution context
   * \param[out] a: arguments of the kernel
   * \param[out] error: true on error
   * \param[in] name: name of the calling function, for the messages
   * \param[in] qs: quadrature space of the material
   * \param[in] d: derivatives of the unknowns
   * \param[in] g: geometric data
   * \param[in] k: kind of the gradients
   * \param[in] use_device: if true, the pointers are located on the device
   */
  static bool getGradientsKernelArguments(
      Context& ctx,
      internals::GradientsKernelArguments& a,
      bool& error,
      const std::string& name,
      const PartialQuadratureSpace& qs,
      const UnknownsDerivatives& d,
      const ElementsGeometry& g,
      const GradientsKind k,
      const bool use_device) {
    error = false;
    const auto& elements = qs.getElements();
    const auto& offsets = qs.getElementOffsets();
    if (elements.Size() == 0) {
      return false;
    }
    error = true;
    if ((d.rule == nullptr) || (qs.getUniformIntegrationRule() != d.rule) ||
        (g.rule != d.rule) || (g.jacobians == nullptr) || (g.nq != d.nq) ||
        (g.ne != d.ne)) {
      return ctx.registerErrorMessage(
          name +
          ": the elements of the material do not use the integration rule "
          "of the derivatives of the unknowns and of the geometric data");
    }
    const auto dim = getSpaceDimension(k);
    const auto vdim = getUnknownsSize(k);
    if ((d.vdim != vdim) || (d.dim != dim)) {
      return ctx.registerErrorMessage(
          name +
          ": unexpected number of components of the unknowns or dimension "
          "of the space (" +
          std::to_string(d.vdim) + ", " + std::to_string(d.dim) +
          "), the expected values being " + std::to_string(vdim) + " and " +
          std::to_string(dim));
    }
    const auto n = std::int64_t{elements.Size()} * d.nq;
    if (n != std::int64_t{qs.getNumberOfIntegrationPoints()}) {
      return ctx.registerErrorMessage(
          name + ": unexpected number of integration points");
    }
    a.number_of_elements = elements.Size();
    a.nq = d.nq;
    a.vdim = d.vdim;
    a.dim = d.dim;
    a.number_of_points = static_cast<int>(n);
    a.gradients_size = static_cast<int>(getGradientsSize(k));
    a.elements = elements.Read(use_device);
    a.offsets = offsets.Read(use_device);
    a.derivatives = d.values.Read(use_device);
    if (gradientsRequireUnknownsValues(k)) {
      if (!d.has_point_values) {
        return ctx.registerErrorMessage(
            name +
            ": the values of the unknowns at the points are not computed");
      }
      a.unknowns_values = d.point_values.Read(use_device);
    }
    a.inverse_jacobians = g.jacobians->invJ.Read(use_device);
    error = false;
    return true;
  }  // end of getGradientsKernelArguments

#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

  bool computeGradients(Context& ctx,
                        Material& m,
                        const UnknownsDerivatives& d,
                        const ElementsGeometry& g,
                        const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(m);
    static_cast<void>(d);
    static_cast<void>(g);
    static_cast<void>(k);
    return ctx.registerErrorMessage(
        "computeGradients: the kernels computing the gradients are not "
        "compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    const auto gsize = getGradientsSize(k);
    if (m.s1.gradients_stride != gsize) {
      return ctx.registerErrorMessage(
          "computeGradients: unexpected number of components of the "
          "gradients (" +
          std::to_string(m.s1.gradients_stride) + " instead of " +
          std::to_string(gsize) + ")");
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::GradientsKernelArguments{};
      auto error = false;
      if (!getGradientsKernelArguments(ctx, a, error, "computeGradients",
                                       m.getPartialQuadratureSpace(), d, g, k,
                                       use_device)) {
        return !error;
      }
      a.point_by_point =
          m.getArrays().array_layout == mgis::FieldLayout::POINT_BY_POINT;
      const auto macroscopic_gradients = m.getMacroscopicGradients();
      for (size_type c = 0; c != gsize; ++c) {
        a.macroscopic_gradients[c] = macroscopic_gradients[c];
      }
      a.gradients = m.getArrays().write(MaterialArrays::StateArray::GRADIENTS,
                                        ets, use_device);
      internals::launchGradientsKernel(k, a, use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"computeGradients: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage("computeGradients: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of computeGradients

  bool getUnknownsValues(Context& ctx,
                         mfem::Vector& v,
                         const PartialQuadratureSpace& qs,
                         const UnknownsDerivatives& d) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(v);
    static_cast<void>(qs);
    static_cast<void>(d);
    return ctx.registerErrorMessage(
        "getUnknownsValues: the kernels computing the gradients are not "
        "compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    const auto& elements = qs.getElements();
    const auto& offsets = qs.getElementOffsets();
    const auto n = std::int64_t{qs.getNumberOfIntegrationPoints()};
    if (!resizeArray(ctx, v, n * d.vdim)) {
      return false;
    }
    if (elements.Size() == 0) {
      return true;
    }
    if ((!d.has_point_values) || (d.rule == nullptr) ||
        (qs.getUniformIntegrationRule() != d.rule)) {
      return ctx.registerErrorMessage(
          "getUnknownsValues: the values of the unknowns are not computed "
          "for the integration rule of the material");
    }
    if (std::int64_t{elements.Size()} * d.nq != n) {
      return ctx.registerErrorMessage(
          "getUnknownsValues: unexpected number of integration points");
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::UnknownsValuesKernelArguments{};
      a.elements = elements.Read(use_device);
      a.offsets = offsets.Read(use_device);
      a.point_values = d.point_values.Read(use_device);
      a.values = v.Write(use_device);
      a.number_of_elements = elements.Size();
      a.nq = d.nq;
      a.vdim = d.vdim;
      internals::launchUnknownsValuesKernel(a, use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"getUnknownsValues: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage("getUnknownsValues: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of getUnknownsValues

}  // end of namespace mfem_mgis::batch
