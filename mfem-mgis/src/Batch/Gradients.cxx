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
#include "mfem/fem/quadinterpolator.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES
#include "GradientsKernels.hxx"
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

namespace mfem_mgis::batch {

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
                                  const mfem::IntegrationRule& ir) noexcept {
    try {
      const auto& mesh = *(fes.GetMesh());
      if (mesh.GetNumGeometries(mesh.Dimension()) > 1) {
        return ctx.registerErrorMessage(
            "computeUnknownsDerivatives: meshes mixing several geometries "
            "are not supported");
      }
      const auto* const qi = fes.GetQuadratureInterpolator(ir);
      if (qi == nullptr) {
        return ctx.registerErrorMessage(
            "computeUnknownsDerivatives: the finite element space is not "
            "supported by the quadrature interpolator of MFEM");
      }
      const auto* const r =
          fes.GetElementRestriction(mfem::ElementDofOrdering::NATIVE);
      if (!resizeArray(ctx, d.element_values, r->Height())) {
        return false;
      }
      r->Mult(u, d.element_values);
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
      // the interpolator is shared: its options are set before each use
      qi->SetOutputLayout(mfem::QVectorLayout::byVDIM);
      qi->DisableTensorProducts();
      qi->Derivatives(d.element_values, d.values);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"computeUnknownsDerivatives: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "computeUnknownsDerivatives: unknown exception");
    }
    return true;
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
        break;
    }
    return 5;
  }  // end of getGradientsSize

  int getSpaceDimension(const GradientsKind k) noexcept {
    return ((k == GradientsKind::SMALL_STRAIN_3D) ||
            (k == GradientsKind::DEFORMATION_GRADIENT_3D))
               ? 3
               : 2;
  }  // end of getSpaceDimension

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
    const auto& qs = m.getPartialQuadratureSpace();
    const auto& elements = qs.getElements();
    const auto& offsets = qs.getElementOffsets();
    if (elements.Size() == 0) {
      return true;
    }
    if ((d.rule == nullptr) || (qs.getUniformIntegrationRule() != d.rule) ||
        (g.rule != d.rule) || (g.factors == nullptr) || (g.nq != d.nq) ||
        (g.ne != d.ne)) {
      return ctx.registerErrorMessage(
          "computeGradients: the elements of the material do not use the "
          "integration rule of the derivatives of the unknowns and of the "
          "geometric data");
    }
    const auto gsize = getGradientsSize(k);
    const auto dim = getSpaceDimension(k);
    if ((d.vdim != dim) || (d.dim != dim)) {
      return ctx.registerErrorMessage(
          "computeGradients: unexpected number of components of the "
          "unknowns or dimension of the space (" +
          std::to_string(d.vdim) + ", " + std::to_string(d.dim) +
          "), the expected value being " + std::to_string(dim));
    }
    if (m.s1.gradients_stride != gsize) {
      return ctx.registerErrorMessage(
          "computeGradients: unexpected number of components of the "
          "gradients (" +
          std::to_string(m.s1.gradients_stride) + " instead of " +
          std::to_string(gsize) + ")");
    }
    if (static_cast<std::int64_t>(elements.Size()) * d.nq !=
        static_cast<std::int64_t>(m.n)) {
      return ctx.registerErrorMessage(
          "computeGradients: unexpected number of integration points");
    }
    const auto macroscopic_gradients = m.getMacroscopicGradients();
    auto a = internals::GradientsKernelArguments{};
    a.number_of_elements = elements.Size();
    a.nq = d.nq;
    a.vdim = d.vdim;
    a.dim = d.dim;
    a.number_of_points = static_cast<int>(m.n);
    a.gradients_size = static_cast<int>(gsize);
    a.point_by_point =
        m.getArrays().array_layout == mgis::FieldLayout::POINT_BY_POINT;
    for (size_type c = 0; c != gsize; ++c) {
      a.macroscopic_gradients[c] = macroscopic_gradients[c];
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      a.elements = elements.Read(use_device);
      a.offsets = offsets.Read(use_device);
      a.derivatives = d.values.Read(use_device);
      a.jacobians = g.factors->J.Read(use_device);
      a.determinants = g.factors->detJ.Read(use_device);
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

}  // end of namespace mfem_mgis::batch
