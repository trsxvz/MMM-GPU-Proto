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

  namespace {

    //! \brief restores the options of a quadrature interpolator at its
    //! destruction
    struct QuadratureInterpolatorOptionsGuard {
      //! \param[in] q: quadrature interpolator
      explicit QuadratureInterpolatorOptionsGuard(
          const mfem::QuadratureInterpolator& q)
          : qi(q),
            layout(q.GetOutputLayout()),
            tensor_products(q.UsesTensorProducts()) {}
      QuadratureInterpolatorOptionsGuard(
          const QuadratureInterpolatorOptionsGuard&) = delete;
      QuadratureInterpolatorOptionsGuard& operator=(
          const QuadratureInterpolatorOptionsGuard&) = delete;
      //! \brief destructor
      ~QuadratureInterpolatorOptionsGuard() {
        this->qi.SetOutputLayout(this->layout);
        this->qi.DisableTensorProducts(!this->tensor_products);
      }

     private:
      //! \brief quadrature interpolator
      const mfem::QuadratureInterpolator& qi;
      //! \brief output layout
      const mfem::QVectorLayout layout;
      //! \brief if true, the interpolator uses tensor products
      const bool tensor_products;
    };

  }  // end of namespace

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
      // the interpolator is shared: its options are set before its use and
      // restored after, even if an exception is thrown
      const auto options = QuadratureInterpolatorOptionsGuard{*qi};
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
        (g.rule != d.rule) || (g.factors == nullptr) || (g.nq != d.nq) ||
        (g.ne != d.ne)) {
      return ctx.registerErrorMessage(
          name +
          ": the elements of the material do not use the integration rule "
          "of the derivatives of the unknowns and of the geometric data");
    }
    const auto dim = getSpaceDimension(k);
    if ((d.vdim != dim) || (d.dim != dim)) {
      return ctx.registerErrorMessage(
          name +
          ": unexpected number of components of the unknowns or dimension "
          "of the space (" +
          std::to_string(d.vdim) + ", " + std::to_string(d.dim) +
          "), the expected value being " + std::to_string(dim));
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
    a.jacobians = g.factors->J.Read(use_device);
    a.determinants = g.factors->detJ.Read(use_device);
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

  bool computeGradientsIncrements(Context& ctx,
                                  mfem::Vector& v,
                                  const PartialQuadratureSpace& qs,
                                  const UnknownsDerivatives& d,
                                  const ElementsGeometry& g,
                                  const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(v);
    static_cast<void>(qs);
    static_cast<void>(d);
    static_cast<void>(g);
    static_cast<void>(k);
    return ctx.registerErrorMessage(
        "computeGradientsIncrements: the kernels computing the gradients are "
        "not compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::GradientsKernelArguments{};
      auto error = false;
      if (!getGradientsKernelArguments(ctx, a, error,
                                       "computeGradientsIncrements", qs, d, g,
                                       k, use_device)) {
        if (!error) {
          v.SetSize(0);
        }
        return !error;
      }
      const auto size =
          std::int64_t{a.number_of_points} * std::int64_t{a.gradients_size};
      if (size > std::numeric_limits<int>::max()) {
        return ctx.registerErrorMessage(
            "computeGradientsIncrements: too many values for MFEM");
      }
      v.UseDevice(true);
      v.SetSize(static_cast<int>(size));
      // stored point by point, without macroscopic gradients
      a.point_by_point = true;
      a.gradients = v.Write(use_device);
      internals::launchGradientsKernel(k, a, use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"computeGradientsIncrements: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "computeGradientsIncrements: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of computeGradientsIncrements

}  // end of namespace mfem_mgis::batch
