/*!
 * \file   src/Batch/Faltus.cxx
 * \brief  This file implements the computations of the regularization of
 * Faltus et al. by kernels of MFEM (see `MFEMMGIS/Batch/Faltus.hxx`)
 * \date   11/10/2026
 */

#include <limits>
#include <string>
#include <cstdint>
#include <exception>
#include "mfem/general/device.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Batch/FBar.hxx"
#include "MFEMMGIS/Batch/Faltus.hxx"
#include "FaltusKernels.hxx"

namespace mfem_mgis::batch {

  /*!
   * \brief copies the values of an array of the memory of MFEM, on the
   * device if a device is enabled
   * \param[out] d: destination
   * \param[in] s: source
   */
  static void copyFaltusArray(mfem::Vector& d, const mfem::Vector& s) {
    d.UseDevice(true);
    d.SetSize(s.Size());
    if (s.Size() != 0) {
      d = s;
    }
  }  // end of copyFaltusArray

  void updateFaltusData(FaltusData& data) {
    auto& b = data.beginning_of_time_step;
    copyFaltusArray(b.centers, data.centers);
    b.valid = data.valid;
    b.defined = true;
  }  // end of updateFaltusData

  void revertFaltusData(FaltusData& data) {
    const auto& b = data.beginning_of_time_step;
    if (b.defined) {
      copyFaltusArray(data.centers, b.centers);
      data.valid = b.valid;
    } else {
      data.centers.SetSize(0);
      data.valid = true;
    }
  }  // end of revertFaltusData

  bool isFaltusSupportedByKernels(const GradientsKind k) noexcept {
    return (k == GradientsKind::DEFORMATION_GRADIENT_3D) ||
           (k == GradientsKind::DEFORMATION_GRADIENT_2D);
  }  // end of isFaltusSupportedByKernels

  /*!
   * \brief checks the geometric data of the integration points and of the
   * center of the elements of a material, and returns the arguments of the
   * kernels
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] a: arguments of the kernels
   * \param[in] qs: quadrature space of the material
   * \param[in] alpha: penalization coefficient
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] g0: geometric data of the center of the elements
   * \param[in] k: kind of the gradients
   * \param[in] use_device: if true, the kernels run on the device
   * \param[in] name: name of the calling function
   * \param[in] threads_per_element: number of threads per element of the
   * kernel
   */
  static bool getFaltusGeometryKernelArguments(
      Context& ctx,
      internals::FaltusGeometryKernelArguments& a,
      const PartialQuadratureSpace& qs,
      const real alpha,
      const ElementsGeometry& g,
      const ElementsGeometry& g0,
      const GradientsKind k,
      const bool use_device,
      const char* const name,
      const std::int64_t threads_per_element) {
    if (!isFaltusSupportedByKernels(k)) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the regularization of Faltus et al. requires a deformation "
          "gradient in 3D or in 2D");
    }
    const auto& elements = qs.getElements();
    if (elements.Size() == 0) {
      return true;
    }
    const auto dim = getSpaceDimension(k);
    if ((g.rule == nullptr) || (g.maps == nullptr) ||
        (g.jacobians == nullptr) ||
        (qs.getUniformIntegrationRule() != g.rule) || (g.dim != dim)) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the elements of the material do not use the integration rule of "
          "the geometric data");
    }
    if ((g0.rule != getElementsCenterIntegrationRule(qs)) || (g0.nq != 1) ||
        (g0.maps == nullptr) || (g0.jacobians == nullptr) || (g0.dim != dim) ||
        (g0.nd != g.nd) || (g0.ne != g.ne)) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": invalid geometric data of the center of the elements");
    }
    if (std::int64_t{elements.Size()} * threads_per_element >
        std::numeric_limits<int>::max()) {
      return ctx.registerErrorMessage(std::string{name} +
                                      ": too many elements for the kernels");
    }
    a.elements = elements.Read(use_device);
    a.shape_functions_gradients = g.maps->G.Read(use_device);
    a.inverse_jacobians = g.jacobians->invJ.Read(use_device);
    a.determinants = g.jacobians->detJ.Read(use_device);
    a.weights = g.weights.Read(use_device);
    a.center_shape_functions_gradients = g0.maps->G.Read(use_device);
    a.center_inverse_jacobians = g0.jacobians->invJ.Read(use_device);
    a.alpha = alpha;
    a.number_of_elements = elements.Size();
    a.nd = g.nd;
    a.nq = g.nq;
    return true;
  }  // end of getFaltusGeometryKernelArguments

  //! \return the number of degrees of freedom of an element
  static std::int64_t getNumberOfElementDegreesOfFreedom(
      const ElementsGeometry& g, const GradientsKind k) {
    return std::int64_t{g.nd} * getUnknownsSize(k);
  }  // end of getNumberOfElementDegreesOfFreedom

  bool computeFaltusCenters(Context& ctx,
                            FaltusData& data,
                            const Material& m,
                            const UnknownsDerivatives* const d0,
                            const ElementsGeometry& g0,
                            const GradientsKind k) noexcept {
    try {
      if (!isFaltusSupportedByKernels(k)) {
        return ctx.registerErrorMessage(
            "computeFaltusCenters: the regularization of Faltus et al. "
            "requires a deformation gradient in 3D or in 2D");
      }
      const auto& qs = m.getPartialQuadratureSpace();
      const auto& elements = qs.getElements();
      const auto gsize = static_cast<std::int64_t>(getGradientsSize(k));
      const auto size = std::int64_t{elements.Size()} * gsize;
      if (size > std::numeric_limits<int>::max()) {
        return ctx.registerErrorMessage(
            "computeFaltusCenters: array too large for MFEM");
      }
      // valid only after a successful computation
      data.valid = false;
      data.centers.UseDevice(true);
      data.centers.SetSize(static_cast<int>(size));
      if (elements.Size() == 0) {
        data.valid = true;
        return true;
      }
      const auto* const r0 = getElementsCenterIntegrationRule(qs);
      const auto vdim = getUnknownsSize(k);
      if ((g0.rule != r0) || (g0.nq != 1) || (g0.maps == nullptr) ||
          (g0.jacobians == nullptr) || (g0.dim != getSpaceDimension(k))) {
        return ctx.registerErrorMessage(
            "computeFaltusCenters: invalid geometric data of the center of "
            "the elements");
      }
      if ((d0 != nullptr) && ((d0->rule != r0) || (d0->nq != 1) ||
                              (d0->vdim != vdim) || (d0->ne != g0.ne))) {
        return ctx.registerErrorMessage(
            "computeFaltusCenters: invalid derivatives of the unknowns at the "
            "center of the elements");
      }
      const auto use_device = mfem::Device::IsEnabled();
      const auto macroscopic_gradients = m.getMacroscopicGradients();
      auto a = internals::FaltusCentersKernelArguments{};
      a.elements = elements.Read(use_device);
      a.derivatives = (d0 != nullptr) ? d0->values.Read(use_device) : nullptr;
      a.inverse_jacobians = g0.jacobians->invJ.Read(use_device);
      a.centers = data.centers.Write(use_device);
      a.number_of_elements = elements.Size();
      for (std::int64_t c = 0; c != gsize; ++c) {
        a.macroscopic_gradients[c] = macroscopic_gradients[c];
      }
      internals::launchFaltusCentersKernel(k, a, use_device);
      data.valid = true;
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"computeFaltusCenters: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "computeFaltusCenters: unknown exception");
    }
    return true;
  }  // end of computeFaltusCenters

  bool addFaltusInnerForces(Context& ctx,
                            mfem::Vector& f,
                            const FaltusData& data,
                            const Material& m,
                            const real alpha,
                            const ElementsGeometry& g,
                            const ElementsGeometry& g0,
                            const GradientsKind k) noexcept {
    try {
      const auto& qs = m.getPartialQuadratureSpace();
      const auto& elements = qs.getElements();
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::FaltusInnerForcesKernelArguments{};
      if (!getFaltusGeometryKernelArguments(ctx, a.geometry, qs, alpha, g, g0,
                                            k, use_device,
                                            "addFaltusInnerForces", 1)) {
        return false;
      }
      if (elements.Size() == 0) {
        return true;
      }
      const auto gsize = static_cast<std::int64_t>(getGradientsSize(k));
      if ((!data.valid) || (static_cast<std::int64_t>(data.centers.Size()) !=
                            std::int64_t{elements.Size()} * gsize)) {
        return ctx.registerErrorMessage(
            "addFaltusInnerForces: the deformation gradients at the center of "
            "the elements are not computed by kernels for the current "
            "gradients (GradientsComputation set to Kernels)");
      }
      if ((m.s1.gradients_stride != getGradientsSize(k)) ||
          (static_cast<std::int64_t>(f.Size()) !=
           getNumberOfElementDegreesOfFreedom(g, k) * g.ne)) {
        return ctx.registerErrorMessage(
            "addFaltusInnerForces: unexpected number of components of the "
            "gradients or size of the values element by element");
      }
      a.offsets = qs.getElementOffsets().Read(use_device);
      a.gradients = m.getArrays().read(MaterialArrays::StateArray::GRADIENTS,
                                       ets, use_device);
      a.point_by_point =
          m.getArrays().array_layout == mgis::FieldLayout::POINT_BY_POINT;
      a.number_of_points = qs.getNumberOfIntegrationPoints();
      a.centers = data.centers.Read(use_device);
      a.values = f.ReadWrite(use_device);
      internals::launchFaltusInnerForcesKernel(k, a, use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"addFaltusInnerForces: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addFaltusInnerForces: unknown exception");
    }
    return true;
  }  // end of addFaltusInnerForces

  bool addFaltusStiffnessMatrices(Context& ctx,
                                  mfem::Vector& m,
                                  const PartialQuadratureSpace& qs,
                                  const real alpha,
                                  const ElementsGeometry& g,
                                  const ElementsGeometry& g0,
                                  const GradientsKind k) noexcept {
    try {
      const auto use_device = mfem::Device::IsEnabled();
      const auto n = getNumberOfElementDegreesOfFreedom(g, k);
      auto a = internals::FaltusGeometryKernelArguments{};
      // one thread per pair of nodes of an element
      if (!getFaltusGeometryKernelArguments(
              ctx, a, qs, alpha, g, g0, k, use_device,
              "addFaltusStiffnessMatrices", std::int64_t{g.nd} * g.nd)) {
        return false;
      }
      if (a.number_of_elements == 0) {
        return true;
      }
      if (static_cast<std::int64_t>(m.Size()) != n * n * g.ne) {
        return ctx.registerErrorMessage(
            "addFaltusStiffnessMatrices: unexpected size of the matrices "
            "element by element");
      }
      internals::launchFaltusMatricesKernel(k, a, m.ReadWrite(use_device),
                                            use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addFaltusStiffnessMatrices: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addFaltusStiffnessMatrices: unknown exception");
    }
    return true;
  }  // end of addFaltusStiffnessMatrices

  bool addFaltusStiffnessMatrices(Context& ctx,
                                  const AssembledStiffnessMatrix& m,
                                  const PartialQuadratureSpace& qs,
                                  const real alpha,
                                  const ElementsGeometry& g,
                                  const ElementsGeometry& g0,
                                  const GradientsKind k) noexcept {
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::FaltusGeometryKernelArguments{};
      // one thread per node of an element
      if (!getFaltusGeometryKernelArguments(
              ctx, a, qs, alpha, g, g0, k, use_device,
              "addFaltusStiffnessMatrices", g.nd)) {
        return false;
      }
      if (a.number_of_elements == 0) {
        return true;
      }
      if (!internals::checkAssembledStiffnessMatrix(
              ctx, m, g.ne, getNumberOfElementDegreesOfFreedom(g, k),
              "addFaltusStiffnessMatrices")) {
        return false;
      }
      auto counters = mfem::Array<int>{};
      auto o = internals::SparseMatrixOutput{};
      internals::getSparseMatrixOutput(o, counters, m, use_device);
      internals::launchFaltusMatricesKernel(k, a, o, use_device);
      if (!internals::checkSparseMatrixOutputCounters(
              ctx, counters, m, "addFaltusStiffnessMatrices")) {
        return false;
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addFaltusStiffnessMatrices: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addFaltusStiffnessMatrices: unknown exception");
    }
    return true;
  }  // end of addFaltusStiffnessMatrices

  bool addFaltusStiffnessMatricesDiagonals(Context& ctx,
                                           mfem::Vector& d,
                                           const PartialQuadratureSpace& qs,
                                           const real alpha,
                                           const ElementsGeometry& g,
                                           const ElementsGeometry& g0,
                                           const GradientsKind k) noexcept {
    try {
      const auto use_device = mfem::Device::IsEnabled();
      const auto n = getNumberOfElementDegreesOfFreedom(g, k);
      auto a = internals::FaltusGeometryKernelArguments{};
      // one thread per node of an element
      if (!getFaltusGeometryKernelArguments(
              ctx, a, qs, alpha, g, g0, k, use_device,
              "addFaltusStiffnessMatricesDiagonals", g.nd)) {
        return false;
      }
      if (a.number_of_elements == 0) {
        return true;
      }
      if (static_cast<std::int64_t>(d.Size()) != n * g.ne) {
        return ctx.registerErrorMessage(
            "addFaltusStiffnessMatricesDiagonals: unexpected size of the "
            "values element by element");
      }
      internals::launchFaltusDiagonalsKernel(k, a, d.ReadWrite(use_device),
                                             use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addFaltusStiffnessMatricesDiagonals: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addFaltusStiffnessMatricesDiagonals: unknown exception");
    }
    return true;
  }  // end of addFaltusStiffnessMatricesDiagonals

  bool addFaltusTangentOperatorAction(Context& ctx,
                                      mfem::Vector& y,
                                      const PartialQuadratureSpace& qs,
                                      const mfem::Vector& x,
                                      const real alpha,
                                      const ElementsGeometry& g,
                                      const ElementsGeometry& g0,
                                      const GradientsKind k) noexcept {
    try {
      const auto use_device = mfem::Device::IsEnabled();
      const auto n = getNumberOfElementDegreesOfFreedom(g, k);
      auto a = internals::FaltusGeometryKernelArguments{};
      // one thread per element
      if (!getFaltusGeometryKernelArguments(
              ctx, a, qs, alpha, g, g0, k, use_device,
              "addFaltusTangentOperatorAction", 1)) {
        return false;
      }
      if (a.number_of_elements == 0) {
        return true;
      }
      if ((static_cast<std::int64_t>(x.Size()) != n * g.ne) ||
          (static_cast<std::int64_t>(y.Size()) != n * g.ne)) {
        return ctx.registerErrorMessage(
            "addFaltusTangentOperatorAction: unexpected size of the values "
            "element by element");
      }
      internals::launchFaltusActionKernel(k, a, x.Read(use_device),
                                          y.ReadWrite(use_device), use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addFaltusTangentOperatorAction: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addFaltusTangentOperatorAction: unknown exception");
    }
    return true;
  }  // end of addFaltusTangentOperatorAction

}  // end of namespace mfem_mgis::batch
