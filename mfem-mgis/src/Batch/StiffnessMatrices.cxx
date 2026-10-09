/*!
 * \file   src/Batch/StiffnessMatrices.cxx
 * \brief  This file implements the computation of the stiffness matrices of
 * the elements of the standard mechanical behaviour integrators by kernels
 * (host part: checks and arguments of the kernels)
 * \date   09/10/2026
 */

#include <string>
#include <cstdint>
#include <exception>
#include "mfem/general/device.hpp"
#include "mfem/fem/fe/fe_base.hpp"
#include "mfem/mesh/mesh.hpp"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Batch/StiffnessMatrices.hxx"
#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES
#include "StiffnessMatricesKernels.hxx"
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

namespace mfem_mgis::batch {

  bool addStiffnessMatrices(Context& ctx,
                            mfem::Vector& m,
                            const PartialQuadratureSpace& qs,
                            const mgis::FieldView<const real>& K,
                            const ElementsGeometry& g,
                            const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(m);
    static_cast<void>(qs);
    static_cast<void>(K);
    static_cast<void>(g);
    static_cast<void>(k);
    return ctx.registerErrorMessage(
        "addStiffnessMatrices: the kernels computing the stiffness matrices "
        "are not compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    const auto& elements = qs.getElements();
    const auto& offsets = qs.getElementOffsets();
    if (elements.Size() == 0) {
      return true;
    }
    if ((g.rule == nullptr) || (qs.getUniformIntegrationRule() != g.rule)) {
      return ctx.registerErrorMessage(
          "addStiffnessMatrices: the elements of the material do not use the "
          "integration rule of the geometric data");
    }
    const auto dim = getSpaceDimension(k);
    if (g.dim != dim) {
      return ctx.registerErrorMessage(
          "addStiffnessMatrices: unexpected dimension of the space (" +
          std::to_string(g.dim) + " instead of " + std::to_string(dim) + ")");
    }
    const auto ndofs = std::int64_t{g.nd} * dim;
    if (static_cast<std::int64_t>(m.Size()) != ndofs * ndofs * g.ne) {
      return ctx.registerErrorMessage(
          "addStiffnessMatrices: unexpected size of the matrices element by "
          "element");
    }
    const auto gsize = getGradientsSize(k);
    if ((K.getNumberOfPoints() != qs.getNumberOfIntegrationPoints()) ||
        (K.getNumberOfComponents() != gsize * gsize)) {
      return ctx.registerErrorMessage(
          "addStiffnessMatrices: unexpected number of points or of components "
          "of the tangent operator blocks");
    }
    auto a = internals::StiffnessMatricesKernelArguments{};
    a.number_of_elements = elements.Size();
    a.nd = g.nd;
    a.nq = g.nq;
    a.gradients_size = static_cast<int>(gsize);
    a.tangent_operator = K.data();
    a.point_stride = static_cast<long long>(K.getPointStride());
    a.component_stride = static_cast<long long>(K.getComponentStride());
    try {
      const auto use_device = mfem::Device::IsEnabled();
      a.elements = elements.Read(use_device);
      a.offsets = offsets.Read(use_device);
      a.shape_functions_gradients = g.maps->G.Read(use_device);
      a.jacobians = g.factors->J.Read(use_device);
      a.determinants = g.factors->detJ.Read(use_device);
      a.weights = g.weights.Read(use_device);
      a.matrices = m.ReadWrite(use_device);
      internals::launchStiffnessMatricesKernel(k, a, use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"addStiffnessMatrices: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addStiffnessMatrices: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of addStiffnessMatrices

}  // end of namespace mfem_mgis::batch
