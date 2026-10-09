/*!
 * \file   src/Batch/InnerForces.cxx
 * \brief  This file implements the computation of the inner forces of the
 * standard mechanical behaviour integrators by kernels (host part: checks
 * and arguments of the kernels)
 * \date   09/10/2026
 */

#include <limits>
#include <string>
#include <cstdint>
#include <exception>
#include "mfem/general/device.hpp"
#include "mfem/fem/fespace.hpp"
#include "mfem/mesh/mesh.hpp"
#include "mfem/fem/fe/fe_base.hpp"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Batch/InnerForces.hxx"
#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES
#include "InnerForcesKernels.hxx"
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

namespace mfem_mgis::batch {

  bool addInnerForces(Context& ctx,
                      mfem::Vector& f,
                      const PartialQuadratureSpace& qs,
                      const mgis::FieldView<const real>& s,
                      const ElementsGeometry& g,
                      const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(f);
    static_cast<void>(qs);
    static_cast<void>(s);
    static_cast<void>(g);
    static_cast<void>(k);
    return ctx.registerErrorMessage(
        "addInnerForces: the kernels computing the inner forces are not "
        "compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    const auto& elements = qs.getElements();
    const auto& offsets = qs.getElementOffsets();
    if (elements.Size() == 0) {
      return true;
    }
    if ((g.rule == nullptr) || (qs.getUniformIntegrationRule() != g.rule)) {
      return ctx.registerErrorMessage(
          "addInnerForces: the elements of the material do not use the "
          "integration rule of the geometric data");
    }
    const auto dim = getSpaceDimension(k);
    if (g.dim != dim) {
      return ctx.registerErrorMessage(
          "addInnerForces: unexpected dimension of the space (" +
          std::to_string(g.dim) + " instead of " + std::to_string(dim) + ")");
    }
    const auto size = std::int64_t{g.nd} * dim * g.ne;
    if (static_cast<std::int64_t>(f.Size()) != size) {
      return ctx.registerErrorMessage(
          "addInnerForces: unexpected size of the values element by element");
    }
    if ((s.getNumberOfPoints() != qs.getNumberOfIntegrationPoints()) ||
        (s.getNumberOfComponents() != getGradientsSize(k))) {
      return ctx.registerErrorMessage(
          "addInnerForces: unexpected number of points or of components of "
          "the thermodynamic forces");
    }
    auto a = internals::InnerForcesKernelArguments{};
    a.number_of_elements = elements.Size();
    a.nd = g.nd;
    a.nq = g.nq;
    a.dim = g.dim;
    a.thermodynamic_forces = s.data();
    a.point_stride = static_cast<long long>(s.getPointStride());
    a.component_stride = static_cast<long long>(s.getComponentStride());
    try {
      const auto use_device = mfem::Device::IsEnabled();
      a.elements = elements.Read(use_device);
      a.offsets = offsets.Read(use_device);
      a.shape_functions_gradients = g.maps->G.Read(use_device);
      a.jacobians = g.factors->J.Read(use_device);
      a.determinants = g.factors->detJ.Read(use_device);
      a.weights = g.weights.Read(use_device);
      a.values = f.ReadWrite(use_device);
      internals::launchInnerForcesKernel(k, a, use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"addInnerForces: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage("addInnerForces: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of addInnerForces

}  // end of namespace mfem_mgis::batch
