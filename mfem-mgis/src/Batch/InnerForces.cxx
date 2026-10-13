/*!
 * \file   src/Batch/InnerForces.cxx
 * \brief  This file implements the computation of the inner forces of the
 * standard behaviour integrators (mechanics, heat transfer, micromorphic
 * damage) by kernels (host part: checks and arguments of the kernels)
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
#include "InnerForcesKernels.hxx"

namespace mfem_mgis::batch {

  /*!
   * \brief implementation of `addInnerForces` and of
   * `addMaterialInnerForces`
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] f: values element by element, or by element of the
   * material
   * \param[in] qs: quadrature space of the material
   * \param[in] s: thermodynamic forces in the global frame
   * \param[in] g: geometric data
   * \param[in] k: kind of the gradients
   * \param[in] by_material_element: if true, the values are stored by
   * element of the material
   * \param[in] name: name of the calling function, for the messages
   */
  static bool addInnerForcesImplementation(Context& ctx,
                                           mfem::Vector& f,
                                           const PartialQuadratureSpace& qs,
                                           const mgis::FieldView<const real>& s,
                                           const ElementsGeometry& g,
                                           const GradientsKind k,
                                           const bool by_material_element,
                                           const char* const name) noexcept {
    const auto& elements = qs.getElements();
    const auto& offsets = qs.getElementOffsets();
    if (elements.Size() == 0) {
      return true;
    }
    if ((g.rule == nullptr) || (g.maps == nullptr) ||
        (g.jacobians == nullptr) ||
        (qs.getUniformIntegrationRule() != g.rule)) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the elements of the material do not use the integration rule of "
          "the geometric data");
    }
    const auto dim = getSpaceDimension(k);
    if (g.dim != dim) {
      return ctx.registerErrorMessage(
          std::string{name} + ": unexpected dimension of the space (" +
          std::to_string(g.dim) + " instead of " + std::to_string(dim) + ")");
    }
    const auto size = std::int64_t{g.nd} * getUnknownsSize(k) *
                      (by_material_element ? elements.Size() : g.ne);
    if (static_cast<std::int64_t>(f.Size()) != size) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": unexpected size of the values element by element");
    }
    if ((s.getNumberOfPoints() != qs.getNumberOfIntegrationPoints()) ||
        (s.getNumberOfComponents() != getGradientsSize(k))) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": unexpected number of points or of components of the "
          "thermodynamic forces");
    }
    auto a = internals::InnerForcesKernelArguments{};
    a.number_of_elements = elements.Size();
    a.nd = g.nd;
    a.nq = g.nq;
    a.dim = g.dim;
    a.by_material_element = by_material_element;
    a.thermodynamic_forces = s.data();
    a.point_stride = static_cast<long long>(s.getPointStride());
    a.component_stride = static_cast<long long>(s.getComponentStride());
    try {
      const auto use_device = mfem::Device::IsEnabled();
      a.elements = elements.Read(use_device);
      a.offsets = offsets.Read(use_device);
      a.shape_functions_gradients = g.maps->G.Read(use_device);
      a.shape_functions = g.maps->B.Read(use_device);
      a.inverse_jacobians = g.jacobians->invJ.Read(use_device);
      a.determinants = g.jacobians->detJ.Read(use_device);
      a.weights = g.weights.Read(use_device);
      a.values = f.ReadWrite(use_device);
      internals::launchInnerForcesKernel(k, a, use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{name} + ": " + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(std::string{name} +
                                      ": unknown exception");
    }
    return true;
  }  // end of addInnerForcesImplementation

  bool addInnerForces(Context& ctx,
                      mfem::Vector& f,
                      const PartialQuadratureSpace& qs,
                      const mgis::FieldView<const real>& s,
                      const ElementsGeometry& g,
                      const GradientsKind k) noexcept {
    return addInnerForcesImplementation(ctx, f, qs, s, g, k, false,
                                        "addInnerForces");
  }  // end of addInnerForces

  bool addMaterialInnerForces(Context& ctx,
                              mfem::Vector& f,
                              const PartialQuadratureSpace& qs,
                              const mgis::FieldView<const real>& s,
                              const ElementsGeometry& g,
                              const GradientsKind k) noexcept {
    return addInnerForcesImplementation(ctx, f, qs, s, g, k, true,
                                        "addMaterialInnerForces");
  }  // end of addMaterialInnerForces

}  // end of namespace mfem_mgis::batch
