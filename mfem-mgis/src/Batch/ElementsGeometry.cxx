/*!
 * \file   src/Batch/ElementsGeometry.cxx
 * \brief  This file implements the computation of the geometric data of the
 * elements used by the kernels of the batch path
 * \date   09/10/2026
 */

#include <limits>
#include <memory>
#include <string>
#include <cstdint>
#include <utility>
#include <algorithm>
#include <exception>
#include "mfem/general/device.hpp"
#include "mfem/linalg/densemat.hpp"
#include "mfem/fem/fespace.hpp"
#include "mfem/mesh/mesh.hpp"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"
#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES
#include "ElementsGeometryKernels.hxx"
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

namespace mfem_mgis::batch {

  /*!
   * \brief computes the jacobians of the element transformations of a mesh
   * without nodes from the coordinates of its vertices
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] g: geometric data. The number of elements, the number of
   * points and the dimension must be set.
   * \param[in] mesh: mesh, without nodes
   * \param[in] ir: integration rule
   */
  static bool computeJacobiansFromVertices(Context& ctx,
                                           ElementsGeometry& g,
                                           const mfem::Mesh& mesh,
                                           const mfem::IntegrationRule& ir) {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(g);
    static_cast<void>(mesh);
    static_cast<void>(ir);
    return ctx.registerErrorMessage(
        "computeElementsGeometry: the kernel computing the jacobians is not "
        "compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    // finite element of the transformations of the elements, whose degrees
    // of freedom are the vertices of the elements, in their order
    // (`mfem::Mesh::GetElementTransformation`)
    const auto* const fe =
        mfem::Mesh::GetTransformationFEforElementType(mesh.GetElementType(0));
    const auto nv = fe->GetDof();
    if ((std::int64_t{g.ne} * g.nq * g.dim * g.dim >
         std::numeric_limits<int>::max()) ||
        (std::int64_t{g.ne} * nv > std::numeric_limits<int>::max())) {
      return ctx.registerErrorMessage(
          "computeElementsGeometry: too many points for the kernel");
    }
    auto element_vertices = mfem::Array<int>(g.ne * nv);
    auto* const ev = element_vertices.HostWrite();
    for (int e = 0; e != g.ne; ++e) {
      const auto* const element = mesh.GetElement(e);
      if (element->GetNVertices() != nv) {
        return ctx.registerErrorMessage(
            "computeElementsGeometry: unexpected number of vertices");
      }
      std::copy_n(element->GetVertices(), nv, ev + nv * e);
    }
    auto vertices = mfem::Vector{};
    mesh.GetVertices(vertices);
    // gradients of the shape functions of the finite element of the
    // transformations at the points. They are computed here: the cache of
    // `mfem::FiniteElement::GetDofToQuad` is indexed by the address of the
    // integration rule, and this finite element is a global object, which
    // would keep the gradients of a destroyed rule for another rule
    // allocated at the same address
    auto gradients = mfem::Vector{};
    gradients.UseDevice(true);
    gradients.SetSize(g.nq * g.dim * nv);
    auto* const G = gradients.HostWrite();
    auto dshape = mfem::DenseMatrix(nv, g.dim);
    for (int q = 0; q != g.nq; ++q) {
      fe->CalcDShape(ir.IntPoint(q), dshape);
      for (int v = 0; v != nv; ++v) {
        for (int c = 0; c != g.dim; ++c) {
          G[q + g.nq * (c + g.dim * v)] = dshape(v, c);
        }
      }
    }
    auto jacobians = std::make_shared<ElementsJacobians>();
    jacobians->J.UseDevice(true);
    jacobians->J.SetSize(g.nq * g.dim * g.dim * g.ne);
    jacobians->detJ.UseDevice(true);
    jacobians->detJ.SetSize(g.nq * g.ne);
    const auto use_device = mfem::Device::IsEnabled();
    auto a = internals::JacobiansKernelArguments{};
    a.element_vertices = element_vertices.Read(use_device);
    a.vertices = vertices.Read(use_device);
    a.shape_functions_gradients = gradients.Read(use_device);
    a.jacobians = jacobians->J.Write(use_device);
    a.determinants = jacobians->detJ.Write(use_device);
    a.number_of_elements = g.ne;
    a.number_of_vertices = mesh.GetNV();
    a.vertices_per_element = nv;
    a.nq = g.nq;
    a.dim = g.dim;
    internals::launchJacobiansKernel(a, use_device);
    // the temporary arrays used by the kernel are released when it is
    // finished: `cudaFree`, called by the memory manager of MFEM, waits for
    // the end of the kernels
    g.jacobians = std::move(jacobians);
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of computeJacobiansFromVertices

  bool computeElementsGeometry(Context& ctx,
                               ElementsGeometry& g,
                               const mfem::FiniteElementSpace& fes,
                               const mfem::IntegrationRule& ir) noexcept {
    try {
      auto& mesh = *(fes.GetMesh());
      if (mesh.GetNumGeometries(mesh.Dimension()) > 1) {
        return ctx.registerErrorMessage(
            "computeElementsGeometry: meshes mixing several geometries are "
            "not supported");
      }
      if (mesh.SpaceDimension() != mesh.Dimension()) {
        return ctx.registerErrorMessage(
            "computeElementsGeometry: the dimension of the space must be the "
            "dimension of the mesh");
      }
      g.rule = &ir;
      g.dim = mesh.Dimension();
      g.nq = ir.GetNPoints();
      g.ne = fes.GetNE();
      if (g.ne == 0) {
        return true;
      }
      const auto* const fe = fes.GetTypicalFE();
      g.maps = &(fe->GetDofToQuad(ir, mfem::DofToQuad::FULL));
      g.nd = fe->GetDof();
      if (mesh.GetNodes() != nullptr) {
        // the geometric factors of the mesh (`mfem::Mesh::GetGeometricFactors`)
        // are destroyed if its nodes are modified: they are computed here
        auto factors =
            mfem::GeometricFactors(&mesh, ir,
                                   mfem::GeometricFactors::JACOBIANS |
                                       mfem::GeometricFactors::DETERMINANTS);
        auto jacobians = std::make_shared<ElementsJacobians>();
        jacobians->J.Swap(factors.J);
        jacobians->detJ.Swap(factors.detJ);
        g.jacobians = std::move(jacobians);
      } else {
        // the mesh is not modified (see `mfem::Mesh::EnsureNodes`): the
        // jacobians are computed from its vertices
        if (!computeJacobiansFromVertices(ctx, g, mesh, ir)) {
          return false;
        }
      }
      g.weights.UseDevice(true);
      g.weights.SetSize(g.nq);
      auto* const w = g.weights.HostWrite();
      for (int q = 0; q != g.nq; ++q) {
        w[q] = ir.IntPoint(q).weight;
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"computeElementsGeometry: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "computeElementsGeometry: unknown exception");
    }
    return true;
  }  // end of computeElementsGeometry

}  // end of namespace mfem_mgis::batch
