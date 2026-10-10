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
#include "mfem/fem/gridfunc.hpp"
#include "mfem/mesh/mesh.hpp"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"
#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES
#include "ElementsGeometryKernels.hxx"
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

namespace mfem_mgis::batch {

  /*!
   * \brief nodes of the transformations of the elements: the nodes of the
   * mesh if it has some, its vertices otherwise (see
   * `mfem::Mesh::GetElementTransformation`)
   */
  struct TransformationsNodes {
    //! \brief finite element of the transformations
    const mfem::FiniteElement* fe = nullptr;
    //! \brief nodes of the elements, in the order of the degrees of freedom
    mfem::Array<int> element_nodes;
    //! \brief copy of the coordinates of the vertices (mesh without nodes)
    mfem::Vector vertices;
    //! \brief coordinates of the nodes
    const mfem::Vector* coordinates = nullptr;
    //! \brief stride between two nodes in the coordinates
    int node_stride = 0;
    //! \brief stride between two components in the coordinates
    int component_stride = 0;
  };

  /*!
   * \brief gathers the nodes of the transformations of the elements
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] n: nodes
   * \param[in] mesh: mesh
   * \param[in] dim: dimension of the space
   */
  static bool getTransformationsNodes(Context& ctx,
                                      TransformationsNodes& n,
                                      const mfem::Mesh& mesh,
                                      const int dim) {
    const auto ne = mesh.GetNE();
    const auto* const nodes = mesh.GetNodes();
    if (nodes == nullptr) {
      // the finite element of the transformations, whose degrees of freedom
      // are the vertices of the elements, in their order
      n.fe =
          mfem::Mesh::GetTransformationFEforElementType(mesh.GetElementType(0));
    } else {
      const auto& nfes = *(nodes->FESpace());
      if (nfes.GetVDim() != dim) {
        return ctx.registerErrorMessage(
            "computeElementsGeometry: unexpected number of components of the "
            "nodes of the mesh");
      }
      n.fe = nfes.GetTypicalFE();
    }
    const auto nd = n.fe->GetDof();
    if (std::int64_t{ne} * nd > std::numeric_limits<int>::max()) {
      return ctx.registerErrorMessage(
          "computeElementsGeometry: too many nodes for the kernel");
    }
    n.element_nodes.SetSize(ne * nd);
    auto* const en = n.element_nodes.HostWrite();
    if (nodes == nullptr) {
      for (int e = 0; e != ne; ++e) {
        const auto* const element = mesh.GetElement(e);
        if (element->GetNVertices() != nd) {
          return ctx.registerErrorMessage(
              "computeElementsGeometry: unexpected number of vertices");
        }
        std::copy_n(element->GetVertices(), nd, en + nd * e);
      }
      mesh.GetVertices(n.vertices);
      n.coordinates = &(n.vertices);
      n.node_stride = 1;
      n.component_stride = mesh.GetNV();
      return true;
    }
    const auto& nfes = *(nodes->FESpace());
    auto dofs = mfem::Array<int>{};
    auto doftrans = mfem::DofTransformation{};
    for (int e = 0; e != ne; ++e) {
      nfes.GetElementDofs(e, dofs, doftrans);
      if ((dofs.Size() != nd) || (!doftrans.IsIdentity())) {
        return ctx.registerErrorMessage(
            "computeElementsGeometry: unsupported nodes of the mesh");
      }
      for (int i = 0; i != nd; ++i) {
        // negative degrees of freedom denote a change of orientation
        if (dofs[i] < 0) {
          return ctx.registerErrorMessage(
              "computeElementsGeometry: unsupported nodes of the mesh");
        }
        en[i + nd * e] = dofs[i];
      }
    }
    n.coordinates = nodes;
    if (nfes.GetOrdering() == mfem::Ordering::byNODES) {
      n.node_stride = 1;
      n.component_stride = nfes.GetNDofs();
    } else {
      n.node_stride = dim;
      n.component_stride = 1;
    }
    return true;
  }  // end of getTransformationsNodes

  /*!
   * \brief computes the jacobians of the element transformations by a kernel
   * from the coordinates of the nodes of the transformations
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] g: geometric data. The number of elements, the number of
   * points and the dimension must be set.
   * \param[in] mesh: mesh
   * \param[in] ir: integration rule
   */
  static bool computeJacobians(Context& ctx,
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
    if (std::int64_t{g.ne} * g.nq * g.dim * g.dim >
        std::numeric_limits<int>::max()) {
      return ctx.registerErrorMessage(
          "computeElementsGeometry: too many points for the kernel");
    }
    auto n = TransformationsNodes{};
    if (!getTransformationsNodes(ctx, n, mesh, g.dim)) {
      return false;
    }
    const auto nd = n.fe->GetDof();
    // gradients of the shape functions of the finite element of the
    // transformations at the points. They are computed here: the cache of
    // `mfem::FiniteElement::GetDofToQuad` is indexed by the address of the
    // integration rule, and would keep the gradients of a destroyed rule for
    // another rule allocated at the same address (the finite element of the
    // transformations of a mesh without nodes is a global object)
    auto gradients = mfem::Vector{};
    gradients.UseDevice(true);
    gradients.SetSize(g.nq * g.dim * nd);
    auto* const G = gradients.HostWrite();
    auto dshape = mfem::DenseMatrix(nd, g.dim);
    for (int q = 0; q != g.nq; ++q) {
      n.fe->CalcDShape(ir.IntPoint(q), dshape);
      for (int i = 0; i != nd; ++i) {
        for (int c = 0; c != g.dim; ++c) {
          G[q + g.nq * (c + g.dim * i)] = dshape(i, c);
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
    a.element_nodes = n.element_nodes.Read(use_device);
    a.coordinates = n.coordinates->Read(use_device);
    a.node_stride = n.node_stride;
    a.component_stride = n.component_stride;
    a.shape_functions_gradients = gradients.Read(use_device);
    a.jacobians = jacobians->J.Write(use_device);
    a.determinants = jacobians->detJ.Write(use_device);
    a.number_of_elements = g.ne;
    a.nodes_per_element = nd;
    a.nq = g.nq;
    a.dim = g.dim;
    internals::launchJacobiansKernel(a, use_device);
    // the temporary arrays used by the kernel are released when it is
    // finished: `cudaFree`, called by the memory manager of MFEM, waits for
    // the end of the kernels
    g.jacobians = std::move(jacobians);
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of computeJacobians

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
      // the transformations of the elements of a NURBS mesh are not given by
      // the degrees of freedom of the elements of the nodes
      if (mesh.NURBSext != nullptr) {
        return ctx.registerErrorMessage(
            "computeElementsGeometry: NURBS meshes are not supported");
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
      // the values and gradients of the shape functions of the unknowns are
      // cached by the finite element for the address of the integration
      // rule: a different number of points reveals a rule destroyed and
      // replaced by another at the same address
      if ((g.maps->FE != fe) || (g.maps->nqpt != g.nq) ||
          (g.maps->ndof != g.nd)) {
        return ctx.registerErrorMessage(
            "computeElementsGeometry: the values of the shape functions "
            "cached by the finite element do not match the integration rule");
      }
      // the mesh is not modified (see `mfem::Mesh::EnsureNodes`)
      if (!computeJacobians(ctx, g, mesh, ir)) {
        return false;
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
