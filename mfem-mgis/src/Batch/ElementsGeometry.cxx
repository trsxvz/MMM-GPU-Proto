/*!
 * \file   src/Batch/ElementsGeometry.cxx
 * \brief  This file implements the computation of the geometric data of the
 * elements used by the kernels of the batch path
 * \date   09/10/2026
 */

#include <memory>
#include <string>
#include <exception>
#include "mfem/fem/fespace.hpp"
#include "mfem/mesh/mesh.hpp"
#include "mfem/fem/fe_coll.hpp"
#include "mfem/fem/gridfunc.hpp"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"

namespace mfem_mgis::batch {

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
      const auto flags = mfem::GeometricFactors::JACOBIANS |
                         mfem::GeometricFactors::DETERMINANTS;
      if (mesh.GetNodes() != nullptr) {
        g.factors = mesh.GetGeometricFactors(ir, flags);
      } else {
        // the mesh is not modified (see `mfem::Mesh::EnsureNodes`): the
        // jacobians are computed from a copy of its vertices
        auto fec = mfem::H1_FECollection(1, g.dim);
        auto nfes = mfem::FiniteElementSpace(&mesh, &fec, mesh.SpaceDimension(),
                                             mfem::Ordering::byVDIM);
        auto nodes = mfem::GridFunction(&nfes);
        mesh.GetNodes(nodes);
        g.owned_factors =
            std::make_shared<mfem::GeometricFactors>(nodes, ir, flags);
        g.factors = g.owned_factors.get();
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
