/*!
 * \file   ElementsGeometryTest.cxx
 * \brief  Test of the geometric data of the elements used by the kernels of
 * the batch path (`mfem_mgis::batch::computeElementsGeometry`): the
 * jacobians and their determinants at the integration points are the ones
 * of the element transformations of MFEM, for meshes without nodes (kernel
 * computing the jacobians from the vertices) of all the geometries, with
 * distorted elements, and for a curved mesh (geometric factors of MFEM).
 * \date   10/10/2026
 */

#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <utility>
#include <iostream>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "mfem/fem/fe_coll.hpp"
#include "mfem/fem/fespace.hpp"
#include "mfem/fem/gridfunc.hpp"
#include "mfem/fem/intrules.hpp"
#include "mfem/mesh/mesh.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \brief moves the vertices (or the nodes) of a mesh by a smooth
 * displacement, small compared to the size of the elements, so that the
 * elements are distorted but valid
 * \param[in] x: coordinates (`point + n * component`)
 * \param[in] n: number of points
 * \param[in] dim: dimension of the space
 * \param[in] h: size of the elements
 */
static void distort(mfem_mgis::real* const x,
                    const int n,
                    const int dim,
                    const mfem_mgis::real h) {
  for (int i = 0; i != n; ++i) {
    auto s = mfem_mgis::real{0};
    for (int c = 0; c != dim; ++c) {
      s += (c + 1) * x[i + n * c];
    }
    for (int c = 0; c != dim; ++c) {
      x[i + n * c] += 0.08 * h * std::sin(7 * s + 2 * c);
    }
  }
}  // end of distort

/*!
 * \brief compares the geometric data of the elements to the element
 * transformations of MFEM
 * \param[in] mesh: mesh
 * \param[in] order: order of the finite element space of the unknowns
 * \param[in] name: name of the comparison
 */
static void compare(mfem::Mesh& mesh, const int order, const std::string& name) {
  using namespace mfem_mgis;
  const auto has_nodes = mesh.GetNodes() != nullptr;
  const auto dim = mesh.Dimension();
  auto fec = mfem::H1_FECollection(order, dim);
  auto fes = mfem::FiniteElementSpace(&mesh, &fec);
  for (const auto iorder : {1, 2, 4}) {
    const auto& ir = mfem::IntRules.Get(mesh.GetElementBaseGeometry(0), iorder);
    const auto cname = name + ", order " + std::to_string(order) +
                       ", rule of order " + std::to_string(iorder);
    auto ctx = Context{};
    auto g = batch::ElementsGeometry{};
    if (!batch::computeElementsGeometry(ctx, g, fes, ir)) {
      check(false, cname + ": " + ctx.getErrorMessage());
      continue;
    }
    check((mesh.GetNodes() != nullptr) == has_nodes,
          cname + ": the nodes of the mesh were modified");
    check((g.nq == ir.GetNPoints()) && (g.dim == dim) &&
              (g.ne == mesh.GetNE()) && (g.nd == fes.GetTypicalFE()->GetDof()),
          cname + ": sizes");
    check((g.jacobians != nullptr) && (g.determinants != nullptr),
          cname + ": no jacobians");
    if ((g.jacobians == nullptr) || (g.determinants == nullptr)) {
      continue;
    }
    check(has_nodes == (g.owned_jacobians == nullptr),
          cname + ": the jacobians are owned only for a mesh without nodes");
    check((g.jacobians->Size() == g.nq * dim * dim * g.ne) &&
              (g.determinants->Size() == g.nq * g.ne),
          cname + ": sizes of the jacobians");
    const auto* const J = g.jacobians->HostRead();
    const auto* const detJ = g.determinants->HostRead();
    const auto* const w = g.weights.HostRead();
    auto dJ = real{0};
    auto ddet = real{0};
    auto dw = real{0};
    for (int e = 0; e != mesh.GetNE(); ++e) {
      auto& T = *(mesh.GetElementTransformation(e));
      for (int q = 0; q != g.nq; ++q) {
        const auto& ip = ir.IntPoint(q);
        T.SetIntPoint(&ip);
        const auto& Jr = T.Jacobian();
        auto m = real{0};
        for (int c = 0; c != dim; ++c) {
          for (int r = 0; r != dim; ++r) {
            m = std::max(m, std::abs(Jr(r, c)));
          }
        }
        for (int c = 0; c != dim; ++c) {
          for (int r = 0; r != dim; ++r) {
            const auto v = J[q + g.nq * (r + dim * (c + dim * e))];
            dJ = std::max(dJ, std::abs(v - Jr(r, c)) / m);
          }
        }
        const auto det = Jr.Det();
        ddet = std::max(ddet, std::abs(detJ[q + g.nq * e] - det) /
                                  std::pow(m, dim));
        dw = std::max(dw, std::abs(w[q] - ip.weight));
      }
    }
    std::cout << cname << ": " << g.ne << " elements, " << g.nq
              << " points, relative differences " << dJ << " (jacobians), "
              << ddet << " (determinants), " << dw << " (weights)\n";
    check(dJ < 1e-13, cname + ": jacobians");
    check(ddet < 1e-13, cname + ": determinants");
    check(dw < 1e-15, cname + ": weights");
  }
}  // end of compare

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  using Type = mfem::Element::Type;
  auto meshes = std::vector<std::pair<std::string, mfem::Mesh>>{};
  meshes.emplace_back("segments", mfem::Mesh::MakeCartesian1D(6, 2.));
  for (const auto& [n, t] : {std::pair{"triangles", Type::TRIANGLE},
                             std::pair{"quadrangles", Type::QUADRILATERAL}}) {
    meshes.emplace_back(n, mfem::Mesh::MakeCartesian2D(3, 4, t, true, 2., 1.5));
  }
  for (const auto& [n, t] : {std::pair{"tetrahedra", Type::TETRAHEDRON},
                             std::pair{"hexahedra", Type::HEXAHEDRON},
                             std::pair{"wedges", Type::WEDGE},
                             std::pair{"pyramids", Type::PYRAMID}}) {
    meshes.emplace_back(
        n, mfem::Mesh::MakeCartesian3D(2, 3, 2, t, 1., 1.5, 1.2, true));
  }
  for (auto& [name, mesh] : meshes) {
    check(mesh.GetNodes() == nullptr, name + ": the mesh has nodes");
    compare(mesh, 1, name);
    // distorted elements (non affine transformations of the quadrangles,
    // hexahedra, wedges and pyramids)
    auto vertices = mfem::Vector{};
    mesh.GetVertices(vertices);
    distort(vertices.HostReadWrite(), mesh.GetNV(), mesh.SpaceDimension(),
            0.5);
    mesh.SetVertices(vertices);
    compare(mesh, 1, "distorted " + name);
    // unknowns of order 2: the jacobians only depend on the vertices
    compare(mesh, 2, "distorted " + name);
  }
  // curved mesh: geometric factors of MFEM
  auto curved = mfem::Mesh::MakeCartesian3D(2, 2, 2, Type::HEXAHEDRON);
  curved.SetCurvature(2, false, -1, mfem::Ordering::byNODES);
  auto& nodes = *(curved.GetNodes());
  distort(nodes.HostReadWrite(), nodes.Size() / 3, 3, 0.5);
  curved.NodesUpdated();
  compare(curved, 1, "curved hexahedra");
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
