/*!
 * \file   ElementsGeometryTest.cxx
 * \brief  Test of the geometric data of the elements used by the kernels of
 * the batch path (`mfem_mgis::batch::computeElementsGeometry`): the
 * jacobians and their determinants at the integration points are the ones
 * of the element transformations of MFEM, for meshes without nodes (kernel
 * computing the jacobians from the vertices) of all the geometries, with
 * distorted and mirrored elements (negative determinants), for parallel
 * meshes refined non conformingly, for two integration rules allocated
 * successively at the same address, and for a curved mesh (geometric
 * factors of MFEM). Meshes mixing several geometries and surface meshes
 * are refused.
 * \date   10/10/2026
 */

#include <new>
#include <array>
#include <cmath>
#include <memory>
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
#ifdef MFEM_USE_MPI
#include "mfem/general/communication.hpp"
#include "mfem/fem/pfespace.hpp"
#include "mfem/mesh/pmesh.hpp"
#endif /* MFEM_USE_MPI */
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

//! \brief if true, the results are printed (main process)
static bool verbose = true;

/*!
 * \return a finite element space on the given mesh, parallel if the mesh is
 * parallel
 * \param[in] mesh: mesh
 * \param[in] fec: finite element collection
 */
static std::unique_ptr<mfem::FiniteElementSpace> makeFiniteElementSpace(
    mfem::Mesh& mesh, const mfem::FiniteElementCollection& fec) {
#ifdef MFEM_USE_MPI
  if (auto* const pmesh = dynamic_cast<mfem::ParMesh*>(&mesh)) {
    return std::make_unique<mfem::ParFiniteElementSpace>(
        pmesh, const_cast<mfem::FiniteElementCollection*>(&fec));
  }
#endif /* MFEM_USE_MPI */
  return std::make_unique<mfem::FiniteElementSpace>(
      &mesh, const_cast<mfem::FiniteElementCollection*>(&fec));
}  // end of makeFiniteElementSpace

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
 * \brief compares the geometric data of the elements for one integration
 * rule to the element transformations of MFEM
 * \param[in] mesh: mesh
 * \param[in] fes: finite element space of the unknowns
 * \param[in] ir: integration rule
 * \param[in] name: name of the comparison
 * \param[in] sign: expected sign of the determinants of the jacobians (no
 * check if null)
 */
static void compare(mfem::Mesh& mesh,
                    const mfem::FiniteElementSpace& fes,
                    const mfem::IntegrationRule& ir,
                    const std::string& name,
                    const int sign) {
  using namespace mfem_mgis;
  const auto has_nodes = mesh.GetNodes() != nullptr;
  const auto dim = mesh.Dimension();
  auto ctx = Context{};
  auto g = batch::ElementsGeometry{};
  if (!batch::computeElementsGeometry(ctx, g, fes, ir)) {
    check(false, name + ": " + ctx.getErrorMessage());
    return;
  }
  check((mesh.GetNodes() != nullptr) == has_nodes,
        name + ": the nodes of the mesh were modified");
  check((g.nq == ir.GetNPoints()) && (g.dim == dim) && (g.ne == mesh.GetNE()) &&
            (g.nd == fes.GetTypicalFE()->GetDof()),
        name + ": sizes");
  check(g.jacobians != nullptr, name + ": no jacobians");
  if (g.jacobians == nullptr) {
    return;
  }
  check((g.jacobians->J.Size() == g.nq * dim * dim * g.ne) &&
            (g.jacobians->detJ.Size() == g.nq * g.ne),
        name + ": sizes of the jacobians");
  const auto* const J = g.jacobians->J.HostRead();
  const auto* const detJ = g.jacobians->detJ.HostRead();
  const auto* const w = g.weights.HostRead();
  auto dJ = real{0};
  auto ddet = real{0};
  auto dw = real{0};
  auto signs = true;
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
      const auto det = detJ[q + g.nq * e];
      ddet = std::max(ddet, std::abs(det - Jr.Det()) / std::pow(m, dim));
      dw = std::max(dw, std::abs(w[q] - ip.weight));
      signs = signs && (sign * det >= 0);
    }
  }
  if (verbose) {
    std::cout << name << ": " << g.ne << " elements, " << g.nq
              << " points, relative differences " << dJ << " (jacobians), "
              << ddet << " (determinants), " << dw << " (weights)\n";
  }
  check(dJ < 1e-13, name + ": jacobians");
  check(ddet < 1e-13, name + ": determinants");
  check(dw < 1e-15, name + ": weights");
  check(signs, name + ": signs of the determinants");
}  // end of compare

/*!
 * \brief compares the geometric data of the elements to the element
 * transformations of MFEM for integration rules of several orders
 * \param[in] mesh: mesh
 * \param[in] order: order of the finite element space of the unknowns
 * \param[in] name: name of the comparison
 * \param[in] sign: expected sign of the determinants of the jacobians (no
 * check if null)
 */
static void compare(mfem::Mesh& mesh,
                    const int order,
                    const std::string& name,
                    const int sign = 0) {
  auto fec = mfem::H1_FECollection(order, mesh.Dimension());
  auto fes = makeFiniteElementSpace(mesh, fec);
  for (const auto iorder : {1, 2, 4}) {
    const auto& ir = mfem::IntRules.Get(mesh.GetElementBaseGeometry(0), iorder);
    compare(mesh, *fes, ir,
            name + ", order " + std::to_string(order) + ", rule of order " +
                std::to_string(iorder),
            sign);
  }
}  // end of compare

/*!
 * \brief checks the geometric data computed for two integration rules with
 * different numbers of points, successively allocated at the same address,
 * the first one being destroyed: the data must not depend on the address
 * of the rule (the caches of the finite elements of MFEM are indexed by
 * the address of the rule)
 * \param[in] mesh: mesh, without nodes
 * \param[in] name: name of the comparison
 */
static void checkRulesAtTheSameAddress(mfem::Mesh& mesh,
                                       const std::string& name) {
  alignas(mfem::IntegrationRule) unsigned char
      storage[sizeof(mfem::IntegrationRule)];
  for (const auto iorder : {1, 4}) {
    // a new space for each rule: the finite element of the unknowns, whose
    // cache is also indexed by the address of the rule, is destroyed
    auto fec = mfem::H1_FECollection(1, mesh.Dimension());
    auto fes = makeFiniteElementSpace(mesh, fec);
    auto* const ir = new (storage) mfem::IntegrationRule(
        mfem::IntRules.Get(mesh.GetElementBaseGeometry(0), iorder));
    compare(mesh, *fes, *ir,
            name + ", rule of order " + std::to_string(iorder) +
                " at a reused address",
            0);
    ir->~IntegrationRule();
  }
}  // end of checkRulesAtTheSameAddress

/*!
 * \brief checks that the meshes which are not supported are refused
 */
static void checkRefusals() {
  using namespace mfem_mgis;
  auto refused = [](mfem::Mesh& mesh, const std::string& name) {
    auto fec = mfem::H1_FECollection(1, mesh.Dimension());
    auto fes = mfem::FiniteElementSpace(&mesh, &fec);
    auto ctx = Context{};
    auto g = batch::ElementsGeometry{};
    const auto& ir = mfem::IntRules.Get(mesh.GetElementBaseGeometry(0), 2);
    check(!batch::computeElementsGeometry(ctx, g, fes, ir),
          name + ": not refused");
    if (verbose) {
      std::cout << name << ": " << ctx.getErrorMessage() << '\n';
    }
  };
  // a triangle and a quadrangle
  auto mixed = mfem::Mesh(2, 5, 2);
  for (const auto& [x, y] :
       {std::pair{0., 0.}, std::pair{1., 0.}, std::pair{1., 1.},
        std::pair{0., 1.}, std::pair{2., 0.5}}) {
    mixed.AddVertex(x, y);
  }
  mixed.AddQuad(0, 1, 2, 3);
  mixed.AddTriangle(1, 4, 2);
  mixed.FinalizeMesh();
  refused(mixed, "mesh mixing triangles and quadrangles");
  // a quadrangle in the three-dimensional space
  auto surface = mfem::Mesh(2, 4, 1, 0, 3);
  for (const auto& [x, y, z] :
       {std::array{0., 0., 0.}, std::array{1., 0., 0.}, std::array{1., 1., 1.},
        std::array{0., 1., 1.}}) {
    surface.AddVertex(x, y, z);
  }
  surface.AddQuad(0, 1, 2, 3);
  surface.FinalizeMesh();
  refused(surface, "surface mesh in the three-dimensional space");
}  // end of checkRefusals

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
#ifdef MFEM_USE_MPI
  verbose = mfem::Mpi::Root();
#endif /* MFEM_USE_MPI */
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
    compare(mesh, 1, name, 1);
    // distorted elements (non affine transformations of the quadrangles,
    // hexahedra, wedges and pyramids)
    auto vertices = mfem::Vector{};
    mesh.GetVertices(vertices);
    distort(vertices.HostReadWrite(), mesh.GetNV(), mesh.SpaceDimension(), 0.5);
    mesh.SetVertices(vertices);
    compare(mesh, 1, "distorted " + name, 1);
    // unknowns of order 2: the jacobians only depend on the vertices
    compare(mesh, 2, "distorted " + name, 1);
    checkRulesAtTheSameAddress(mesh, "distorted " + name);
#ifdef MFEM_USE_MPI
    // parallel mesh, refined non conformingly (local elements and vertices)
    if ((name == "tetrahedra") || (name == "hexahedra")) {
      auto nc = mfem::Mesh(mesh);
      nc.EnsureNCMesh(true);
      auto pmesh = mfem::ParMesh(MPI_COMM_WORLD, nc);
      auto refined = mfem::Array<int>{};
      for (int e = 0; e < std::min(pmesh.GetNE(), 3); ++e) {
        refined.Append(e);
      }
      pmesh.GeneralRefinement(refined, 1);
      check(pmesh.Nonconforming(),
            name + ": the parallel mesh is not non conforming");
      compare(pmesh, 1, "parallel non conforming distorted " + name, 1);
    }
#endif /* MFEM_USE_MPI */
    // mirrored elements: negative determinants
    for (int i = 0; i != mesh.GetNV(); ++i) {
      vertices(i) = -vertices(i);
    }
    mesh.SetVertices(vertices);
    compare(mesh, 1, "mirrored distorted " + name, -1);
  }
  // curved mesh: geometric factors of MFEM
  auto curved = mfem::Mesh::MakeCartesian3D(2, 2, 2, Type::HEXAHEDRON);
  curved.SetCurvature(2, false, -1, mfem::Ordering::byNODES);
  auto& nodes = *(curved.GetNodes());
  distort(nodes.HostReadWrite(), nodes.Size() / 3, 3, 0.5);
  curved.NodesUpdated();
  compare(curved, 1, "curved hexahedra", 1);
  checkRefusals();
  if (verbose) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
