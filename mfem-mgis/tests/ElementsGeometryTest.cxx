/*!
 * \file   ElementsGeometryTest.cxx
 * \brief  Test of the geometric data of the elements used by the kernels of
 * the batch path (`mfem_mgis::batch::computeElementsGeometry`): the
 * jacobians and their determinants at the integration points are the ones
 * of the element transformations of MFEM, for meshes without nodes (kernel
 * computing the jacobians from the vertices) of all the geometries, with
 * distorted and mirrored elements (negative determinants), for parallel
 * meshes refined non conformingly, for two integration rules allocated
 * successively at the same address, for curved meshes (kernel computing
 * the jacobians from the nodes, ordered by nodes or by components, of
 * order 2 or 3, discontinuous for a periodic mesh) with a rule whose
 * points are not ordered as a tensor product, and the jacobians are kept
 * when the nodes of the mesh are modified. The values of the shape
 * functions of a destroyed rule, cached by the finite element, are
 * detected. Meshes mixing several geometries, surface meshes and NURBS
 * meshes are refused.
 * \date   10/10/2026
 */

#include <new>
#include <array>
#include <sstream>
#include <cmath>
#include <memory>
#include <tuple>
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
#include "MFEMMGIS/Batch/Gradients.hxx"
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
 * \param[in] x: coordinates (`point * ps + component * cs`)
 * \param[in] n: number of points
 * \param[in] dim: dimension of the space
 * \param[in] h: size of the elements
 * \param[in] ps: stride between two points
 * \param[in] cs: stride between two components
 */
static void distort(mfem_mgis::real* const x,
                    const int n,
                    const int dim,
                    const mfem_mgis::real h,
                    const int ps,
                    const int cs) {
  for (int i = 0; i != n; ++i) {
    auto s = mfem_mgis::real{0};
    for (int c = 0; c != dim; ++c) {
      s += (c + 1) * x[i * ps + c * cs];
    }
    for (int c = 0; c != dim; ++c) {
      x[i * ps + c * cs] += 0.08 * h * std::sin(7 * s + 2 * c);
    }
  }
}  // end of distort

/*!
 * \brief applies a transformation to the nodes of a curved mesh
 * \param[in, out] mesh: mesh with nodes
 * \param[in] f: transformation of the coordinates `(x, n, dim, ps, cs)` (see
 * `distort`)
 */
template <typename Transformation>
static void transformNodes(mfem::Mesh& mesh, const Transformation& f) {
  auto& nodes = *(mesh.GetNodes());
  const auto& nfes = *(nodes.FESpace());
  const auto n = nfes.GetNDofs();
  const auto dim = nfes.GetVDim();
  const auto by_nodes = nfes.GetOrdering() == mfem::Ordering::byNODES;
  f(nodes.HostReadWrite(), n, dim, by_nodes ? 1 : dim, by_nodes ? n : 1);
  mesh.NodesUpdated();
}  // end of transformNodes

//! \brief changes the sign of the first coordinate of the points
static void mirror(
    mfem_mgis::real* const x, const int n, const int, const int ps, const int) {
  for (int i = 0; i != n; ++i) {
    x[i * ps] = -x[i * ps];
  }
}  // end of mirror

//! \return the integration rule with the points in the reverse order
static mfem::IntegrationRule reverse(const mfem::IntegrationRule& ir) {
  const auto n = ir.GetNPoints();
  auto r = mfem::IntegrationRule(n);
  for (int i = 0; i != n; ++i) {
    r.IntPoint(i) = ir.IntPoint(n - 1 - i);
  }
  return r;
}  // end of reverse

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
  check((g.maps != nullptr) && (g.maps->nqpt == g.nq) && (g.maps->ndof == g.nd),
        name + ": shape functions of the unknowns");
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
 * \brief checks that the geometric data are refused if the values of the
 * shape functions cached by the finite element of the unknowns belong to a
 * destroyed integration rule: two rules with different numbers of points
 * allocated successively at the same address, for the same finite element
 * space
 * \param[in] mesh: mesh
 * \param[in] name: name of the check
 */
static void checkStaleShapeFunctions(mfem::Mesh& mesh,
                                     const std::string& name) {
  using namespace mfem_mgis;
  alignas(mfem::IntegrationRule) unsigned char
      storage[sizeof(mfem::IntegrationRule)];
  auto fec = mfem::H1_FECollection(1, mesh.Dimension());
  auto fes = makeFiniteElementSpace(mesh, fec);
  for (const auto iorder : {1, 4}) {
    auto* const ir = new (storage) mfem::IntegrationRule(
        mfem::IntRules.Get(mesh.GetElementBaseGeometry(0), iorder));
    auto ctx = Context{};
    auto g = batch::ElementsGeometry{};
    const auto success = batch::computeElementsGeometry(ctx, g, *fes, *ir);
    // derivatives of the unknowns, which also use the values of the shape
    // functions cached by the finite element
    auto u = mfem::Vector(fes->GetVSize());
    u = 1;
    auto d = batch::UnknownsDerivatives{};
    auto ctx2 = Context{};
    const auto success2 =
        batch::computeUnknownsDerivatives(ctx2, d, *fes, u, *ir, true);
    if (iorder == 1) {
      check(success, name + ": " + ctx.getErrorMessage());
      check(success2, name + ": " + ctx2.getErrorMessage());
    } else {
      check(!success, name + ": stale shape functions not detected");
      check(!success2, name +
                           ": stale shape functions not detected by "
                           "computeUnknownsDerivatives");
      if (verbose) {
        std::cout << name << ": " << ctx.getErrorMessage() << '\n'
                  << name << ": " << ctx2.getErrorMessage() << '\n';
      }
    }
    ir->~IntegrationRule();
  }
}  // end of checkStaleShapeFunctions

/*!
 * \brief checks that the jacobians computed for a curved mesh are kept when
 * its nodes are modified (the geometric factors of the mesh are then
 * destroyed by MFEM)
 * \param[in] mesh: mesh with nodes
 * \param[in] name: name of the check
 */
static void checkJacobiansLifetime(mfem::Mesh& mesh, const std::string& name) {
  using namespace mfem_mgis;
  auto fec = mfem::H1_FECollection(1, mesh.Dimension());
  auto fes = makeFiniteElementSpace(mesh, fec);
  const auto& ir = mfem::IntRules.Get(mesh.GetElementBaseGeometry(0), 2);
  auto ctx = Context{};
  auto g = batch::ElementsGeometry{};
  if (!batch::computeElementsGeometry(ctx, g, *fes, ir)) {
    check(false, name + ": " + ctx.getErrorMessage());
    return;
  }
  const auto* const J = g.jacobians->J.HostRead();
  const auto before = std::vector<real>(J, J + g.jacobians->J.Size());
  // the factors of the mesh are computed, then destroyed
  static_cast<void>(
      mesh.GetGeometricFactors(ir, mfem::GeometricFactors::JACOBIANS));
  transformNodes(mesh, [](auto* const x, const int n, const int dim,
                          const int ps, const int cs) {
    for (int i = 0; i != n; ++i) {
      for (int c = 0; c != dim; ++c) {
        x[i * ps + c * cs] *= 2;
      }
    }
  });
  const auto* const J2 = g.jacobians->J.HostRead();
  check(std::equal(before.begin(), before.end(), J2),
        name + ": jacobians modified by the modification of the mesh");
}  // end of checkJacobiansLifetime

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
  // a NURBS mesh: the unit square of MFEM (data/square-nurbs.mesh)
  auto nurbs_stream = std::istringstream{
      "MFEM NURBS mesh v1.0\n"
      "dimension\n2\n"
      "elements\n1\n1 3 0 1 2 3\n"
      "boundary\n4\n1 1 0 1\n1 1 2 3\n1 1 3 0\n1 1 1 2\n"
      "edges\n4\n0 0 1\n0 3 2\n1 0 3\n1 1 2\n"
      "vertices\n4\n"
      "knotvectors\n2\n1 2 0 0 1 1\n1 2 0 0 1 1\n"
      "weights\n1\n1\n1\n1\n"
      "FiniteElementSpace\n"
      "FiniteElementCollection: NURBS1\n"
      "VDim: 2\n"
      "Ordering: 1\n"
      "0 0\n1 0\n1 1\n0 1\n"};
  auto nurbs = mfem::Mesh(nurbs_stream);
  refused(nurbs, "NURBS mesh");
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
    distort(vertices.HostReadWrite(), mesh.GetNV(), mesh.SpaceDimension(), 0.5,
            1, mesh.GetNV());
    mesh.SetVertices(vertices);
    compare(mesh, 1, "distorted " + name, 1);
    // unknowns of order 2: the jacobians only depend on the vertices
    compare(mesh, 2, "distorted " + name, 1);
    checkRulesAtTheSameAddress(mesh, "distorted " + name);
    if (name == "tetrahedra") {
      checkStaleShapeFunctions(mesh, "distorted " + name);
    }
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
  // curved meshes: nodes of order 2 or 3, ordered by nodes or by
  // components, and a periodic mesh, whose nodes are discontinuous
  auto curved = std::vector<std::pair<std::string, mfem::Mesh>>{};
  for (const auto& [n, t, order, ordering] :
       {std::tuple{"curved hexahedra", Type::HEXAHEDRON, 2,
                   mfem::Ordering::byNODES},
        std::tuple{"curved tetrahedra", Type::TETRAHEDRON, 3,
                   mfem::Ordering::byVDIM}}) {
    auto& mesh =
        curved.emplace_back(n, mfem::Mesh::MakeCartesian3D(2, 2, 2, t)).second;
    mesh.SetCurvature(order, false, -1, ordering);
  }
  {
    auto& mesh = curved
                     .emplace_back(
                         "curved quadrangles",
                         mfem::Mesh::MakeCartesian2D(3, 2, Type::QUADRILATERAL))
                     .second;
    mesh.SetCurvature(3, false, -1, mfem::Ordering::byVDIM);
    auto square = mfem::Mesh::MakeCartesian2D(4, 3, Type::QUADRILATERAL);
    const auto translations = std::vector<mfem::Vector>{mfem::Vector({1., 0.})};
    curved.emplace_back(
        "periodic quadrangles",
        mfem::Mesh::MakePeriodic(
            square, square.CreatePeriodicVertexMapping(translations)));
  }
  for (auto& [name, mesh] : curved) {
    check(mesh.GetNodes() != nullptr, name + ": the mesh has no nodes");
    transformNodes(mesh, [](auto* const x, const int n, const int dim,
                            const int ps, const int cs) {
      // smaller displacement than for the vertices: the nodes of order 3
      // of an element are closer
      distort(x, n, dim, 0.2, ps, cs);
    });
    compare(mesh, 1, name, 1);
    checkRulesAtTheSameAddress(mesh, name);
    // the geometric factors of MFEM assume tensor product rules for the
    // quadrangles and the hexahedra
    {
      auto fec = mfem::H1_FECollection(1, mesh.Dimension());
      auto fes = makeFiniteElementSpace(mesh, fec);
      const auto ir =
          reverse(mfem::IntRules.Get(mesh.GetElementBaseGeometry(0), 4));
      compare(mesh, *fes, ir, name + ", points in the reverse order", 1);
    }
    transformNodes(mesh, mirror);
    compare(mesh, 1, "mirrored " + name, -1);
    checkJacobiansLifetime(mesh, name);
  }
  checkRefusals();
  if (verbose) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
