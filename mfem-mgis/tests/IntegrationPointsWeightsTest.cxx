/*!
 * \file   IntegrationPointsWeightsTest.cxx
 * \brief  Test of the weights of the integration points kept by the non
 * linear evolution problems (`getIntegrationPointsWeights`) and of the mean
 * values of the thermodynamic forces, which use them
 * (`computeMeanThermodynamicForcesValues`).
 *
 * The unit cube made of two materials (1 for x < 1/2, 2 for x > 1/2) is
 * refined once, then its vertices are moved by x -> x (1 + x) / 2 and
 * y -> y (3 + y) / 4: the elements are boxes of different sizes, the volumes
 * of the materials are 3/8 and 5/8, and the integrals of x over them are
 * 9/128 and 55/128. For each material:
 * - the sum of the weights is the volume of the material;
 * - each weight is the one computed from the transformation of its element
 *   (`AbstractBehaviourIntegrator::getIntegrationPointWeight`), at the
 *   offset of its point in the partial quadrature space;
 * - a second call returns the same array, with the same values;
 * - the integrals of thermodynamic forces which are uniform, then which
 *   depend on x, and the volumes, are the analytical ones: the mean values
 *   follow the states, the weights being the kept ones.
 *
 * The weights are not computed again by the following calls, which is only
 * visible if the mesh is modified (the documentation of
 * `getIntegrationPointsWeights` forbids it): the vertices are moved again
 * at the end of the test, and the weights and the volumes are still the
 * ones of the former mesh.
 *
 * With the `--parallel` option, each process checks its own weights, and
 * the sums are the ones of all the processes.
 * \date   14/11/2026
 */

#include <cmath>
#include <array>
#include <string>
#include <vector>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/vector.hpp"
#include "mfem/mesh/mesh.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "Tolerance.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \return the sum of the values of all the processes
 * \tparam parallel: if true, the computations are parallel
 * \param[in] v: value of this process
 */
template <bool parallel>
static double sum(const double v) {
  if constexpr (parallel) {
#ifdef MFEM_USE_MPI
    auto r = double{};
    MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return r;
#else  /* MFEM_USE_MPI */
    return v;
#endif /* MFEM_USE_MPI */
  } else {
    return v;
  }
}  // end of sum

/*!
 * \brief moves the vertices of a mesh of the unit cube made of two
 * materials: the elements are boxes of different sizes, the interface
 * between the materials being moved to x = 3/8
 * \param[in, out] m: mesh
 */
static void moveVertices(mfem::Mesh& m) {
  m.Transform([](const mfem::Vector& x, mfem::Vector& y) {
    y = x;
    y[0] = x[0] * (1 + x[0]) / 2;
    y[1] = x[1] * (3 + x[1]) / 4;
  });
}  // end of moveVertices

//! \return true if two values are equal, up to a relative tolerance
static bool areClose(const double a, const double b, const double eps) {
  return std::abs(a - b) <= eps * std::abs(b);
}  // end of areClose

/*!
 * \brief weights of the integration points and mean values of the
 * thermodynamic forces of the two materials
 * \tparam parallel: if true, the computations are parallel
 * \param[in, out] ctx: execution context
 * \param[in] mesh: mesh of the unit cube made of two materials
 * \param[in] library: behaviour library
 */
template <bool parallel>
static void test(mfem_mgis::Context& ctx,
                 const char* const mesh,
                 const char* const library) {
  using namespace mfem_mgis;
  auto or_die = ctx.getFatalFailureHandler();
  auto problem = construct<NonLinearEvolutionProblem>(
                     ctx, Parameters{{"MeshFileName", mesh},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", 1},
                                     {"UnknownsSize", 3},
                                     {"NumberOfUniformRefinements", 1},
                                     {"Hypothesis", "Tridimensional"},
                                     {"Parallel", parallel}}) |
                 or_die;
  auto fed = problem.getFiniteElementDiscretization();
  moveVertices(fed.template getMesh<parallel>());
  constexpr auto volumes = std::array<double, 2u>{3. / 8, 5. / 8};
  constexpr auto x_integrals = std::array<double, 2u>{9. / 128, 55. / 128};
  // sums of the values at the points: in float, the roundings of the
  // weights and of the sums made by the library
  constexpr auto eps = double{mfem_mgis::unit_tests::tolerance(1e-12, 1e-4)};
  for (const size_type m : {1, 2}) {
    problem.addBehaviourIntegrator(ctx, "Mechanics", m, library,
                                   "IsotropicLinearHardeningPlasticity") |
        or_die;
  }
  auto& impl = problem.template getImplementation<parallel>();
  const auto& fes = fed.template getFiniteElementSpace<parallel>();
  // first coordinate of the integration points of each material
  auto coordinates = std::array<std::vector<real>, 2u>{};
  // addresses of the weights of each material
  auto addresses = std::array<const std::vector<real>*, 2u>{};
  for (const size_type m : {1, 2}) {
    const auto im = static_cast<std::size_t>(m - 1);
    const auto name = "material " + std::to_string(m);
    const auto& bi = problem.getBehaviourIntegrator(ctx, m, 0) | or_die;
    const auto& qspace = bi.getPartialQuadratureSpace();
    const auto& w = impl.getIntegrationPointsWeights(ctx, m, 0) | or_die;
    addresses[im] = &w;
    const auto np =
        static_cast<std::size_t>(qspace.getNumberOfIntegrationPoints());
    check(w.size() == np, name + ": number of weights");
    coordinates[im].assign(np, real{0});
    // weights computed from the transformations of the elements
    auto same = w.size() == np;
    auto visited = std::size_t{};
    auto position = mfem::Vector(3);
    for (int i = 0; same && (i != fes.GetNE()); ++i) {
      if (fes.GetAttribute(i) != m) {
        continue;
      }
      const auto& e = *(fes.GetFE(i));
      auto& tr = *(fes.GetElementTransformation(i));
      const auto& ir = bi.getIntegrationRule(e, tr);
      const auto offset = qspace.getOffset(i);
      for (int j = 0; j != ir.GetNPoints(); ++j) {
        const auto o = static_cast<std::size_t>(offset + j);
        if (o >= np) {
          same = false;
          break;
        }
        const auto& ip = ir.IntPoint(j);
        tr.SetIntPoint(&ip);
        const auto expected = bi.getIntegrationPointWeight(tr, ip);
        same = same && areClose(w[o], expected,
                                mfem_mgis::unit_tests::tolerance(1e-14));
        tr.Transform(ip, position);
        coordinates[im][o] = position[0];
        ++visited;
      }
    }
    check(same && (visited == np),
          name + ": weights computed from the transformations of the elements");
    // volume of the material
    auto volume = 0.;
    for (const auto v : w) {
      volume += v;
    }
    volume = sum<parallel>(volume);
    std::cout << name << ": sum of the weights " << volume << " (volume "
              << volumes[im] << ")\n";
    check(areClose(volume, volumes[im], eps), name + ": sum of the weights");
    // second call: the same array, with the same values
    const auto copy = w;
    const auto* const values = w.data();
    const auto& w2 = impl.getIntegrationPointsWeights(ctx, m, 0) | or_die;
    check((&w2 == &w) && (w2.data() == values) && (w2 == copy),
          name + ": weights computed again by a second call");
  }
  // integrals of the thermodynamic forces and volumes of the materials: the
  // first component of the thermodynamic forces is `a + b x`, the other ones
  // are uniform
  const auto check_integrals = [&](const double a, const double b,
                                   const std::string& state) {
    for (const size_type m : {1, 2}) {
      const auto im = static_cast<std::size_t>(m - 1);
      auto& mat = problem.getMaterial(ctx, m, 0) | or_die;
      // the thermodynamic forces are written through their view, in all the
      // layouts
      mat.synchronizeOnHost();
      const auto sv = mgis::behaviour::getThermodynamicForcesView(mat.s1);
      const auto nc =
          static_cast<size_type>(mat.s1.thermodynamic_forces_stride);
      for (size_type o = 0; o != static_cast<size_type>(mat.n); ++o) {
        const auto x = coordinates[im][static_cast<std::size_t>(o)];
        sv(o, 0) = static_cast<real>(m * (a + b * x));
        for (size_type k = 1; k != nc; ++k) {
          sv(o, k) = static_cast<real>(m * a * (k + 1));
        }
      }
    }
    const auto ores = computeMeanThermodynamicForcesValues<parallel>(ctx, impl);
    if (isInvalid(ores)) {
      std::cerr << ctx.getErrorMessage() << '\n';
      check(false, state + ": mean values of the thermodynamic forces");
      return;
    }
    const auto& [integrals, material_volumes] = *ores;
    for (const size_type m : {1, 2}) {
      const auto im = static_cast<std::size_t>(m - 1);
      const auto name = state + ", material " + std::to_string(m);
      const auto mi = static_cast<std::size_t>(m);
      if ((integrals.size() <= mi) || (material_volumes.size() <= mi)) {
        check(false, name + ": no integral");
        continue;
      }
      check(areClose(sum<parallel>(material_volumes[mi]), volumes[im], eps),
            name + ": volume");
      auto ok = true;
      for (std::size_t k = 0; k != integrals[mi].size(); ++k) {
        const auto expected =
            k == 0 ? m * (a * volumes[im] + b * x_integrals[im])
                   : m * a * static_cast<double>(k + 1) * volumes[im];
        ok = ok && areClose(sum<parallel>(integrals[mi][k]), expected, eps);
      }
      check(ok && (!integrals[mi].empty()),
            name + ": integrals of the thermodynamic forces");
      // the weights are the kept ones
      const auto& w = impl.getIntegrationPointsWeights(ctx, m, 0) | or_die;
      check(&w == addresses[im], name + ": weights computed again");
    }
  };
  check_integrals(1e6, 0, "uniform thermodynamic forces");
  check_integrals(-3e6, 5e6, "thermodynamic forces depending on x");
  // the weights are kept: they are not computed again, even once the mesh
  // is modified, which must not be done after the first call (see
  // `getIntegrationPointsWeights`). The volumes of the materials would
  // become 33/128 and 95/128.
  auto kept = std::array<std::vector<real>, 2u>{};
  for (const auto im : {std::size_t{0}, std::size_t{1}}) {
    kept[im] = *(addresses[im]);
  }
  moveVertices(fed.template getMesh<parallel>());
  for (const size_type m : {1, 2}) {
    const auto im = static_cast<std::size_t>(m - 1);
    const auto& w = impl.getIntegrationPointsWeights(ctx, m, 0) | or_die;
    check((&w == addresses[im]) && (w == kept[im]),
          "material " + std::to_string(m) +
              ": weights computed again after a modification of the mesh");
  }
  check_integrals(2e6, 0, "kept weights, uniform thermodynamic forces");
}  // end of test

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* mesh = "cube-two-materials.mesh";
  const char* library = "src/libBehaviourTest.so";
  auto parallel = false;
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&mesh, "-m", "--mesh",
                 "mesh of the unit cube made of two materials (1 for "
                 "x < 1/2, 2 for x > 1/2)");
  args.AddOption(&library, "-l", "--library", "behaviour library");
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "parallel computations");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  auto ctx = Context{};
  if (parallel) {
#ifdef MFEM_USE_MPI
    test<true>(ctx, mesh, library);
#else  /* MFEM_USE_MPI */
    check(false, "parallel computations require MPI");
#endif /* MFEM_USE_MPI */
  } else {
    test<false>(ctx, mesh, library);
  }
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
