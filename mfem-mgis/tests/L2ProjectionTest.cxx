/*!
 * \file   tests/L2ProjectionTest.cxx
 * \brief  Tests of the L2 projection
 * \author Thomas Helfer
 * \date   20/01/2026
 */

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/L2Projection.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "UnitTestingUtilities.hxx"
#include "Tolerance.hxx"

struct TestParameters {
  const char* mesh_file = nullptr;
  int linearsolver = 0;
  int order = 1;
  int parallel = 0;
};  // end of struct TestParameters

static void parseCommandLineOptions(TestParameters& params,
                                    int argc,
                                    char** argv) {
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&params.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&params.linearsolver, "-ls", "--linearsolver",
                 "identifier of the linear solver: 0 -> CG, 1 -> GMRES, 2 -> "
                 "UMFPack (serial), 3-> MUMPS(serial), 2 -> HypreFGMRES "
                 "(//), 3 -> HyprePCG (//), 4 -> HypreGMRES (//)");
  args.AddOption(&params.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&params.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (params.mesh_file == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
}  // end of parseCommandLineOptions

// single component test
template <bool parallel>
bool test(mfem_mgis::Context& ctx, const TestParameters& params) {
  using namespace mfem_mgis;
  auto fed = FiniteElementDiscretization{
      ctx,
      {{"MeshFileName", params.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", params.order},
       {"UnknownsSize", 1},
       {"NumberOfUniformRefinements", parallel ? 1 : 0},
       {"Parallel", parallel}}};
  auto space = std::make_shared<PartialQuadratureSpace>(
      fed, MaterialIdentifier{.id = 5},
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation& tr) noexcept
      -> const mfem::IntegrationRule& {
        const auto order = 2 * tr.OrderGrad(&e);
        return mfem::IntRules.Get(e.GetGeomType(), order);
      });
  auto fct = PartialQuadratureFunction::evaluate(
      space, [](const real x, const real y) noexcept { return cos(x) * y; });
  auto s = mfem_mgis::unit_tests::getLinearSolver<parallel>(
      ctx, fed.getFiniteElementSpace<parallel>(), params);
  if (isInvalid(s)) {
    return false;
  }
  const auto oresult = computeL2Projection<parallel>(ctx, s, {*fct});
  if (isInvalid(oresult)) {
    return false;
  }
  mfem::ParaViewDataCollection exporter("Result-test1");
  if constexpr (parallel) {
#ifdef MFEM_USE_MPI
    exporter.SetMesh(&(fed.getMesh<true>()));
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  } else {
    exporter.SetMesh(&(fed.getMesh<false>()));
  }
  exporter.SetDataFormat(mfem::VTKFormat::BINARY);
  exporter.RegisterField("Result", oresult->result.get());
  exporter.SetCycle(1);
  exporter.SetTime(1);
  exporter.Save();
  return true;
}  // end of test

// multi-component test
template <bool parallel>
bool test2(mfem_mgis::Context& ctx, const TestParameters& params) {
  using namespace mfem_mgis;
  auto fed = FiniteElementDiscretization{
      ctx,
      {{"MeshFileName", params.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", params.order},
       {"UnknownsSize", 1},
       {"NumberOfUniformRefinements", parallel ? 1 : 0},
       {"Parallel", parallel}}};
  auto space = std::make_shared<PartialQuadratureSpace>(
      fed, MaterialIdentifier{.id = 5},
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation& tr) noexcept
      -> const mfem::IntegrationRule& {
        const auto order = 2 * tr.OrderGrad(&e);
        return mfem::IntRules.Get(e.GetGeomType(), order);
      });
  auto f1 = PartialQuadratureFunction::evaluate(
      space, [](const real x, const real y) noexcept { return cos(x) * y; });
  auto f2 = PartialQuadratureFunction::evaluate(
      space,
      [](const real x, const real y) noexcept { return exp(-x) * y * y; });
  auto fct = std::make_shared<PartialQuadratureFunction>(space, 2);
  auto f1_values = f1->getValues();
  auto f2_values = f2->getValues();
  auto fct_values = fct->getValues();
  for (size_type idx = 0; idx != space->getNumberOfIntegrationPoints(); ++idx) {
    fct_values[2 * idx] = f1_values[idx];
    fct_values[2 * idx + 1] = f2_values[idx];
  }
  auto s = mfem_mgis::unit_tests::getLinearSolver<parallel>(
      ctx, fed.getFiniteElementSpace<parallel>(), params);
  if (isInvalid(s)) {
    return false;
  }
  const auto oresult = computeL2Projection<parallel>(ctx, s, {*fct});
  if (isInvalid(oresult)) {
    return false;
  }
  mfem::ParaViewDataCollection exporter("Result-test2");
  if constexpr (parallel) {
#ifdef MFEM_USE_MPI
    exporter.SetMesh(&(fed.getMesh<true>()));
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  } else {
    exporter.SetMesh(&(fed.getMesh<false>()));
  }
  exporter.SetDataFormat(mfem::VTKFormat::BINARY);
  exporter.RegisterField("Result", oresult->result.get());
  exporter.SetCycle(1);
  exporter.SetTime(1);
  exporter.Save();
  return true;
}  // end of test2

// single component test on a boundary
template <bool parallel>
bool test3(mfem_mgis::Context& ctx, const TestParameters& params) {
  using namespace mfem_mgis;
  auto fed = FiniteElementDiscretization{
      ctx,
      {{"MeshFileName", params.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", params.order},
       {"UnknownsSize", 1},
       {"NumberOfUniformRefinements", parallel ? 1 : 0},
       {"Parallel", parallel}}};
  auto space = std::make_shared<PartialQuadratureSpace>(
      fed, BoundaryIdentifier{.id = 6},
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation& tr) noexcept
      -> const mfem::IntegrationRule& {
        const auto order = 2 * tr.OrderGrad(&e);
        return mfem::IntRules.Get(e.GetGeomType(), order);
      });
  auto fct = PartialQuadratureFunction::evaluate(
      space, [](const real x, const real y) noexcept { return cos(x) * y; });
  auto s = mfem_mgis::unit_tests::getLinearSolver<parallel>(
      ctx, fed.getFiniteElementSpace<parallel>(), params);
  if (isInvalid(s)) {
    return false;
  }
  const auto oresult = computeL2Projection<parallel>(ctx, s, {*fct});
  if (isInvalid(oresult)) {
    return false;
  }
  // the projection is done on the submesh of the boundary
  if (oresult->submesh.get() == nullptr) {
    return ctx.registerErrorMessage("no submesh created");
  }
  mfem::ParaViewDataCollection exporter("Result-test3");
  exporter.SetMesh(oresult->submesh.get());
  exporter.SetDataFormat(mfem::VTKFormat::BINARY);
  exporter.RegisterField("Result", oresult->result.get());
  exporter.SetCycle(1);
  exporter.SetTime(1);
  exporter.Save();
  return true;
}  // end of test3

/*!
 * \brief projection of linear functions with one and three components: the
 * L2 projection of a function of the finite element space is this function,
 * whose values at the nodes are the ones of its interpolant. The values are
 * compared on the host, whatever the device.
 */
template <bool parallel>
bool test4(mfem_mgis::Context& ctx, const TestParameters& params) {
  using namespace mfem_mgis;
  auto fed = FiniteElementDiscretization{
      ctx,
      {{"MeshFileName", params.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", params.order},
       {"UnknownsSize", 1},
       {"NumberOfUniformRefinements", parallel ? 1 : 0},
       {"Parallel", parallel}}};
  auto space = std::make_shared<PartialQuadratureSpace>(
      fed, MaterialIdentifier{.id = 5},
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation& tr) noexcept
          -> const mfem::IntegrationRule& {
        const auto order = 2 * tr.OrderGrad(&e);
        return mfem::IntRules.Get(e.GetGeomType(), order);
      });
  // coefficients a, b and c of the components a + b x + c y
  const auto coefficients = std::array<std::array<real, 3>, 3>{
      {{1, 2, -1}, {-0.5, 0, 3}, {2, -1, 0.5}}};
  for (const auto n : {1, 3}) {
    auto fct = std::make_shared<PartialQuadratureFunction>(space, n);
    auto values = fct->getValues();
    for (int c = 0; c != n; ++c) {
      const auto [a, b, d] = coefficients[c];
      const auto f = PartialQuadratureFunction::evaluate(
          space, [a, b, d](const real x, const real y) noexcept {
            return a + b * x + d * y;
          });
      const auto f_values = f->getValues();
      for (size_type idx = 0; idx != space->getNumberOfIntegrationPoints();
           ++idx) {
        values[n * idx + c] = f_values[idx];
      }
    }
    auto s = mfem_mgis::unit_tests::getLinearSolver<parallel>(
        ctx, fed.getFiniteElementSpace<parallel>(), params);
    if (isInvalid(s)) {
      return false;
    }
    const auto oresult = computeL2Projection<parallel>(ctx, s, {*fct});
    if (isInvalid(oresult)) {
      return false;
    }
    auto& result = *(oresult->result);
    auto reference = [&result] {
      if constexpr (parallel) {
        return GridFunction<true>(result.ParFESpace());
      } else {
        return GridFunction<false>(result.FESpace());
      }
    }();
    auto coefficient = mfem::VectorFunctionCoefficient(
        n, [&coefficients, n](const mfem::Vector& x, mfem::Vector& v) {
          for (int c = 0; c != n; ++c) {
            const auto [a, b, d] = coefficients[c];
            v[c] = a + b * x[0] + d * x[1];
          }
        });
    reference.ProjectCoefficient(coefficient);
    const auto* const presult = result.HostRead();
    const auto* const preference = reference.HostRead();
    // largest difference and largest value
    auto e = std::array<real, 2>{};
    for (int i = 0; i != result.Size(); ++i) {
      e[0] = std::max(e[0], std::abs(presult[i] - preference[i]));
      e[1] = std::max(e[1], std::abs(preference[i]));
    }
#ifdef MFEM_USE_MPI
    if constexpr (parallel) {
      MPI_Allreduce(MPI_IN_PLACE, e.data(), 2, mfem_mgis::mpi_type<real>,
                    MPI_MAX, MPI_COMM_WORLD);
    }
#endif /* MFEM_USE_MPI */
    // in float, the largest relative difference measured is 6.2e-4: the
    // conjugate gradient stops at a relative residual of 1.2e-5 (bound of
    // mfem-mgis), multiplied by the conditioning of the mass matrix
    if (!(e[0] <= mfem_mgis::unit_tests::tolerance(1e-8, 5e-3) * e[1])) {
      auto msg = std::ostringstream{};
      msg << "test4: invalid projection of a linear function with " << n
          << " component(s) (largest difference " << e[0] << ", largest value "
          << e[1] << ")";
      return ctx.registerErrorMessage(msg.str());
    }
  }
  return true;
}  // end of test4

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  // options treatment
  auto params = TestParameters{};
  initialize(argc, argv);
  parseCommandLineOptions(params, argc, argv);
  const auto success = [&ctx, &params] {
    if (params.parallel) {
#ifdef MFEM_USE_MPI
      return test<true>(ctx, params);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return test<false>(ctx, params);
  }();
  const auto success2 = [&ctx, &params] {
    if (params.parallel) {
#ifdef MFEM_USE_MPI
      return test2<true>(ctx, params);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return test2<false>(ctx, params);
  }();
  const auto success3 = [&ctx, &params] {
    if (params.parallel) {
#ifdef MFEM_USE_MPI
      return test3<true>(ctx, params);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return test3<false>(ctx, params);
  }();
  const auto success4 = [&ctx, &params] {
    if (params.parallel) {
#ifdef MFEM_USE_MPI
      return test4<true>(ctx, params);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return test4<false>(ctx, params);
  }();
  if (!(success && success2 && success3 && success4)) {
    std::cerr << ctx.getErrorMessage() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
