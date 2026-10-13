/*!
 * \file   tests/PostProcessingsKernelsTest.cxx
 * \brief  Comparison of the computations on functions defined at the
 * integration points made by kernels (environment variable
 * `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels`) and element by element:
 * measures and integrals, values at the nodes (`updateGridFunction`), L2
 * projections, implicit gradient regularizations. The use of the kernels is
 * checked by the counter `getNumberOfPostProcessingsComputationsByKernels`.
 * \date   11/10/2026
 */

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <optional>
#include <string_view>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/L2Projection.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/PostProcessingsComputation.hxx"
#include "MFEMMGIS/Batch/IntegrationPointsFunctions.hxx"

struct TestParameters {
  //! \brief mesh
  const char* mesh_file = nullptr;
  //! \brief library of the behaviours of the tests
  const char* library = nullptr;
  //! \brief modelling hypothesis
  const char* hypothesis = "Tridimensional";
  //! \brief order of the finite elements
  int order = 2;
  //! \brief ordering of the degrees of freedom
  const char* ordering = "byNODES";
  //! \brief parallel computation
  bool parallel = false;
};

static TestParameters parseCommandLineOptions(int& argc, char** argv) {
  auto p = TestParameters{};
  mfem::OptionsParser args(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.hypothesis, "-H", "--hypothesis",
                 "Modelling hypothesis (Tridimensional or PlaneStrain).");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.ordering, "-O", "--ordering",
                 "Ordering of the degrees of freedom (byNODES or byVDIM).");
  args.AddOption(&p.parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "Perform parallel computations.");
  args.Parse();
  if ((!args.Good()) || (p.mesh_file == nullptr) || (p.library == nullptr)) {
    args.PrintUsage(std::cout);
    std::exit(EXIT_FAILURE);
  }
  return p;
}  // end of parseCommandLineOptions

//! \brief selects the computation of the post-processings
static void setPostProcessingsComputation(const char* const m) {
  ::setenv("MFEM_MGIS_POST_PROCESSINGS_COMPUTATION", m, 1);
}  // end of setPostProcessingsComputation

/*!
 * \return the largest absolute value of a vector over all the processes, or
 * an infinite value if a value is not finite
 */
template <bool parallel>
static mfem_mgis::real getMaximumAbsoluteValue(const mfem::Vector& v) {
  const auto* const pv = v.HostRead();
  auto r = mfem_mgis::real{};
  for (int i = 0; i != v.Size(); ++i) {
    if (!std::isfinite(pv[i])) {
      r = std::numeric_limits<mfem_mgis::real>::infinity();
      break;
    }
    r = std::max(r, std::abs(pv[i]));
  }
#ifdef MFEM_USE_MPI
  if constexpr (parallel) {
    MPI_Allreduce(MPI_IN_PLACE, &r, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
#endif /* MFEM_USE_MPI */
  return r;
}  // end of getMaximumAbsoluteValue

/*!
 * \return true if two vectors are equal, up to a relative tolerance of the
 * largest absolute value of the first one
 */
template <bool parallel>
static bool check(mfem_mgis::Context& ctx,
                  const std::string& what,
                  const mfem::Vector& a,
                  const mfem::Vector& b,
                  const mfem_mgis::real tolerance) {
  if (a.Size() != b.Size()) {
    return ctx.registerErrorMessage(what + ": vectors of different sizes");
  }
  auto d = mfem::Vector(a.Size());
  const auto* const pa = a.HostRead();
  const auto* const pb = b.HostRead();
  for (int i = 0; i != a.Size(); ++i) {
    d[i] = pa[i] - pb[i];
  }
  const auto e = getMaximumAbsoluteValue<parallel>(d);
  const auto s = getMaximumAbsoluteValue<parallel>(a);
  if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) {
    std::cout << what << ": largest difference " << e << ", largest value " << s
              << '\n';
  }
  if (!((s > 0) && (std::isfinite(s)) && (e <= tolerance * s))) {
    return ctx.registerErrorMessage(what + ": the results differ");
  }
  return true;
}  // end of check

//! \return true if two values are equal, up to a relative tolerance
static bool check(mfem_mgis::Context& ctx,
                  const std::string& what,
                  const mfem_mgis::real a,
                  const mfem_mgis::real b,
                  const mfem_mgis::real tolerance) {
  if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) {
    std::cout << what << ": " << a << " (element by element), " << b
              << " (kernels)\n";
  }
  if (!((a != 0) && (std::isfinite(a)) &&
        (std::abs(a - b) <= tolerance * std::abs(a)))) {
    return ctx.registerErrorMessage(what + ": the results differ");
  }
  return true;
}  // end of check

//! \brief results of the computations
template <bool parallel>
struct Results {
  //! \brief measures of the materials
  std::vector<mfem_mgis::real> measures;
  //! \brief integrals of the scalar function on the materials
  std::vector<mfem_mgis::real> integrals;
  //! \brief values at the nodes of the scalar and vector functions
  std::vector<mfem::Vector> nodes;
  //! \brief L2 projections of the scalar and vector functions
  std::vector<mfem::Vector> projections;
  //! \brief implicit gradient regularizations of the scalar and vector
  //! functions
  std::vector<mfem::Vector> regularizations;
};

/*!
 * \brief computes the measures, integrals, values at the nodes and L2
 * projections of a scalar function and of a function with three components
 * on all the materials
 */
template <bool parallel>
static std::optional<Results<parallel>> compute(
    mfem_mgis::Context& ctx,
    mfem_mgis::NonLinearEvolutionProblem& p,
    const std::vector<int>& materials,
    const std::vector<mfem_mgis::ImmutablePartialQuadratureFunctionView>&
        scalar_functions,
    const std::vector<mfem_mgis::ImmutablePartialQuadratureFunctionView>&
        vector_functions) {
  using namespace mfem_mgis;
  auto r = Results<parallel>{};
  for (size_type i = 0; i != materials.size(); ++i) {
    const auto obi = p.getBehaviourIntegrator(ctx, materials[i], 0);
    if (isInvalid(obi)) {
      return {};
    }
    r.measures.push_back(computeMeasure(*obi));
    r.integrals.push_back(computeIntegral<real>(*obi, scalar_functions[i]));
  }
  for (const auto& fcts : {scalar_functions, vector_functions}) {
    auto gf = makeGridFunction<parallel>(ctx, fcts);
    if (gf == nullptr) {
      return {};
    }
    if (!updateGridFunction<parallel>(ctx, *gf, fcts)) {
      return {};
    }
    r.nodes.emplace_back(*gf);
    auto& fespace =
        p.getFiniteElementDiscretization().getFiniteElementSpace<parallel>();
    auto& f = LinearSolverFactory<parallel>::getFactory();
    auto s = f.generate(ctx, "CGSolver", fespace,
                        {{"VerbosityLevel", 0},
                         {"AbsoluteTolerance", 1e-30},
                         {"RelativeTolerance", 1e-14},
                         {"MaximumNumberOfIterations", 10000}});
    if (isInvalid(s)) {
      return {};
    }
    const auto oprojection = computeL2Projection<parallel>(ctx, s, fcts);
    if (isInvalid(oprojection)) {
      return {};
    }
    r.projections.emplace_back(*(oprojection->result));
    const auto oregularization =
        computeImplicitGradientRegularization<parallel>(ctx, s, fcts, 0.02);
    if (isInvalid(oregularization)) {
      return {};
    }
    r.regularizations.emplace_back(*(oregularization->result));
  }
  return r;
}  // end of compute

template <bool parallel>
static bool test(mfem_mgis::Context& ctx, const TestParameters& params) {
  using namespace mfem_mgis;
  const auto dim = std::string_view{params.hypothesis} == "PlaneStrain" ? 2 : 3;
  auto oproblem = construct<NonLinearEvolutionProblem>(
      ctx, Parameters{{"MeshFileName", params.mesh_file},
                      {"FiniteElementFamily", "H1"},
                      {"FiniteElementOrder", params.order},
                      {"UnknownsSize", dim},
                      {"FiniteElementSpaceOrdering", params.ordering},
                      {"Hypothesis", params.hypothesis},
                      {"Parallel", parallel}});
  if (isInvalid(oproblem)) {
    return false;
  }
  auto& p = *oproblem;
  // a behaviour integrator on each material, and functions defined at its
  // integration points
  const auto& attributes =
      getMaterialsAttributes(p.getFiniteElementDiscretization());
  auto materials = std::vector<int>{};
  auto scalar_functions =
      std::vector<std::shared_ptr<PartialQuadratureFunction>>{};
  auto vector_functions =
      std::vector<std::shared_ptr<PartialQuadratureFunction>>{};
  for (int i = 0; i != attributes.Size(); ++i) {
    const auto m = attributes[i];
    if (!p.addBehaviourIntegrator(ctx, "Mechanics", m, params.library,
                                  "IsotropicLinearElasticity")) {
      return false;
    }
    const auto omaterial = p.getMaterial(ctx, m, 0);
    if (isInvalid(omaterial)) {
      return false;
    }
    const auto qs = omaterial->getPartialQuadratureSpacePointer();
    // the functions differ from one material to another
    const auto c = static_cast<real>(m);
    auto evaluate = [&qs](const auto f) {
      return PartialQuadratureFunction::evaluate(
          qs, [f](const mfem::FiniteElement&, mfem::ElementTransformation& tr) {
            auto x = mfem::Vector{};
            tr.Transform(tr.GetIntPoint(), x);
            return f(x[0], x[1], x.Size() > 2 ? x[2] : real{0});
          });
    };
    scalar_functions.push_back(
        evaluate([c](const real x, const real y, const real z) {
          return c + std::cos(x) * y + z * z;
        }));
    auto fct = std::make_shared<PartialQuadratureFunction>(qs, 3);
    const auto f0 = evaluate([c](const real x, const real y, const real z) {
      return c * x * y - z;
    });
    const auto f2 = evaluate([](const real x, const real y, const real z) {
      return std::exp(-x) * y * y + std::sin(z);
    });
    for (size_type o = 0; o != qs->getNumberOfIntegrationPoints(); ++o) {
      auto values = fct->getIntegrationPointValues(o);
      values[0] = f0->getIntegrationPointValue(o);
      values[1] = scalar_functions.back()->getIntegrationPointValue(o);
      values[2] = f2->getIntegrationPointValue(o);
    }
    vector_functions.push_back(fct);
    materials.push_back(m);
  }
  auto views = [](const auto& fcts) {
    auto r = std::vector<ImmutablePartialQuadratureFunctionView>{};
    for (const auto& f : fcts) {
      r.push_back(f->view());
    }
    return r;
  };
  // element by element, then by kernels
  setPostProcessingsComputation("ElementByElement");
  const auto n0 = getNumberOfPostProcessingsComputationsByKernels();
  const auto r1 = compute<parallel>(ctx, p, materials, views(scalar_functions),
                                    views(vector_functions));
  if (!r1.has_value()) {
    return false;
  }
  if (getNumberOfPostProcessingsComputationsByKernels() != n0) {
    return ctx.registerErrorMessage(
        "computations made by kernels with ElementByElement");
  }
  // the computations by kernels must not give nodes to the mesh (as the
  // geometric factors of MFEM would)
  const auto& mesh = p.getFiniteElementDiscretization().getMesh<parallel>();
  const auto has_nodes = mesh.GetNodes() != nullptr;
  setPostProcessingsComputation("Kernels");
  const auto r2 = compute<parallel>(ctx, p, materials, views(scalar_functions),
                                    views(vector_functions));
  if (!r2.has_value()) {
    return false;
  }
  if ((mesh.GetNodes() != nullptr) != has_nodes) {
    return ctx.registerErrorMessage(
        "the computations by kernels gave nodes to the mesh");
  }
  // two measures and integrals per material, values at the nodes (if the
  // elements are small enough for the local projections by kernels), L2
  // projections and implicit gradient regularizations of two functions
  const auto nk = getNumberOfPostProcessingsComputationsByKernels() - n0;
  const auto nd = p.getFiniteElementDiscretization()
                      .getFiniteElementSpace<parallel>()
                      .GetTypicalFE()
                      ->GetDof();
  const auto values_at_nodes_by_kernels =
      nd <= batch::maximumNumberOfNodesOfLocalProjections;
  const auto expected =
      2 * materials.size() + (values_at_nodes_by_kernels ? 6 : 4);
  if (nk != expected) {
    return ctx.registerErrorMessage(
        "unexpected number of computations by kernels (" + std::to_string(nk) +
        " instead of " + std::to_string(expected) + ")");
  }
  auto success = true;
  for (size_type i = 0; i != materials.size(); ++i) {
    const auto m = std::to_string(materials[i]);
    success = check(ctx, "measure of material " + m, r1->measures[i],
                    r2->measures[i], 1e-12) &&
              success;
    success = check(ctx, "integral on material " + m, r1->integrals[i],
                    r2->integrals[i], 1e-12) &&
              success;
  }
  for (size_type i = 0; i != 2; ++i) {
    const auto f = std::string{i == 0 ? "scalar" : "vector"} + " function";
    success = check<parallel>(ctx, "values at the nodes of the " + f,
                              r1->nodes[i], r2->nodes[i], 1e-12) &&
              success;
    success = check<parallel>(ctx, "L2 projection of the " + f,
                              r1->projections[i], r2->projections[i], 1e-10) &&
              success;
    success = check<parallel>(
                  ctx, "implicit gradient regularization of the " + f,
                  r1->regularizations[i], r2->regularizations[i], 1e-10) &&
              success;
  }
  return success;
}  // end of test

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const auto params = parseCommandLineOptions(argc, argv);
  auto ctx = Context{};
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
  if (!success) {
    std::cerr << ctx.getErrorMessage() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}  // end of main
