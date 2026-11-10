/*!
 * \file   tests/PartialQuadratureFunctionTest.cxx
 * \brief  Tests of the `PartialQuadratureFunction` class
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#include <cmath>
#include <memory>
#include <optional>
#include <algorithm>
#include <vector>
#include <utility>
#include <cstdlib>
#include <iostream>
#ifdef DO_USE_MPI
#include <mpi.h>
#endif
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/QPEvaluator/AbstractQPEvaluator.hxx"
#include "MFEMMGIS/QPEvaluator/QPEvaluators.hxx"
#if defined(MGIS_FUNCTION_SUPPORT) && defined(MGIS_HAVE_TFEL)
#include "MGIS/Function/TFEL/Mechanics.hxx"
#endif /* defined(MGIS_FUNCTION_SUPPORT) && defined(MGIS_HAVE_TFEL) */
#include "UnitTestingUtilities.hxx"

auto parameters = mfem_mgis::unit_tests::TestParameters{};
auto ctx = mgis::Context{};

struct PartialQuadratureFunctionTest final : public tfel::tests::TestCase {
  PartialQuadratureFunctionTest()
      : tfel::tests::TestCase("MFEMMGIS", "PartialQuadratureFunctionTest") {
  }  // end of PartialQuadratureFunctionTest
  tfel::tests::TestResult execute() override {
    this->setup();
    this->test1<true>();
    this->test1<false>();
    this->test2();
    this->test3();
    // the problem is destroyed before the finalization of MFEM (this test
    // case is a static object, destroyed at the exit of the program, after
    // the release of the device)
    this->problem.reset();
    return this->result;
  }

 private:
  void setup() {
    using namespace mfem_mgis;
    auto or_raise = ctx.getThrowingFailureHandler();
    this->problem = make_unique<NonLinearEvolutionProblem>(
        ctx, ctx,
        Parameters{{"MeshFileName", parameters.mesh_file},
                   {"FiniteElementFamily", "H1"},
                   {"FiniteElementOrder", parameters.order},
                   {"UnknownsSize", 3},
                   {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
                   {"Hypothesis", "Tridimensional"},
                   {"Parallel", bool(parameters.parallel)}});
    if (isInvalid(this->problem)) {
      return;
    }
    problem->addBehaviourIntegrator(ctx, "Mechanics", 1, parameters.library,
                                    parameters.behaviour) |
        or_raise;
    auto& m1 = problem->getMaterial(ctx, 1, 0) | or_raise;
    mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
    mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
  }  // end of setup

  template <bool mutable_version>
  void test1() {
    using namespace mfem_mgis;
    using MaterialReference =
        std::conditional_t<mutable_version, Material&, const Material&>;
    using ExpectedResultType =
        std::conditional_t<mutable_version, PartialQuadratureFunction,
                           ImmutablePartialQuadratureFunctionView>;
    if (isInvalid(problem)) {
      return;
    }
    auto or_raise = ctx.getThrowingFailureHandler();
    auto& b = this->problem->getBehaviourIntegrator(ctx, 1, 0) | or_raise;
    const auto& qspace = b.getPartialQuadratureSpace();
    auto& mutable_material = b.getMaterial(ctx) | or_raise;
    auto& m1 = static_cast<MaterialReference&>(mutable_material);
    auto strain = getGradient(ctx, m1, "Strain") | or_raise;
    auto stress = getThermodynamicForce(ctx, m1, "Stress") | or_raise;
    auto p = getInternalStateVariable(ctx, m1, "EquivalentStrain") | or_raise;
    TFEL_TESTS_STATIC_ASSERT(
        (std::same_as<decltype(strain), ExpectedResultType>));
    TFEL_TESTS_STATIC_ASSERT(
        (std::same_as<decltype(stress), ExpectedResultType>));
    TFEL_TESTS_STATIC_ASSERT((std::same_as<decltype(p), ExpectedResultType>));
    TFEL_TESTS_ASSERT(strain.getNumberOfComponents() == 6);
    TFEL_TESTS_ASSERT(stress.getNumberOfComponents() == 6);
    TFEL_TESTS_ASSERT(p.getNumberOfComponents() == 1);
    // strides of the arrays of the material, which depend on their layouts
    // (6, 6 and 1 point by point without padding)
    {
      using namespace mgis::behaviour;
      // the strides of an empty view may be null
      const auto same_strides = [](const auto& f, const auto& v) {
        constexpr auto one = mgis::size_type{1};
        return (f.getDataStride() == std::max(v.getPointStride(), one)) &&
               (f.getComponentStride() ==
                std::max(v.getComponentStride(), one));
      };
      TFEL_TESTS_ASSERT(same_strides(strain, getGradientsView(m1.s1)));
      TFEL_TESTS_ASSERT(
          same_strides(stress, getThermodynamicForcesView(m1.s1)));
      TFEL_TESTS_ASSERT(same_strides(p, getInternalStateVariablesView(m1.s1)));
    }
    if (getSpaceSize(qspace) != 0) {
      if constexpr (!mutable_version) {
        TFEL_TESTS_ASSERT(strain.getValues().data() == m1.s1.gradients.data());
        TFEL_TESTS_ASSERT(stress.getValues().data() ==
                          m1.s1.thermodynamic_forces.data());
        TFEL_TESTS_ASSERT(p.getValues().data() ==
                          m1.s1.internal_state_variables.data());
      }
    }
    // requesting an unknown internal state variable
    auto oeel = getInternalStateVariable(ctx, m1, "ElasticStrain");
    TFEL_TESTS_ASSERT(isInvalid(oeel));
  }  // end of test1
  //! \brief views on external data described by `ViewSpecifications`
  void test2() {
    using namespace mfem_mgis;
    if (isInvalid(problem)) {
      return;
    }
    auto or_raise = ctx.getThrowingFailureHandler();
    auto& b = this->problem->getBehaviourIntegrator(ctx, 1, 0) | or_raise;
    const auto om = std::as_const(b).getMaterial(ctx);
    TFEL_TESTS_ASSERT(isValid(om));
    const auto qspace = om->getPartialQuadratureSpace();
    const auto n = getSpaceSize(qspace);
    auto values = std::vector<real>(3 * n, real{0});
    for (size_type i = 0; i != n; ++i) {
      values[3 * i + 1] = static_cast<real>(i);
    }
    // the second component of an array of stride 3
    auto of = PartialQuadratureFunction::borrow(
        ctx, qspace, values,
        {.data_begin = 1, .data_size = 1, .data_stride = 3});
    TFEL_TESTS_ASSERT(isValid(of));
    TFEL_TESTS_CHECK(of->getDataStride() == 3);
    TFEL_TESTS_CHECK(of->getNumberOfComponents() == 1);
    for (size_type i = 0; i != n; ++i) {
      TFEL_TESTS_CHECK(std::abs(of->getIntegrationPointValue(i) -
                                static_cast<real>(i)) < 1e-14);
    }
    // the data range is outside the stride
    auto ctx2 = Context{};
    TFEL_TESTS_CHECK(isInvalid(PartialQuadratureFunction::borrow(
        ctx2, qspace, values,
        {.data_begin = 2, .data_size = 2, .data_stride = 3})));
    // the size of the values does not match the stride
    if (n != 0) {
      auto ctx3 = Context{};
      TFEL_TESTS_CHECK(isInvalid(PartialQuadratureFunction::borrow(
          ctx3, qspace, values,
          {.data_begin = 0, .data_size = 1, .data_stride = 2})));
    }
  }  // end of test2
  /*!
   * \brief views on values stored component by component, with a padding
   * filled with NaN: the accessors, the copies and the evaluators give the
   * values of the function stored point by point
   */
  void test3() {
    using namespace mfem_mgis;
    if (isInvalid(problem)) {
      return;
    }
    auto or_raise = ctx.getThrowingFailureHandler();
    auto& b = this->problem->getBehaviourIntegrator(ctx, 1, 0) | or_raise;
    const auto om = std::as_const(b).getMaterial(ctx);
    TFEL_TESTS_ASSERT(isValid(om));
    const auto qspace = om->getPartialQuadratureSpace();
    const auto n = getSpaceSize(qspace);
    constexpr size_type nc = 6;
    const auto ld = n + 2;
    const auto value = [](const size_type i, const size_type c) {
      return real(1 + i) * real(c + 2) / 3 - real(c % 2);
    };
    auto soa = std::vector<real>(ld * nc, std::nan(""));
    auto aos = PartialQuadratureFunction(qspace, nc);
    for (size_type i = 0; i != n; ++i) {
      for (size_type c = 0; c != nc; ++c) {
        soa[c * ld + i] = value(i, c);
        aos.getIntegrationPointValues(i)[c] = value(i, c);
      }
    }
    auto of = PartialQuadratureFunction::borrow(ctx, qspace, soa,
                                                {.data_begin = 0,
                                                 .data_size = nc,
                                                 .data_stride = 1,
                                                 .data_component_stride = ld});
    TFEL_TESTS_ASSERT(isValid(of));
    TFEL_TESTS_CHECK(of->getComponentStride() == ld);
    TFEL_TESTS_CHECK(!of->hasContiguousComponents());
    for (size_type i = 0; i != n; ++i) {
      const auto v = (*of)(i);
      TFEL_TESTS_CHECK(v.size() == nc);
      for (size_type c = 0; c != nc; ++c) {
        TFEL_TESTS_CHECK(v[c] == value(i, c));
      }
    }
    // the accessors returning contiguous values refuse the view
    if (n != 0) {
      auto refused = false;
      try {
        std::ignore = of->getIntegrationPointValues(0);
      } catch (std::exception&) {
        refused = true;
      }
      TFEL_TESTS_CHECK(refused);
    }
    // a scalar view on the third component is contiguous
    auto os = PartialQuadratureFunction::borrow(ctx, qspace, soa,
                                                {.data_begin = 2 * ld,
                                                 .data_size = 1,
                                                 .data_stride = 1,
                                                 .data_component_stride = ld});
    TFEL_TESTS_ASSERT(isValid(os));
    TFEL_TESTS_CHECK(os->hasContiguousComponents());
    for (size_type i = 0; i != n; ++i) {
      TFEL_TESTS_CHECK(os->getIntegrationPointValue(i) == value(i, 2));
    }
    // copy in a function stored point by point
    auto ocopy =
        PartialQuadratureFunction::copy(ctx, std::as_const(*of).view());
    TFEL_TESTS_ASSERT(isValid(ocopy));
    TFEL_TESTS_CHECK(ocopy->getComponentStride() == 1);
    for (size_type i = 0; i != n; ++i) {
      for (size_type c = 0; c != nc; ++c) {
        TFEL_TESTS_CHECK(ocopy->getIntegrationPointValues(i)[c] == value(i, c));
      }
    }
    // assignment of values stored point by point to the view: the padding is
    // left untouched
    for (size_type i = 0; i != n; ++i) {
      for (size_type c = 0; c != nc; ++c) {
        aos.getIntegrationPointValues(i)[c] = 2 * value(i, c);
      }
    }
    TFEL_TESTS_CHECK(assign_values(ctx, of->view(), std::as_const(aos).view()));
    for (size_type c = 0; c != nc; ++c) {
      for (size_type i = 0; i != n; ++i) {
        TFEL_TESTS_CHECK(soa[c * ld + i] == 2 * value(i, c));
      }
      TFEL_TESTS_CHECK(std::isnan(soa[c * ld + n]));
      TFEL_TESTS_CHECK(std::isnan(soa[c * ld + n + 1]));
    }
    // evaluation in a given layout (external state variables with several
    // components): the values are copied only if they are not dense in it
    {
      const auto e = StandardQPEvaluator(
          qspace, nc,
          [&aos](Context&, const real, const real)
              -> std::optional<ImmutablePartialQuadratureFunctionView> {
            return std::as_const(aos).view();
          });
      for (const auto l : {mgis::FieldLayout::POINT_BY_POINT,
                           mgis::FieldLayout::COMPONENT_BY_COMPONENT}) {
        const auto r = evaluate(
            ctx, e, 0, 1, {.ensure_contiguous_storage = true, .layout = l});
        TFEL_TESTS_ASSERT(isValid(r));
        TFEL_TESTS_CHECK(r->isDense(l));
        const auto values = r->getValues();
        if (n != 0) {
          TFEL_TESTS_CHECK(
              (values.data() == std::as_const(aos).getValues().data()) ==
              (l == mgis::FieldLayout::POINT_BY_POINT));
        }
        for (size_type i = 0; i != n; ++i) {
          for (size_type c = 0; c != nc; ++c) {
            const auto p = (l == mgis::FieldLayout::POINT_BY_POINT) ? i * nc + c
                                                                    : c * n + i;
            TFEL_TESTS_CHECK(values[p] == 2 * value(i, c));
          }
        }
      }
    }
#if defined(MGIS_FUNCTION_SUPPORT) && defined(MGIS_HAVE_TFEL)
    // the evaluators give the same results to the bit in both layouts
    {
      using namespace mgis::function;
      auto r1 = PartialQuadratureFunction(qspace, 1);
      auto r2 = PartialQuadratureFunction(qspace, 1);
      TFEL_TESTS_CHECK(std::as_const(aos).view() | as_stensor<3> | vmis | r1);
      TFEL_TESTS_CHECK(std::as_const(*of).view() | as_stensor<3> | vmis | r2);
      for (size_type i = 0; i != n; ++i) {
        TFEL_TESTS_CHECK(r1.getIntegrationPointValue(i) ==
                         r2.getIntegrationPointValue(i));
      }
    }
#endif /* defined(MGIS_FUNCTION_SUPPORT) && defined(MGIS_HAVE_TFEL) */
    // the values of a component must fit in the component stride
    if (n > 1) {
      auto ctx2 = Context{};
      TFEL_TESTS_CHECK(isInvalid(
          PartialQuadratureFunction::borrow(ctx2, qspace, soa,
                                            {.data_begin = 0,
                                             .data_size = nc,
                                             .data_stride = 1,
                                             .data_component_stride = n - 1})));
    }
  }  // end of test3
  //
  std::unique_ptr<mfem_mgis::NonLinearEvolutionProblem> problem;
};

TFEL_TESTS_GENERATE_PROXY(PartialQuadratureFunctionTest,
                          "PartialQuadratureFunctionTest");

int main(int argc, char** argv) {
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(parameters, argc, argv);
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);  // mfem_mgis::getDefaultLogStream());
  m.addXMLTestOutput(
      "PartialQuadratureFunctionTest-" + std::string(parameters.behaviour) +
      "-" + std::to_string(mfem_mgis::getMPISizeForDefaultCommunicator()) +
      "-" + std::to_string(mfem_mgis::getMPIRankForDefaultCommunicator()) +
      ".xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
