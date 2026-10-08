/*!
 * \file   test-cxx.cxx
 * \brief
 * \author Thomas Helfer
 * \date   22/04/2026
 */

#include <array>
#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <type_traits>
#ifdef MFRONT_GPMP_TEST_USE_KOKKOS
#include <Kokkos_Core.hpp>
#ifdef KOKKOS_ENABLE_CUDA
#include <cuda_runtime.h>
#endif /* KOKKOS_ENABLE_CUDA */
#endif /* MFRONT_GPMP_TEST_USE_KOKKOS */

#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"

#include "MFront/GenericMaterialProperty/MaterialProperty.h"

// This test is shared by the backends based on the C++ programming model
// (`stlpar`, `openmp`, `kokkos`), selected by the `MFRONT_GPMP_TEST_BACKEND`
// macro (`stlpar` by default). The arrays are allocated by `std::vector`
// (on the host, or in managed memory with `nvc++ -gpu=mem:managed`), or in
// the shared memory space of Kokkos if the `MFRONT_GPMP_TEST_USE_KOKKOS`
// macro is defined.
#ifndef MFRONT_GPMP_TEST_BACKEND
#define MFRONT_GPMP_TEST_BACKEND stlpar
#endif /* MFRONT_GPMP_TEST_BACKEND */
#define MFRONT_GPMP_TEST_STRINGIFY2(x) #x
#define MFRONT_GPMP_TEST_STRINGIFY(x) MFRONT_GPMP_TEST_STRINGIFY2(x)
// clang-format off: the spaces would appear in the name of the header
#define MFRONT_GPMP_TEST_HEADER(name) \
  MFRONT_GPMP_TEST_STRINGIFY(name-generic-parallel-MFRONT_GPMP_TEST_BACKEND.hxx)
// clang-format on

#include MFRONT_GPMP_TEST_HEADER(PoissonRatioTest)
#include MFRONT_GPMP_TEST_HEADER(Inconel600_YoungModulus)
#include MFRONT_GPMP_TEST_HEADER(Inconel600_YoungModulus_qt)
#include MFRONT_GPMP_TEST_HEADER(UO2_YoungModulus_Martin1989)
#include MFRONT_GPMP_TEST_HEADER(LinearDataInterpolationTest)
#include MFRONT_GPMP_TEST_HEADER(MaterialLawCallTest)
#include "symbols.hxx"

#ifdef MFRONT_GPMP_TEST_USE_KOKKOS
// execution space of the material properties (`execution_space` option of
// the backend)
#ifndef MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE
#define MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE Kokkos::DefaultExecutionSpace
#endif /* MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE */
//! \brief allocator in the shared memory space of Kokkos
template <typename T>
struct TestAllocator {
  using value_type = T;
  TestAllocator() = default;
  template <typename T2>
  TestAllocator(const TestAllocator<T2>&) noexcept {}
  T* allocate(const std::size_t n) {
    return static_cast<T*>(
        Kokkos::kokkos_malloc<Kokkos::SharedSpace>(n * sizeof(T)));
  }
  void deallocate(T* const p, const std::size_t) noexcept {
    Kokkos::kokkos_free<Kokkos::SharedSpace>(p);
  }
  template <typename T2>
  bool operator==(const TestAllocator<T2>&) const noexcept {
    return true;
  }
};
#else  /* MFRONT_GPMP_TEST_USE_KOKKOS */
template <typename T>
using TestAllocator = std::allocator<T>;
#endif /* MFRONT_GPMP_TEST_USE_KOKKOS */

//! \brief array of values given to the material properties
template <typename T>
using test_vector = std::vector<T, TestAllocator<T>>;

struct Iconel600YoungModulusTest : public tfel::tests::TestCase {
  Iconel600YoungModulusTest()
      : tfel::tests::TestCase("MFront/MaterialProperty/GenericParallel",
                              "Iconel600YoungModulusTest") {
  }  // end of Iconel600YoungModulusTest

  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    this->test4();
    this->test5();
    this->test6();
    this->test7();
    this->test8();
    this->test9();
    this->test10();
    this->test11();
    this->test12();
    this->test13();
    this->test14();
    this->test15();
    this->test16();
    this->test17();
    TFEL_TESTS_ASSERT(checkGenericParallelSymbols(
        MFRONT_GPMP_TEST_STRINGIFY(MFRONT_GPMP_TEST_BACKEND)));
    return this->result;
  }  // end of execute

 private:
  [[nodiscard]] static double fct(const double TK) noexcept {
    const double TC = TK - 273.15;
    return (-3.1636e-3 * TC * TC - 3.8654 * TC + 2.1421e+4) * 1e7;
  }  // end of fct

  void test1() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    auto E_qt = test_vector<double>(4);
    const auto T = test_vector<double>{300, 500, 300, 800};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, E.data(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt(&output_qt, E_qt.data(), 1, args.data(),
                               args_strides.data(), 1, 4, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i] - fct(T[i])) < eps);
    }
  }

  void test2() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    auto E_qt = test_vector<double>(4);
    const auto T = test_vector<double>{300};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{0};
    Inconel600_YoungModulus(&output, E.data(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt(&output_qt, E_qt.data(), 1, args.data(),
                               args_strides.data(), 1, 4, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[0])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i] - fct(T[0])) < eps);
    }
  }

  void test3() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(1);
    auto E_qt = test_vector<double>(1);
    const auto T = test_vector<double>{300};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{0};
    Inconel600_YoungModulus(&output, E.data(), 0, args.data(),
                            args_strides.data(), 1, 1, policy, nullptr);
    Inconel600_YoungModulus_qt(&output_qt, E_qt.data(), 0, args.data(),
                               args_strides.data(), 1, 1, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    TFEL_TESTS_ASSERT(std::abs(E[0] - fct(T[0])) < eps);
    TFEL_TESTS_ASSERT(std::abs(E_qt[0] - fct(T[0])) < eps);
  }

  void test4() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    auto E_qt = test_vector<double>(4);
    const auto T = test_vector<double>{300, 500, 300, 800};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.data()};
    Inconel600_YoungModulus2(&output, E.data(), args.data(), 1, 4, policy,
                             nullptr);
    Inconel600_YoungModulus_qt2(&output_qt, E_qt.data(), args.data(), 1, 4,
                                policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i] - fct(T[i])) < eps);
    }
  }
  // check detection of physical bounds violation
  void test5() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    auto E_qt = test_vector<double>(4);
    const auto T = test_vector<double>{-300, 500, 300, 800};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, E.data(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt(&output_qt, E_qt.data(), 1, args.data(),
                               args_strides.data(), 1, 4, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output.first_failure_index, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, -1);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output_qt.first_failure_index, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, -1);
    TFEL_TESTS_ASSERT(std::abs(E[0]) < eps);
    TFEL_TESTS_ASSERT(std::abs(E_qt[0]) < eps);
    for (std::size_t i = 1; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i] - fct(T[i])) < eps);
    }
  }

  void test6() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    auto E_qt = test_vector<double>(4);
    const auto T = test_vector<double>{-300, 500, 300, 800};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.data()};
    Inconel600_YoungModulus2(&output, E.data(), args.data(), 1, 4, policy,
                             nullptr);
    Inconel600_YoungModulus_qt2(&output_qt, E_qt.data(), args.data(), 1, 4,
                                policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output.first_failure_index, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, -1);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output_qt.first_failure_index, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, -1);
    TFEL_TESTS_ASSERT(std::abs(E[0]) < eps);
    TFEL_TESTS_ASSERT(std::abs(E_qt[0]) < eps);
    for (std::size_t i = 1; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i] - fct(T[i])) < eps);
    }
  }

  void test7() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4 * 3);
    auto E_qt = test_vector<double>(4 * 3);
    const auto T =
        test_vector<double>{300, -300, 500, -500, 300, -300, 800, -800};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{2};
    Inconel600_YoungModulus(&output, E.data(), 3, args.data(),
                            args_strides.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt(&output_qt, E_qt.data(), 3, args.data(),
                               args_strides.data(), 1, 4, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    for (std::size_t i = 0; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i * 3] - fct(T[i * 2])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i * 3] - fct(T[i * 2])) < eps);
    }
  }

  void test8() {
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(1);
    const auto T = test_vector<double>{300};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, E.data(), 1, args.data(),
                            args_strides.data(), 3, 1, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -5);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "invalid number of arguments (3 given, 1 expected)");
    //
    Inconel600_YoungModulus(&output, E.data(), 0, args.data(),
                            args_strides.data(), 1, 1, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -7);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "if the output is uniform, all the arguments must "
                           "be uniform (i.e. their strides must be null)");
    //
    const auto args_strides2 = std::array<mfront_gmp_size_type, 1u>{0};
    Inconel600_YoungModulus(&output, E.data(), 0, args.data(),
                            args_strides2.data(), 1, 10, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -7);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "if the output is uniform, the number of points "
                           "must be equal to one (10 given)");
  }

  void test9() {
    constexpr auto eps = double{1e-14};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto nu = test_vector<double>(4);
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    PoissonRatioTest(&output, nu.data(), 1, nullptr, nullptr, 0, 4, policy,
                     nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != nu.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(nu[i] - 0.39991) < eps);
    }
    nu.resize(1);
    nu.at(0) = 0;
    PoissonRatioTest(&output, nu.data(), 0, nullptr, nullptr, 0, 1, policy,
                     nullptr);
    TFEL_TESTS_ASSERT(std::abs(nu[0] - 0.39991) < eps);
  }

  void test10() {
    constexpr auto eps = double{1e-14};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto nu = test_vector<double>(4);
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    PoissonRatioTest2(&output, nu.data(), nullptr, 0, 4, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != nu.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(nu[i] - 0.39991) < eps);
    }
    nu.resize(1);
    nu.at(0) = 0;
    PoissonRatioTest2(&output, nu.data(), nullptr, 0, 1, policy, nullptr);
    TFEL_TESTS_ASSERT(std::abs(nu[0] - 0.39991) < eps);
  }

  void test11() {
    constexpr auto eps = double{1e-2};
    //
    auto young = [](const double T, const double f) {
      constexpr auto E0 = 2.2693e11;
      constexpr auto dE_dT = -1.53994698e7;
      constexpr auto d2E_dT2 = -1.9198278e4;
      constexpr auto f0 = 0.4;
      return (1 - f / f0) * (E0 + dE_dT * T + (d2E_dT2 / 2) * T * T);
    };
    //
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    const auto T = test_vector<double>{300, 500, 300, 800};
    const auto f = test_vector<double>{0, 0.1, 0.15, 0.12};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 2u>{T.data(), f.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 2u>{1, 1};
    UO2_YoungModulus_Martin1989(&output, E.data(), 1, args.data(),
                                args_strides.data(), args.size(), 4, policy,
                                nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - young(T[i], f[i])) < eps);
    }
  }

  void test15() {
    // execution options given explicitly: instance of the execution space
    // with Kokkos, no queue otherwise
    constexpr auto eps = double{1e-2};
#ifdef MFRONT_GPMP_TEST_USE_KOKKOS
    auto space = MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE{};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &space};
#else  /* MFRONT_GPMP_TEST_USE_KOKKOS */
    const auto options = mfront_gpmp_ExecutionOptions{};
#endif /* MFRONT_GPMP_TEST_USE_KOKKOS */
    auto output = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    const auto T = test_vector<double>{300, 500, 300, 800};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, E.data(), 1, args.data(),
                            args_strides.data(), 1, 4,
                            mfront_gmp_OutOfBoundsPolicy{}, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
    }
  }

  void test16() {
    // asynchronous evaluation: the statistics are given by the execution
    // options and reported once the computations are completed. The second
    // temperature is out of its physical bounds.
    constexpr auto eps = double{1e-2};
    auto statistics = test_vector<mfront_gpmp_Statistics>(1);
    mfront_gpmp_initializeStatistics(statistics.data());
#ifdef MFRONT_GPMP_TEST_USE_KOKKOS
    auto space = MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE{};
    const auto options = mfront_gpmp_ExecutionOptions{
        .queue = &space, .statistics = statistics.data()};
#else  /* MFRONT_GPMP_TEST_USE_KOKKOS */
    const auto options =
        mfront_gpmp_ExecutionOptions{.statistics = statistics.data()};
#endif /* MFRONT_GPMP_TEST_USE_KOKKOS */
    auto output = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    const auto T = test_vector<double>{300, -1, 500, 800};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    Inconel600_YoungModulus(&output, E.data(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
#ifdef MFRONT_GPMP_TEST_USE_KOKKOS
    space.fence();
#endif /* MFRONT_GPMP_TEST_USE_KOKKOS */
    Inconel600_YoungModulus_generic_parallel_report(&output, statistics.data(),
                                                    policy);
    TFEL_TESTS_CHECK_EQUAL(output.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output.first_failure_index, mfront_gmp_size_type{1});
    for (const auto i : {0, 2, 3}) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
    }
  }

  void test17() {
    // stack size of the threads of the device preserved by the execution
    // options. With a CUDA execution space of Kokkos, an error is reported
    // if it is too small for the kernel. On the host, the option has no
    // effect.
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(2);
    const auto T = test_vector<double>{300, 500};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto preserve = mfront_gpmp_ExecutionOptions{
        .queue = nullptr, .statistics = nullptr, .preserve_stack_size = 1};
#if defined(MFRONT_GPMP_TEST_USE_KOKKOS) && defined(KOKKOS_ENABLE_CUDA)
    if constexpr (std::is_same_v<MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE,
                                 Kokkos::Cuda>) {
      auto initial_stack_size = std::size_t{};
      TFEL_TESTS_ASSERT(cudaDeviceGetLimit(&initial_stack_size,
                                           cudaLimitStackSize) == cudaSuccess);
      // the kernels require at least 1024 bytes (see
      // mfront::generic_parallel::gpu::ensureStackSize)
      TFEL_TESTS_ASSERT(cudaDeviceSetLimit(cudaLimitStackSize, 256) ==
                        cudaSuccess);
      Inconel600_YoungModulus(&output, E.data(), 1, args.data(),
                              args_strides.data(), 1, 2, policy, &preserve);
      TFEL_TESTS_CHECK_EQUAL(output.status, -6);
      TFEL_TESTS_ASSERT(std::string(output.msg).starts_with(
          "the stack size of the threads of the device must be at least"));
      TFEL_TESTS_ASSERT(cudaDeviceSetLimit(cudaLimitStackSize,
                                           initial_stack_size) == cudaSuccess);
    }
#endif /* defined(MFRONT_GPMP_TEST_USE_KOKKOS) && \
          defined(KOKKOS_ENABLE_CUDA) */
    // the stack size has been increased by the previous evaluations
    Inconel600_YoungModulus(&output, E.data(), 1, args.data(),
                            args_strides.data(), 1, 2, policy, &preserve);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
    }
  }

  // bound violated with the warning policy: the result is computed at all
  // points and must be used with caution (status 1), as with the `generic`
  // interface
  void test14() {
    constexpr auto eps = double{1e-2};
    auto young = [](const double T, const double f) {
      constexpr auto E0 = 2.2693e11;
      constexpr auto dE_dT = -1.53994698e7;
      constexpr auto d2E_dT2 = -1.9198278e4;
      constexpr auto f0 = 0.4;
      return (1 - f / f0) * (E0 + dE_dT * T + (d2E_dT2 / 2) * T * T);
    };
    auto output = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    const auto T = test_vector<double>{300, 3000, 300, 800};
    const auto f = test_vector<double>{0, 0.1, 0.15, 0.12};
    const auto policy = GENERIC_MATERIALPROPERTY_WARNING_POLICY;
    const auto args = std::array<const double *, 2u>{T.data(), f.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 2u>{1, 1};
    UO2_YoungModulus_Martin1989(&output, E.data(), 1, args.data(),
                                args_strides.data(), args.size(), 4, policy,
                                nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 1);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 1);
    // only the second temperature is out of its bounds
    TFEL_TESTS_CHECK_EQUAL(output.number_of_unreliable_results,
                           mfront_gmp_size_type{1});
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - young(T[i], f[i])) <
                        eps * std::abs(young(T[i], f[i])));
    }
  }

  void test12() {
    // linear interpolation of data
    this->checkTemperatureDependentMaterialProperty(
        LinearDataInterpolationTest, {293.15, 493.15, 793.15, 893.15},
        {240e9, 210e9, 175e9, 170e9});
  }

  void test13() {
    // material property calling another material property
    const auto T = std::array<double, 4u>{300, 500, 300, 800};
    auto G = std::array<double, 4u>{};
    for (std::size_t i = 0; i != T.size(); ++i) {
      G[i] = fct(T[i]) / (2 * (1 + 0.3));
    }
    this->checkTemperatureDependentMaterialProperty(MaterialLawCallTest, T, G);
  }

  /*!
   * \brief evaluate a material property depending on the temperature at four
   * points and compare the results to reference values
   * \param[in] mp: material property
   * \param[in] T_values: temperatures
   * \param[in] references: reference values
   */
  void checkTemperatureDependentMaterialProperty(
      const auto mp,
      const std::array<double, 4u>& T_values,
      const std::array<double, 4u>& references) {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto v = test_vector<double>(T_values.size());
    const auto T = test_vector<double>(T_values.begin(), T_values.end());
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    mp(&output, v.data(), 1, args.data(), args_strides.data(), args.size(),
       T.size(), policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != v.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(v[i] - references[i]) < eps);
    }
  }
};

TFEL_TESTS_GENERATE_PROXY(Iconel600YoungModulusTest,
                          "Iconel600YoungModulusTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main(const int argc, const char *const *argv) {
  if (argc != 2) {
    std::cerr << "no test name given on the command line (argc: " << argc
              << ")";
    return EXIT_FAILURE;
  }
#ifdef MFRONT_GPMP_TEST_USE_KOKKOS
  auto kokkos = Kokkos::ScopeGuard{};
#endif /* MFRONT_GPMP_TEST_USE_KOKKOS */
  auto &m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput(std::string{argv[1]} + ".xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
