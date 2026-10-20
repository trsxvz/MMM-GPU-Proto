/*!
 * \file   test-sycl.cxx
 * \brief
 * \author Thomas Helfer
 * \date 19/05/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of theses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <array>
#include <vector>
#include <memory>
#include <iostream>
#include <sycl/sycl.hpp>

#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"

#include "MFront/GenericMaterialProperty/MaterialProperty.h"
#include "PoissonRatioTest-generic-parallel-sycl.hxx"
#include "Inconel600_YoungModulus-generic-parallel-sycl.hxx"
#include "Inconel600_YoungModulus_qt-generic-parallel-sycl.hxx"
#include "UO2_YoungModulus_Martin1989-generic-parallel-sycl.hxx"
#include "LinearDataInterpolationTest-generic-parallel-sycl.hxx"
#include "MaterialLawCallTest-generic-parallel-sycl.hxx"
#include "symbols.hxx"

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
    TFEL_TESTS_ASSERT(checkGenericParallelSymbols("sycl"));
    return this->result;
  }  // end of execute

 private:
  //! \brief type of the memory allocator
  using usm_allocator_type =
      sycl::usm_allocator<double, sycl::usm::alloc::host>;
  //! \brief a simple alias for a vectors
  using usm_vector = std::vector<double, usm_allocator_type>;

  [[nodiscard]] auto to_device(const usm_vector &v, ::sycl::queue &queue) {
    auto deleter = [queue](double *const p) { sycl::free(p, queue); };
    auto ptr = std::unique_ptr<double, decltype(deleter)>(
        sycl::malloc_device<double>(v.size(), queue), deleter);
    queue.copy(v.data(), ptr.get(), v.size()).wait();
    return std::move(ptr);
  }

  void to_host(usm_vector &v, const auto &ptr, ::sycl::queue &queue) {
    queue.copy(ptr.get(), v.data(), v.size()).wait();
  }

  [[nodiscard]] static double fct(const double TK) noexcept {
    const double TC = TK - 273.15;
    return (-3.1636e-3 * TC * TC - 3.8654 * TC + 2.1421e+4) * 1e7;
  }  // end of fct

  void test1() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto usm_allocator = usm_allocator_type{Q};
    //
    auto E = usm_vector(4, Q);
    auto E_qt = usm_vector(4, Q);
    const auto T = usm_vector({300, 500, 300, 800}, Q);
    //
    auto d_E = to_device(E, Q);
    auto d_E_qt = to_device(E_qt, Q);
    auto d_T = to_device(T, Q);
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{d_T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, d_E.get(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, &options);
    Inconel600_YoungModulus_qt(&output_qt, d_E_qt.get(), 1, args.data(),
                               args_strides.data(), 1, 4, policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    //
    to_host(E, d_E, Q);
    to_host(E_qt, d_E_qt, Q);
    //
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i] - fct(T[i])) < eps);
    }
  }

  void test2() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto usm_allocator = usm_allocator_type{Q};
    //
    auto E = usm_vector(4, Q);
    auto E_qt = usm_vector(4, Q);
    const auto T = usm_vector({300}, Q);
    //
    auto d_E = to_device(E, Q);
    auto d_E_qt = to_device(E_qt, Q);
    auto d_T = to_device(T, Q);
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{d_T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{0};
    Inconel600_YoungModulus(&output, d_E.get(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, &options);
    Inconel600_YoungModulus_qt(&output_qt, d_E_qt.get(), 1, args.data(),
                               args_strides.data(), 1, 4, policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    //
    to_host(E, d_E, Q);
    to_host(E_qt, d_E_qt, Q);
    //
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[0])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i] - fct(T[0])) < eps);
    }
  }

  void test3() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto usm_allocator = usm_allocator_type{Q};
    //
    auto E = usm_vector(1, Q);
    auto E_qt = usm_vector(1, Q);
    const auto T = usm_vector({300}, Q);
    //
    auto d_E = to_device(E, Q);
    auto d_E_qt = to_device(E_qt, Q);
    auto d_T = to_device(T, Q);
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{d_T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{0};
    Inconel600_YoungModulus(&output, d_E.get(), 0, args.data(),
                            args_strides.data(), 1, 1, policy, &options);
    Inconel600_YoungModulus_qt(&output_qt, d_E_qt.get(), 0, args.data(),
                               args_strides.data(), 1, 1, policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    //
    to_host(E, d_E, Q);
    to_host(E_qt, d_E_qt, Q);
    //
    TFEL_TESTS_ASSERT(std::abs(E[0] - fct(T[0])) < eps);
    TFEL_TESTS_ASSERT(std::abs(E_qt[0] - fct(T[0])) < eps);
  }

  void test4() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto usm_allocator = usm_allocator_type{Q};
    //
    auto E = usm_vector(4, Q);
    auto E_qt = usm_vector(4, Q);
    const auto T = usm_vector({300, 500, 300, 800}, Q);
    //
    auto d_E = to_device(E, Q);
    auto d_E_qt = to_device(E_qt, Q);
    auto d_T = to_device(T, Q);
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{d_T.get()};
    Inconel600_YoungModulus2(&output, d_E.get(), args.data(), 1, 4, policy,
                             &options);
    Inconel600_YoungModulus_qt2(&output_qt, d_E_qt.get(), args.data(), 1, 4,
                                policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    //
    to_host(E, d_E, Q);
    to_host(E_qt, d_E_qt, Q);
    //
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
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto usm_allocator = usm_allocator_type{Q};
    //
    auto E = usm_vector(4, Q);
    auto E_qt = usm_vector(4, Q);
    const auto T = usm_vector({-300, 500, 300, 800}, Q);
    //
    auto d_E = to_device(E, Q);
    auto d_E_qt = to_device(E_qt, Q);
    auto d_T = to_device(T, Q);
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{d_T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, d_E.get(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, &options);
    Inconel600_YoungModulus_qt(&output_qt, d_E_qt.get(), 1, args.data(),
                               args_strides.data(), 1, 4, policy, &options);
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
    //
    to_host(E, d_E, Q);
    to_host(E_qt, d_E_qt, Q);
    //
    for (std::size_t i = 1; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i] - fct(T[i])) < eps);
    }
  }

  void test6() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto usm_allocator = usm_allocator_type{Q};
    //
    auto E = usm_vector(4, Q);
    auto E_qt = usm_vector(4, Q);
    const auto T = usm_vector({-300, 500, 300, 800}, Q);
    //
    auto d_E = to_device(E, Q);
    auto d_E_qt = to_device(E_qt, Q);
    auto d_T = to_device(T, Q);
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{d_T.get()};
    Inconel600_YoungModulus2(&output, d_E.get(), args.data(), 1, 4, policy,
                             &options);
    Inconel600_YoungModulus_qt2(&output_qt, d_E_qt.get(), args.data(), 1, 4,
                                policy, &options);
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
    //
    to_host(E, d_E, Q);
    to_host(E_qt, d_E_qt, Q);
    //
    for (std::size_t i = 1; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i] - fct(T[i])) < eps);
    }
  }

  void test7() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto usm_allocator = usm_allocator_type{Q};
    //
    auto E = usm_vector(4 * 3, Q);
    auto E_qt = usm_vector(4 * 3, Q);
    const auto T = usm_vector({300, -300, 500, -500, 300, -300, 800, -800}, Q);
    //
    auto d_E = to_device(E, Q);
    auto d_E_qt = to_device(E_qt, Q);
    auto d_T = to_device(T, Q);
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{d_T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{2};
    Inconel600_YoungModulus(&output, d_E.get(), 3, args.data(),
                            args_strides.data(), 1, 4, policy, &options);
    Inconel600_YoungModulus_qt(&output_qt, d_E_qt.get(), 3, args.data(),
                               args_strides.data(), 1, 4, policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    //
    to_host(E, d_E, Q);
    to_host(E_qt, d_E_qt, Q);
    //
    for (std::size_t i = 0; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i * 3] - fct(T[i * 2])) < eps);
      TFEL_TESTS_ASSERT(std::abs(E_qt[i * 3] - fct(T[i * 2])) < eps);
    }
  }

  void test8() {
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto usm_allocator = usm_allocator_type{Q};
    //
    auto E = usm_vector(1, Q);
    const auto T = usm_vector({300}, Q);
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    //
    auto d_E = to_device(E, Q);
    auto d_T = to_device(T, Q);
    //
    const auto args = std::array<const double *, 1u>{d_T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, d_E.get(), 1, args.data(),
                            args_strides.data(), 3, 1, policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, -5);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "invalid number of arguments (3 given, 1 expected)");
    //
    Inconel600_YoungModulus(&output, d_E.get(), 0, args.data(),
                            args_strides.data(), 1, 1, policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, -7);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "if the output is uniform, all the arguments must "
                           "be uniform (i.e. their strides must be null)");
    //
    const auto args_strides2 = std::array<mfront_gmp_size_type, 1u>{0};
    Inconel600_YoungModulus(&output, d_E.get(), 0, args.data(),
                            args_strides2.data(), 1, 10, policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, -7);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "if the output is uniform, the number of points "
                           "must be equal to one (10 given)");
  }

  void test9() {
    constexpr auto eps = double{1e-14};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto usm_allocator = usm_allocator_type{Q};
    //
    auto nu = usm_vector(4, Q);
    auto d_nu = to_device(nu, Q);
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    PoissonRatioTest(&output, d_nu.get(), 1, nullptr, nullptr, 0, 4, policy,
                     &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    to_host(nu, d_nu, Q);
    for (std::size_t i = 0; i != nu.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(nu[i] - 0.39991) < eps);
    }
    nu.resize(1);
    nu.at(0) = 0;
    auto d_nu2 = to_device(nu, Q);
    PoissonRatioTest(&output, d_nu2.get(), 0, nullptr, nullptr, 0, 1, policy,
                     &options);
    to_host(nu, d_nu2, Q);
    TFEL_TESTS_ASSERT(std::abs(nu[0] - 0.39991) < eps);
  }

  void test10() {
    constexpr auto eps = double{1e-14};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto usm_allocator = usm_allocator_type{Q};
    //
    auto nu = usm_vector(4, Q);
    auto d_nu = to_device(nu, Q);
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    PoissonRatioTest2(&output, d_nu.get(), nullptr, 0, 4, policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    to_host(nu, d_nu, Q);
    for (std::size_t i = 0; i != nu.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(nu[i] - 0.39991) < eps);
    }
    nu.resize(1);
    nu.at(0) = 0;
    auto d_nu2 = to_device(nu, Q);
    PoissonRatioTest2(&output, d_nu2.get(), nullptr, 0, 1, policy, &options);
    to_host(nu, d_nu2, Q);
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
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    //
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto E = usm_vector(4, Q);
    const auto T = usm_vector({300, 500, 300, 800}, Q);
    const auto f = usm_vector({0, 0.1, 0.15, 0.12}, Q);
    //
    auto d_E = to_device(E, Q);
    auto d_T = to_device(T, Q);
    auto d_f = to_device(f, Q);
    //
    const auto args = std::array<const double *, 2u>{d_T.get(), d_f.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 2u>{1, 1};
    UO2_YoungModulus_Martin1989(&output, d_E.get(), 1, args.data(),
                                args_strides.data(), args.size(), 4, policy,
                                &options);
    to_host(E, d_E, Q);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - young(T[i], f[i])) < eps);
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
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    auto Q = sycl::queue{sycl::gpu_selector_v};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
    auto v = usm_vector(T_values.size(), Q);
    const auto T = usm_vector(T_values.begin(), T_values.end(), Q);
    auto d_v = to_device(v, Q);
    auto d_T = to_device(T, Q);
    const auto args = std::array<const double *, 1u>{d_T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    mp(&output, d_v.get(), 1, args.data(), args_strides.data(), args.size(),
       T.size(), policy, &options);
    to_host(v, d_v, Q);
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
  auto &m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput(std::string{argv[1]} + ".xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
