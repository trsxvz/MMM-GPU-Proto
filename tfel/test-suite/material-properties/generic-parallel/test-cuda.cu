/*!
 * \file   test-cuda.cu
 * \brief
 * \author Thomas Helfer
 * \date   22/04/2026
 */

#include <array>
#include <vector>
#include <iostream>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"

#include "MFront/GenericMaterialProperty/MaterialProperty.h"
#include "PoissonRatioTest-generic-parallel-cuda.hxx"
#include "Inconel600_YoungModulus-generic-parallel-cuda.hxx"
#include "Inconel600_YoungModulus_qt-generic-parallel-cuda.hxx"
#include "UO2_YoungModulus_Martin1989-generic-parallel-cuda.hxx"
#include "LinearDataInterpolationTest-generic-parallel-cuda.hxx"
#include "MaterialLawCallTest-generic-parallel-cuda.hxx"
#include "symbols.hxx"

/*!
 * \brief kernel writing the given values after a delay, so that the
 * computations submitted after it on the same stream must wait for it
 * \param[out] T: written values
 * \param[in] values: values
 * \param[in] cycles: number of clock cycles of the delay
 */
__global__ void writeAfterDelay(double *const T,
                                const double *const values,
                                const long long cycles) {
  const auto start = clock64();
  while (clock64() - start < cycles) {
  }
  T[threadIdx.x] = values[threadIdx.x];
}  // end of writeAfterDelay

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
    TFEL_TESTS_ASSERT(checkGenericParallelSymbols("cuda"));
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
    auto E = thrust::host_vector<double>(4);
    const auto T = thrust::host_vector<double>{300, 500, 300, 800};
    //
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_E(E.size());
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, thrust::raw_pointer_cast(d_E.data()), 1,
                            args.data(), args_strides.data(), 1, 4, policy,
                            nullptr);
    E = d_E;
    //
    auto E_qt = thrust::host_vector<double>(4);
    thrust::device_vector<double> d_E_qt(E_qt.size());
    Inconel600_YoungModulus_qt(
        &output_qt, thrust::raw_pointer_cast(d_E_qt.data()), 1, args.data(),
        args_strides.data(), 1, 4, policy, nullptr);
    E_qt = d_E_qt;
    //
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
    auto E = thrust::host_vector<double>(4);
    auto E_qt = thrust::host_vector<double>(4);
    const auto T = thrust::host_vector<double>{300};
    //
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_E(E.size());
    thrust::device_vector<double> d_E_qt(E_qt.size());
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{0};
    //
    Inconel600_YoungModulus(&output, thrust::raw_pointer_cast(d_E.data()), 1,
                            args.data(), args_strides.data(), 1, 4, policy,
                            nullptr);
    Inconel600_YoungModulus_qt(
        &output_qt, thrust::raw_pointer_cast(d_E_qt.data()), 1, args.data(),
        args_strides.data(), 1, 4, policy, nullptr);
    E = d_E;
    E_qt = d_E_qt;
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
    auto E = thrust::host_vector<double>(1);
    auto E_qt = thrust::host_vector<double>(1);
    const auto T = thrust::host_vector<double>{300};
    //
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_E(E.size());
    thrust::device_vector<double> d_E_qt(E_qt.size());
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{0};
    Inconel600_YoungModulus(&output, thrust::raw_pointer_cast(d_E.data()), 0,
                            args.data(), args_strides.data(), 1, 1, policy,
                            nullptr);
    Inconel600_YoungModulus_qt(
        &output_qt, thrust::raw_pointer_cast(d_E_qt.data()), 0, args.data(),
        args_strides.data(), 1, 1, policy, nullptr);
    E = d_E;
    E_qt = d_E_qt;
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
    auto E = thrust::host_vector<double>(4);
    auto E_qt = thrust::host_vector<double>(4);
    const auto T = thrust::host_vector<double>{300, 500, 300, 800};
    //
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_E(E.size());
    thrust::device_vector<double> d_E_qt(E_qt.size());
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    Inconel600_YoungModulus2(&output, thrust::raw_pointer_cast(d_E.data()),
                             args.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt2(&output_qt,
                                thrust::raw_pointer_cast(d_E_qt.data()),
                                args.data(), 1, 4, policy, nullptr);
    E = d_E;
    E_qt = d_E_qt;
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
    auto E = thrust::host_vector<double>(4);
    auto E_qt = thrust::host_vector<double>(4);
    const auto T = thrust::host_vector<double>{-300, 500, 300, 800};
    //
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_E(E.size());
    thrust::device_vector<double> d_E_qt(E_qt.size());
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, thrust::raw_pointer_cast(d_E.data()), 1,
                            args.data(), args_strides.data(), 1, 4, policy,
                            nullptr);
    Inconel600_YoungModulus_qt(
        &output_qt, thrust::raw_pointer_cast(d_E_qt.data()), 1, args.data(),
        args_strides.data(), 1, 4, policy, nullptr);
    E = d_E;
    E_qt = d_E_qt;
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
    auto E = thrust::host_vector<double>(4);
    auto E_qt = thrust::host_vector<double>(4);
    const auto T = thrust::host_vector<double>{-300, 500, 300, 800};
    //
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_E(E.size());
    thrust::device_vector<double> d_E_qt(E_qt.size());
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    Inconel600_YoungModulus2(&output, thrust::raw_pointer_cast(d_E.data()),
                             args.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt2(&output_qt,
                                thrust::raw_pointer_cast(d_E_qt.data()),
                                args.data(), 1, 4, policy, nullptr);
    E = d_E;
    E_qt = d_E_qt;
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
    auto E = thrust::host_vector<double>(4 * 3);
    auto E_qt = thrust::host_vector<double>(4 * 3);
    const auto T =
        thrust::host_vector<double>{300, -300, 500, -500, 300, -300, 800, -800};
    //
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_E(E.size());
    thrust::device_vector<double> d_E_qt(E_qt.size());
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{2};
    Inconel600_YoungModulus(&output, thrust::raw_pointer_cast(d_E.data()), 3,
                            args.data(), args_strides.data(), 1, 4, policy,
                            nullptr);
    Inconel600_YoungModulus_qt(
        &output_qt, thrust::raw_pointer_cast(d_E_qt.data()), 3, args.data(),
        args_strides.data(), 1, 4, policy, nullptr);
    E = d_E;
    E_qt = d_E_qt;
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
    auto E = thrust::host_vector<double>(1);
    const auto T = thrust::host_vector<double>{300};
    //
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_E(E.size());
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, thrust::raw_pointer_cast(d_E.data()), 1,
                            args.data(), args_strides.data(), 3, 1, policy,
                            nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -5);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "invalid number of arguments (3 given, 1 expected)");
    //
    Inconel600_YoungModulus(&output, thrust::raw_pointer_cast(d_E.data()), 0,
                            args.data(), args_strides.data(), 1, 1, policy,
                            nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -7);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "if the output is uniform, all the arguments must "
                           "be uniform (i.e. their strides must be null)");
    //
    const auto args_strides2 = std::array<mfront_gmp_size_type, 1u>{0};
    Inconel600_YoungModulus(&output, thrust::raw_pointer_cast(d_E.data()), 0,
                            args.data(), args_strides2.data(), 1, 10, policy,
                            nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -7);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "if the output is uniform, the number of points "
                           "must be equal to one (10 given)");
  }

  void test9() {
    constexpr auto eps = double{1e-14};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto nu = thrust::host_vector<double>(4);
    //
    thrust::device_vector<double> d_nu(nu.size());
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    PoissonRatioTest(&output, thrust::raw_pointer_cast(d_nu.data()), 1, nullptr,
                     nullptr, 0, 4, policy, nullptr);
    nu = d_nu;
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != nu.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(nu[i] - 0.39991) < eps);
    }
    nu = thrust::host_vector<double>(1);
    d_nu = thrust::device_vector<double>(nu.size());
    PoissonRatioTest(&output, thrust::raw_pointer_cast(d_nu.data()), 0, nullptr,
                     nullptr, 0, 1, policy, nullptr);
    nu = d_nu;
    TFEL_TESTS_ASSERT(std::abs(nu[0] - 0.39991) < eps);
  }

  void test10() {
    constexpr auto eps = double{1e-14};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto nu = thrust::host_vector<double>(4);
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    //
    thrust::device_vector<double> d_nu(nu.size());
    //
    PoissonRatioTest2(&output, thrust::raw_pointer_cast(d_nu.data()), nullptr,
                      0, 4, policy, nullptr);
    nu = d_nu;
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != nu.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(nu[i] - 0.39991) < eps);
    }
    nu = thrust::host_vector<double>(1);
    d_nu = thrust::device_vector<double>(nu.size());
    PoissonRatioTest2(&output, thrust::raw_pointer_cast(d_nu.data()), nullptr,
                      0, 1, policy, nullptr);
    nu = d_nu;
    TFEL_TESTS_ASSERT(std::abs(nu[0] - 0.39991) < eps);
  }
  void test11() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    //
    auto young = [](const double T, const double f) {
      constexpr auto E0 = 2.2693e11;
      constexpr auto dE_dT = -1.53994698e7;
      constexpr auto d2E_dT2 = -1.9198278e4;
      constexpr auto f0 = 0.4;
      return (1 - f / f0) * (E0 + dE_dT * T + (d2E_dT2 / 2) * T * T);
    };
    //
    auto E = thrust::host_vector<double>(4);
    const auto T = thrust::host_vector<double>{300, 500, 300, 800};
    const auto f = thrust::host_vector<double>{0, 0.1, 0.15, 0.12};
    //
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_f(f);
    thrust::device_vector<double> d_E(E.size());
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args =
        std::array<const double *, 2u>{thrust::raw_pointer_cast(d_T.data()),
                                       thrust::raw_pointer_cast(d_f.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 2u>{1, 1};
    UO2_YoungModulus_Martin1989(&output, thrust::raw_pointer_cast(d_E.data()),
                                1, args.data(), args_strides.data(),
                                args.size(), 4, policy, nullptr);
    E = d_E;
    //
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - young(T[i], f[i])) < eps);
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
    auto E = thrust::host_vector<double>(4);
    const auto T = thrust::host_vector<double>{300, 3000, 300, 800};
    const auto f = thrust::host_vector<double>{0, 0.1, 0.15, 0.12};
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_f(f);
    thrust::device_vector<double> d_E(E.size());
    const auto policy = GENERIC_MATERIALPROPERTY_WARNING_POLICY;
    const auto args =
        std::array<const double *, 2u>{thrust::raw_pointer_cast(d_T.data()),
                                       thrust::raw_pointer_cast(d_f.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 2u>{1, 1};
    UO2_YoungModulus_Martin1989(&output, thrust::raw_pointer_cast(d_E.data()),
                                1, args.data(), args_strides.data(),
                                args.size(), 4, policy, nullptr);
    E = d_E;
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

  void test15() {
    // computations on a stream given by the execution options, which is not
    // the default stream
    constexpr auto eps = double{1e-2};
    auto stream = cudaStream_t{};
    TFEL_TESTS_ASSERT(cudaStreamCreateWithFlags(
                          &stream, cudaStreamNonBlocking) == cudaSuccess);
    const auto T = thrust::host_vector<double>{300, 500, 300, 800};
    const thrust::device_vector<double> d_values(T);
    thrust::device_vector<double> d_T(T.size(), 0.);
    thrust::device_vector<double> d_E(T.size());
    // the temperatures are written on the stream, after a delay: the
    // material property must be evaluated after this kernel, without
    // synchronization by the host
    writeAfterDelay<<<1, 4, 0, stream>>>(
        thrust::raw_pointer_cast(d_T.data()),
        thrust::raw_pointer_cast(d_values.data()), 100'000'000);
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = stream};
    auto output = mfront_gpmp_OutputStatus{};
    Inconel600_YoungModulus(&output, thrust::raw_pointer_cast(d_E.data()), 1,
                            args.data(), args_strides.data(), 1, 4,
                            mfront_gmp_OutOfBoundsPolicy{}, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    const auto E = thrust::host_vector<double>(d_E);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
    }
    TFEL_TESTS_ASSERT(cudaStreamDestroy(stream) == cudaSuccess);
  }

  void test16() {
    // asynchronous evaluation on a stream: the statistics, in managed memory,
    // are given by the execution options and reported once the stream is
    // synchronized. The second temperature is out of its physical bounds.
    constexpr auto eps = double{1e-2};
    auto stream = cudaStream_t{};
    TFEL_TESTS_ASSERT(cudaStreamCreateWithFlags(
                          &stream, cudaStreamNonBlocking) == cudaSuccess);
    mfront_gpmp_Statistics* statistics = nullptr;
    TFEL_TESTS_ASSERT(cudaMallocManaged(&statistics,
                                        sizeof(mfront_gpmp_Statistics)) ==
                      cudaSuccess);
    mfront_gpmp_initializeStatistics(statistics);
    const auto T = thrust::host_vector<double>{300, -1, 500, 800};
    const thrust::device_vector<double> d_values(T);
    thrust::device_vector<double> d_T(T.size(), 0.);
    thrust::device_vector<double> d_E(T.size());
    // the temperatures are written on the stream after a delay
    writeAfterDelay<<<1, 4, 0, stream>>>(
        thrust::raw_pointer_cast(d_T.data()),
        thrust::raw_pointer_cast(d_values.data()), 100'000'000);
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto options =
        mfront_gpmp_ExecutionOptions{.queue = stream, .statistics = statistics};
    auto output = mfront_gpmp_OutputStatus{};
    Inconel600_YoungModulus(&output, thrust::raw_pointer_cast(d_E.data()), 1,
                            args.data(), args_strides.data(), 1, 4, policy,
                            &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    // the function returns once the computations are submitted
    TFEL_TESTS_ASSERT(cudaStreamQuery(stream) == cudaErrorNotReady);
    TFEL_TESTS_ASSERT(cudaStreamSynchronize(stream) == cudaSuccess);
    Inconel600_YoungModulus_generic_parallel_report(&output, statistics,
                                                    policy);
    TFEL_TESTS_CHECK_EQUAL(output.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output.first_failure_index, mfront_gmp_size_type{1});
    const auto E = thrust::host_vector<double>(d_E);
    for (const auto i : {0, 2, 3}) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
    }
    TFEL_TESTS_ASSERT(cudaFree(statistics) == cudaSuccess);
    TFEL_TESTS_ASSERT(cudaStreamDestroy(stream) == cudaSuccess);
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
    const auto T =
        thrust::host_vector<double>(T_values.begin(), T_values.end());
    thrust::device_vector<double> d_T(T);
    thrust::device_vector<double> d_v(T.size());
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args =
        std::array<const double *, 1u>{thrust::raw_pointer_cast(d_T.data())};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    mp(&output, thrust::raw_pointer_cast(d_v.data()), 1, args.data(),
       args_strides.data(), args.size(), T.size(), policy, nullptr);
    const auto v = thrust::host_vector<double>(d_v);
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
