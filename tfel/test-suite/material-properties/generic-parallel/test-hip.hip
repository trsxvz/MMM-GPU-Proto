/*!
 * \file   test-hip.hip
 * \brief
 * \author Thomas Helfer
 * \date   22/04/2026
 */

#include <array>
#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>
#include <initializer_list>

#include "hip/hip_runtime.h"
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"

#include "MFront/GenericMaterialProperty/MaterialProperty.h"
#include "PoissonRatioTest-generic-parallel-hip.hxx"
#include "Inconel600_YoungModulus-generic-parallel-hip.hxx"
#include "Inconel600_YoungModulus_qt-generic-parallel-hip.hxx"
#include "UO2_YoungModulus_Martin1989-generic-parallel-hip.hxx"
#include "LinearDataInterpolationTest-generic-parallel-hip.hxx"
#include "MaterialLawCallTest-generic-parallel-hip.hxx"
#include "symbols.hxx"

[[nodiscard]] std::unique_ptr<double, void (*)(double *const) noexcept>
createManagedData(const std::size_t s) {
  double *ptr;
  const auto e = hipMallocManaged(&ptr, s * sizeof(double));
  if (e != hipSuccess) {
    std::cerr << hipGetErrorString(e) << '\n';
    std::exit(EXIT_FAILURE);
  }
  // the memory returned by hipMallocManaged is not initialized
  std::fill_n(ptr, s, 0.);
  return {ptr, +[](double *p) noexcept { std::ignore = hipFree(p); }};
}  // end of createManagedData

[[nodiscard]] std::unique_ptr<double, void (*)(double *const) noexcept>
createManagedData(const std::initializer_list<double> l) {
  auto ptr = createManagedData(l.size());
  std::copy(l.begin(), l.end(), ptr.get());
  return {std::move(ptr)};
}  // end of createManagedData

[[nodiscard]] double getValue(
    const std::unique_ptr<double, void (*)(double *const) noexcept> &ptr,
    const std::size_t i) {
  return (ptr.get())[i];
}  // end of getValue

/*!
 * \brief kernel writing the given values after a delay, so that the
 * computations submitted after it on the same stream must wait for it
 * \param[out] T: written values
 * \param[in] values: values
 * \param[in] cycles: number of clock cycles of the delay
 */
__global__ void writeAfterDelay(double* const T,
                                const double* const values,
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
    TFEL_TESTS_ASSERT(checkGenericParallelSymbols("hip"));
    return this->result;
  }  // end of execute

 private:
  [[nodiscard]] static double fct(const double TK) noexcept {
    const double TC = TK - 273.15;
    return (-3.1636e-3 * TC * TC - 3.8654 * TC + 2.1421e+4) * 1e7;
  }  // end of fct

  void test14() {
    // computations on a stream given by the execution options, which is not
    // the default stream
    constexpr auto eps = double{1e-2};
    auto stream = hipStream_t{};
    TFEL_TESTS_ASSERT(hipStreamCreateWithFlags(
                          &stream, hipStreamNonBlocking) == hipSuccess);
    auto E = createManagedData(4);
    const auto values = createManagedData({300, 500, 300, 800});
    const auto T = createManagedData(4);
    // the temperatures are written on the stream, after a delay: the
    // material property must be evaluated after this kernel, without
    // synchronization by the host
    writeAfterDelay<<<1, 4, 0, stream>>>(T.get(), values.get(), 100'000'000);
    const auto args = std::array<const double *, 1u>{T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    const auto options = mfront_gpmp_ExecutionOptions{.queue = stream};
    auto output = mfront_gpmp_OutputStatus{};
    Inconel600_YoungModulus(&output, E.get(), 1, args.data(),
                            args_strides.data(), 1, 4,
                            mfront_gmp_OutOfBoundsPolicy{}, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    for (std::size_t i = 0; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(E, i) - fct(getValue(values, i))) <
                        eps);
    }
    TFEL_TESTS_ASSERT(hipStreamDestroy(stream) == hipSuccess);
  }

  void test16() {
    // stack size of the threads preserved by the execution options: an
    // error is reported if it is too small for the kernel. Otherwise, the
    // stack size is increased if required.
    constexpr auto eps = double{1e-2};
    auto initial_stack_size = std::size_t{};
    TFEL_TESTS_ASSERT(hipDeviceGetLimit(&initial_stack_size,
                                        hipLimitStackSize) == hipSuccess);
    // the kernels require at least 1024 bytes (see
    // mfront::generic_parallel::gpu::ensureStackSize)
    TFEL_TESTS_ASSERT(hipDeviceSetLimit(hipLimitStackSize, 256) ==
                      hipSuccess);
    auto stack_size = std::size_t{};
    TFEL_TESTS_ASSERT(hipDeviceGetLimit(&stack_size, hipLimitStackSize) ==
                      hipSuccess);
    TFEL_TESTS_ASSERT(stack_size < 1024);
    auto E = createManagedData(2);
    const auto T = createManagedData({300, 500});
    const auto args = std::array<const double *, 1u>{T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    auto output = mfront_gpmp_OutputStatus{};
    const auto preserve = mfront_gpmp_ExecutionOptions{
        .queue = nullptr, .statistics = nullptr, .preserve_stack_size = 1};
    Inconel600_YoungModulus(&output, E.get(), 1, args.data(),
                            args_strides.data(), 1, 2, policy, &preserve);
    TFEL_TESTS_CHECK_EQUAL(output.status, -6);
    TFEL_TESTS_ASSERT(std::string(output.msg).starts_with(
        "the stack size of the threads of the device must be at least"));
    TFEL_TESTS_ASSERT(hipDeviceGetLimit(&stack_size, hipLimitStackSize) ==
                      hipSuccess);
    TFEL_TESTS_ASSERT(stack_size < 1024);
    // default options: the stack size is increased
    Inconel600_YoungModulus(&output, E.get(), 1, args.data(),
                            args_strides.data(), 1, 2, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_ASSERT(hipDeviceGetLimit(&stack_size, hipLimitStackSize) ==
                      hipSuccess);
    TFEL_TESTS_ASSERT(stack_size >= 1024);
    for (std::size_t i = 0; i != 2; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(E, i) - fct(getValue(T, i))) < eps);
    }
    TFEL_TESTS_ASSERT(hipDeviceSetLimit(hipLimitStackSize,
                                        initial_stack_size) == hipSuccess);
  }

  void test15() {
    // asynchronous evaluation on a stream: the statistics, in managed memory,
    // are given by the execution options and reported once the stream is
    // synchronized. The second temperature is out of its physical bounds.
    constexpr auto eps = double{1e-2};
    auto stream = hipStream_t{};
    TFEL_TESTS_ASSERT(hipStreamCreateWithFlags(
                          &stream, hipStreamNonBlocking) == hipSuccess);
    mfront_gpmp_Statistics* statistics = nullptr;
    TFEL_TESTS_ASSERT(hipMallocManaged(&statistics,
                                        sizeof(mfront_gpmp_Statistics)) ==
                      hipSuccess);
    mfront_gpmp_initializeStatistics(statistics);
    auto E = createManagedData(4);
    const auto values = createManagedData({300, -1, 500, 800});
    const auto T = createManagedData(4);
    // the temperatures are written on the stream after a delay
    writeAfterDelay<<<1, 4, 0, stream>>>(T.get(), values.get(), 100'000'000);
    const auto args = std::array<const double *, 1u>{T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto options =
        mfront_gpmp_ExecutionOptions{.queue = stream, .statistics = statistics};
    auto output = mfront_gpmp_OutputStatus{};
    Inconel600_YoungModulus(&output, E.get(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, &options);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    // the function returns once the computations are submitted
    TFEL_TESTS_ASSERT(hipStreamQuery(stream) == hipErrorNotReady);
    TFEL_TESTS_ASSERT(hipStreamSynchronize(stream) == hipSuccess);
    Inconel600_YoungModulus_generic_parallel_report(&output, statistics,
                                                    policy);
    TFEL_TESTS_CHECK_EQUAL(output.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output.first_failure_index, mfront_gmp_size_type{1});
    for (const auto i : {0, 2, 3}) {
      TFEL_TESTS_ASSERT(std::abs(getValue(E, i) - fct(getValue(values, i))) <
                        eps);
    }
    TFEL_TESTS_ASSERT(hipFree(statistics) == hipSuccess);
    TFEL_TESTS_ASSERT(hipStreamDestroy(stream) == hipSuccess);
  }

  void test1() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = createManagedData(4);
    const auto T = createManagedData({300, 500, 300, 800});
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, E.get(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, nullptr);
    //
    auto E_qt = createManagedData(4);
    Inconel600_YoungModulus_qt(&output_qt, E_qt.get(), 1, args.data(),
                               args_strides.data(), 1, 4, policy, nullptr);
    //
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    for (std::size_t i = 0; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(E, i) - fct(getValue(T, i))) < eps);
      TFEL_TESTS_ASSERT(std::abs(getValue(E_qt, i) - fct(getValue(T, i))) <
                        eps);
    }
  }

  void test2() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = createManagedData(4);
    auto E_qt = createManagedData(4);
    const auto T = createManagedData({300});
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{0};
    //
    Inconel600_YoungModulus(&output, E.get(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt(&output_qt, E_qt.get(), 1, args.data(),
                               args_strides.data(), 1, 4, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    for (std::size_t i = 0; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(E, i) - fct(*T)) < eps);
      TFEL_TESTS_ASSERT(std::abs(getValue(E_qt, i) - fct(*T)) < eps);
    }
  }

  void test3() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = createManagedData(1);
    auto E_qt = createManagedData(1);
    const auto T = createManagedData({300});
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{0};
    Inconel600_YoungModulus(&output, E.get(), 0, args.data(),
                            args_strides.data(), 1, 1, policy, nullptr);
    Inconel600_YoungModulus_qt(&output_qt, E_qt.get(), 0, args.data(),
                               args_strides.data(), 1, 1, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    TFEL_TESTS_ASSERT(std::abs(*E - fct(*T)) < eps);
    TFEL_TESTS_ASSERT(std::abs(*E_qt - fct(*T)) < eps);
  }

  void test4() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = createManagedData(4);
    auto E_qt = createManagedData(4);
    const auto T = createManagedData({300, 500, 300, 800});
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.get()};
    Inconel600_YoungModulus2(&output, E.get(), args.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt2(&output_qt, E_qt.get(), args.data(), 1, 4,
                                policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    for (std::size_t i = 0; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(E, i) - fct(getValue(T, i))) < eps);
      TFEL_TESTS_ASSERT(std::abs(getValue(E_qt, i) - fct(getValue(T, i))) <
                        eps);
    }
  }

  // check detection of physical bounds violation
  void test5() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = createManagedData(4);
    auto E_qt = createManagedData(4);
    const auto T = createManagedData({-300, 500, 300, 800});
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, E.get(), 1, args.data(),
                            args_strides.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt(&output_qt, E_qt.get(), 1, args.data(),
                               args_strides.data(), 1, 4, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output.first_failure_index, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, -1);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output_qt.first_failure_index, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, -1);
    TFEL_TESTS_ASSERT(std::abs(*E) < eps);
    TFEL_TESTS_ASSERT(std::abs(*E_qt) < eps);
    for (std::size_t i = 1; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(E, i) - fct(getValue(T, i))) < eps);
      TFEL_TESTS_ASSERT(std::abs(getValue(E_qt, i) - fct(getValue(T, i))) <
                        eps);
    }
  }

  void test6() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = createManagedData(4);
    auto E_qt = createManagedData(4);
    const auto T = createManagedData({-300, 500, 300, 800});
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.get()};
    Inconel600_YoungModulus2(&output, E.get(), args.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt2(&output_qt, E_qt.get(), args.data(), 1, 4,
                                policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output.first_failure_index, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, -1);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output_qt.first_failure_index, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, -1);
    TFEL_TESTS_ASSERT(std::abs(*E) < eps);
    TFEL_TESTS_ASSERT(std::abs(*E_qt) < eps);
    for (std::size_t i = 1; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(E, i) - fct(getValue(T, i))) < eps);
      TFEL_TESTS_ASSERT(std::abs(getValue(E_qt, i) - fct(getValue(T, i))) <
                        eps);
    }
  }

  void test7() {
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = createManagedData(4 * 3);
    auto E_qt = createManagedData(4 * 3);
    const auto T =
        createManagedData({300, -300, 500, -500, 300, -300, 800, -800});
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{2};
    Inconel600_YoungModulus(&output, E.get(), 3, args.data(),
                            args_strides.data(), 1, 4, policy, nullptr);
    Inconel600_YoungModulus_qt(&output_qt, E_qt.get(), 3, args.data(),
                               args_strides.data(), 1, 4, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output_qt.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output_qt.bounds_status, 0);
    for (std::size_t i = 0; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(E, i * 3) - fct(getValue(T, i * 2))) <
                        eps);
      TFEL_TESTS_ASSERT(
          std::abs(getValue(E_qt, i * 3) - fct(getValue(T, i * 2))) < eps);
    }
  }

  void test8() {
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto E = createManagedData(1);
    const auto T = createManagedData({300});
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    Inconel600_YoungModulus(&output, E.get(), 1, args.data(),
                            args_strides.data(), 3, 1, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -5);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "invalid number of arguments (3 given, 1 expected)");
    //
    Inconel600_YoungModulus(&output, E.get(), 0, args.data(),
                            args_strides.data(), 1, 1, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, -7);
    TFEL_TESTS_CHECK_EQUAL(std::string(output.msg),
                           "if the output is uniform, all the arguments must "
                           "be uniform (i.e. their strides must be null)");
    //
    const auto args_strides2 = std::array<mfront_gmp_size_type, 1u>{0};
    Inconel600_YoungModulus(&output, E.get(), 0, args.data(),
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
    auto nu = createManagedData(4);
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    PoissonRatioTest(&output, nu.get(), 1, nullptr, nullptr, 0, 4, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(nu, i) - 0.39991) < eps);
    }
    nu = createManagedData(1);
    PoissonRatioTest(&output, nu.get(), 0, nullptr, nullptr, 0, 1, policy, nullptr);
    TFEL_TESTS_ASSERT(std::abs(*nu - 0.39991) < eps);
  }

  void test10() {
    constexpr auto eps = double{1e-14};
    auto output = mfront_gpmp_OutputStatus{};
    auto output_qt = mfront_gpmp_OutputStatus{};
    auto nu = createManagedData(4);
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    //
    PoissonRatioTest2(&output, nu.get(), nullptr, 0, 4, policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(nu, i) - 0.39991) < eps);
    }
    nu = createManagedData(1);
    PoissonRatioTest2(&output, nu.get(), nullptr, 0, 1, policy, nullptr);
    TFEL_TESTS_ASSERT(std::abs(*nu - 0.39991) < eps);
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
    auto E = createManagedData(4);
    const auto T = createManagedData({300, 500, 300, 800});
    const auto f = createManagedData({0, 0.1, 0.15, 0.12});
    //
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 2u>{T.get(), f.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 2u>{1, 1};
    UO2_YoungModulus_Martin1989(&output, E.get(), 1, args.data(),
                                args_strides.data(), args.size(), 4, policy, nullptr);
    //
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != 4; ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(E, i) -
                                 young(getValue(T, i), getValue(f, i))) < eps);
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
    auto v = createManagedData(T_values.size());
    const auto T = createManagedData(T_values.size());
    std::copy(T_values.begin(), T_values.end(), T.get());
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto args = std::array<const double *, 1u>{T.get()};
    const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
    mp(&output, v.get(), 1, args.data(), args_strides.data(), args.size(),
       T_values.size(), policy, nullptr);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{0});
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, 0);
    for (std::size_t i = 0; i != T_values.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(getValue(v, i) - references[i]) < eps);
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
