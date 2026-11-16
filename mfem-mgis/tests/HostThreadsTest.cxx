/*!
 * \file   HostThreadsTest.cxx
 * \brief  Test of the loops of the host shared between the threads of
 * OpenMP (`batch::internals::forEachThreadRange`, internal header of the
 * library, compiled here as the library compiles it): whatever the number
 * of threads (`OMP_NUM_THREADS`), the ranges given to the function are not
 * empty, follow each other, have the same size up to one index and cover
 * the indices once; nothing is called without index; an exception thrown by
 * a call is thrown again. The expected number of ranges of a large loop,
 * which is the number of threads when the library is built with the option
 * `enable-openmp` and one otherwise, is given as argument.
 * \date   16/11/2026
 */

#include <mutex>
#include <string>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "HostThreads.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \return the ranges given by `forEachThreadRange` for a number of indices,
 * sorted, after checking that each index is treated once
 * \param[in] n: number of indices
 */
template <typename Index>
static std::vector<std::pair<Index, Index>> getRanges(const Index n) {
  using mfem_mgis::batch::internals::forEachThreadRange;
  const auto name = std::to_string(n) + " indices";
  auto ranges = std::vector<std::pair<Index, Index>>{};
  auto calls = std::vector<char>(static_cast<std::size_t>(n), 0);
  auto mutex = std::mutex{};
  forEachThreadRange(n, [&](const Index begin, const Index end) {
    for (auto i = begin; i != end; ++i) {
      ++calls[static_cast<std::size_t>(i)];
    }
    const auto lock = std::lock_guard{mutex};
    ranges.emplace_back(begin, end);
  });
  check(std::all_of(calls.begin(), calls.end(),
                    [](const char c) { return c == 1; }),
        name + ": indices not treated once");
  std::sort(ranges.begin(), ranges.end());
  auto next = Index{};
  auto smallest = n;
  auto largest = Index{};
  for (const auto& [begin, end] : ranges) {
    check(begin == next, name + ": ranges which do not follow each other");
    check(begin < end, name + ": empty range");
    smallest = std::min(smallest, static_cast<Index>(end - begin));
    largest = std::max(largest, static_cast<Index>(end - begin));
    next = end;
  }
  check(next == n, name + ": indices not covered");
  check(ranges.empty() || (largest - smallest <= 1),
        name + ": ranges of different sizes");
  return ranges;
}  // end of getRanges

int main(const int argc, const char* const* const argv) {
  using mfem_mgis::batch::internals::forEachThreadRange;
  if (argc != 2) {
    std::cerr << "usage: HostThreadsTest <expected number of ranges>\n";
    return EXIT_FAILURE;
  }
  const auto expected = static_cast<std::size_t>(std::atoi(argv[1]));
  check(getRanges(0).empty(), "function called without index");
  for (const auto n : {1, 2, 3, 7, 16, 17, 1000}) {
    const auto ranges = getRanges(n);
    check(ranges.size() == std::min(expected, static_cast<std::size_t>(n)),
          std::to_string(n) + " indices: " + std::to_string(ranges.size()) +
              " ranges, " + std::to_string(expected) + " expected");
  }
  check(getRanges(std::size_t{100000}).size() == expected,
        "unexpected number of ranges for indices of another type");
  // more indices than an `int` holds: the bounds of the ranges do not
  // overflow
  {
    const auto n = std::int64_t{1} << 40;
    auto total = std::int64_t{};
    auto mutex = std::mutex{};
    forEachThreadRange(n,
                       [&](const std::int64_t begin, const std::int64_t end) {
                         const auto lock = std::lock_guard{mutex};
                         check((begin >= 0) && (begin < end) && (end <= n),
                               "invalid range of a large loop");
                         total += end - begin;
                       });
    check(total == n, "large loop not covered");
  }
  // an exception thrown by the calls is thrown again, once all returned
  {
    auto calls = std::size_t{};
    auto mutex = std::mutex{};
    auto caught = false;
    try {
      forEachThreadRange(1000, [&](const int, const int) {
        {
          const auto lock = std::lock_guard{mutex};
          ++calls;
        }
        throw std::runtime_error("failure of a range");
      });
    } catch (std::runtime_error& e) {
      caught = std::string{e.what()} == "failure of a range";
    }
    check(caught, "exception not thrown again");
    check(calls == expected, "calls missing after an exception");
  }
  std::cout << expected
            << " ranges: " << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
