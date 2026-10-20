// Minimum and maximum reductions of OpenACC on 64 bits integers. One point
// out of 1000 (index k modulo 1000) gives its index, the others the
// identity of the reduction: the expected minimum is k, the expected
// maximum is the last such index.
#include <cstdio>
#include <initializer_list>
#include <limits>

template <typename T>
T minimum(const long n, const long k) {
  auto m = std::numeric_limits<T>::max();
#pragma acc parallel loop reduction(min : m)
  for (long i = 0; i < n; ++i) {
    const auto v =
        (i % 1000 == k) ? static_cast<T>(i) : std::numeric_limits<T>::max();
    m = v < m ? v : m;
  }
  return m;
}

template <typename T>
T maximum(const long n, const long k) {
  auto m = T{0};
#pragma acc parallel loop reduction(max : m)
  for (long i = 0; i < n; ++i) {
    const auto v = (i % 1000 == k) ? static_cast<T>(i) : T{0};
    m = v > m ? v : m;
  }
  return m;
}

int main() {
  constexpr auto n = 100000l;
  constexpr auto k = 7l;
  const auto results = {
      static_cast<double>(minimum<unsigned long long>(n, k)),
      static_cast<double>(minimum<long long>(n, k)),
      static_cast<double>(minimum<unsigned long>(n, k)),
      static_cast<double>(minimum<long>(n, k)),
      static_cast<double>(minimum<unsigned int>(n, k)),
      static_cast<double>(minimum<double>(n, k)),
      static_cast<double>(minimum<float>(n, k)),
      static_cast<double>(maximum<unsigned long long>(n, k))};
  const char* const names[] = {"min unsigned long long", "min long long",
                               "min unsigned long",      "min long",
                               "min unsigned int",       "min double",
                               "min float",              "max unsigned long long"};
  const double expected[] = {7, 7, 7, 7, 7, 7, 7, 99007};
  auto status = 0;
  auto i = 0;
  for (const auto r : results) {
    const auto ok = r == expected[i];
    std::printf("%-24s %8.0f %s\n", names[i], r, ok ? "OK" : "WRONG");
    status = ok ? status : 1;
    ++i;
  }
  return status;
}
