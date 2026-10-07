// nvc++ -mp=gpu: a lambda capturing an array of pointers by value, mapped
// by `map(to : kernel)` in a target region, reads wrong pointers on the
// GPU. A lambda capturing the same pointer directly, or a structure holding
// it, works. An array of integers works.
// nvc++ -O2 -mp=gpu -gpu=mem:managed -std=c++20 omp-lambda-pointer-array.cxx
#include <array>
#include <cstdio>
#include <vector>

int main() {
  auto x = std::vector<double>{1, 2, 3, 4};
  auto y = std::vector<double>(4, -1);
  auto* const out = y.data();
  const auto args = std::array<const double*, 1>{x.data()};
  auto kernel = [out, args](const int i) { out[i] = args[0][i]; };
#pragma omp target teams distribute parallel for map(to : kernel)
  for (int i = 0; i < 4; ++i) {
    kernel(i);
  }
  const auto ok = (y[0] == 1) && (y[1] == 2) && (y[2] == 3) && (y[3] == 4);
  std::printf("%g %g %g %g: %s\n", y[0], y[1], y[2], y[3],
              ok ? "OK" : "FAILED");
  return ok ? 0 : 1;
}
