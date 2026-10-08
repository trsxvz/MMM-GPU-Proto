// nvc++ -mp=gpu: in a lambda mapped by `map(to : kernel)` in a target
// region, a nested lambda capturing by reference a pointer captured by the
// kernel does not write through this pointer, without diagnostic (-O0 to
// -O3). The same nested lambda capturing by copy works.
// nvc++ -O2 -mp=gpu -gpu=mem:managed -std=c++20 omp-lambda-reference-capture.cxx
#include <cstdio>
#include <vector>

int main() {
  auto x = std::vector<double>{1, 2, 3, 4};
  auto y = std::vector<double>(4, -1);
  auto* const out = y.data();
  const auto* const in = x.data();
  auto kernel = [out, in](const int i) {
    // replacing [&] by [=] gives the expected results
    [&]() { out[i] = in[i]; }();
  };
#pragma omp target teams distribute parallel for map(to : kernel)
  for (int i = 0; i < 4; ++i) {
    kernel(i);
  }
  const auto ok = (y[0] == 1) && (y[1] == 2) && (y[2] == 3) && (y[3] == 4);
  std::printf("%g %g %g %g: %s\n", y[0], y[1], y[2], y[3],
              ok ? "OK" : "FAILED");
  return ok ? 0 : 1;
}
