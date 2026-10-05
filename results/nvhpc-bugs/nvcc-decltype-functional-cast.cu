// nvcc 13.0: `decltype(s){}` with `const auto s = c.size();` in a function
// template used in device code is rewritten by cudafe++ for the host
// compiler as `const unsigned long{}`, which is not valid C++.
// nvcc -std=c++20 -arch=native -c decltype-functional-cast.cu
#include <array>

template <typename Container>
__host__ __device__ constexpr auto findIndex(const Container& c) {
  const auto s = c.size();
  auto i = decltype(s){};
  while (i + 1 != s) {
    ++i;
  }
  return i;
}

__global__ void kernel(unsigned long* r) {
  constexpr auto a = std::array<double, 3>{1, 2, 3};
  *r = findIndex(a);
}
