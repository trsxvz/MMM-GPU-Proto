// nvcc (CUDA 13.0): a constexpr function using the `__builtin_constant_p`
// built-in (a host function for nvcc) is silently compiled as a function
// doing nothing in device code (--expt-relaxed-constexpr, -O0 to -O3).
// libstdc++ (GCC 14) uses this built-in in `std::__advance`: `std::advance`,
// `std::next` and `std::prev` do nothing in device code. Expected output:
// 2 2 1 1.
// nvcc -std=c++20 -arch=native --expt-relaxed-constexpr \
//   nvcc-builtin-constant-p.cu && ./a.out
#include <cstdio>
#include <iterator>
// random access version of `std::__advance` of libstdc++ 14
template <typename I, typename D>
constexpr void advance(I& i, const D n) {
  if (__builtin_constant_p(n) && n == 1) {
    ++i;
  } else if (__builtin_constant_p(n) && n == -1) {
    --i;
  } else {
    i += n;
  }
}
__global__ void k(long* r, const int* a, const long n) {
  const int* p = a;
  advance(p, n);
  r[0] = p - a;
  r[1] = std::next(a, n) - a;
  r[2] = std::prev(a + 2) - a;
  r[3] = std::prev(a + n) - a;
}
int main() {
  long* r;
  int* a;
  cudaMallocManaged(&r, 4 * sizeof(long));
  cudaMallocManaged(&a, 4 * sizeof(int));
  k<<<1, 1>>>(r, a, 2);
  cudaDeviceSynchronize();
  std::printf("%ld %ld %ld %ld (expected 2 2 1 1)\n", r[0], r[1], r[2], r[3]);
}
