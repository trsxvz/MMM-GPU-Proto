// nvcc (CUDA 13.0): the device copy of a function scope
// `static constexpr` object, whose constructor is called, is filled with
// zeros, without any diagnostic. Expected output: 1.5 3 4.5 6.
// nvcc -std=c++20 -arch=native --expt-relaxed-constexpr \
//   nvcc-static-constexpr-object.cu && ./a.out
#include <cstdio>
struct Data {
  __host__ __device__ constexpr Data() : v{1.5, 3, 4.5, 6} {}
  double v[4];
  __host__ __device__ static const Data& get() {
    static constexpr Data d;
    return d;
  }
};
__global__ void k(double* r) { r[threadIdx.x] = Data::get().v[threadIdx.x]; }
int main() {
  double* r;
  cudaMallocManaged(&r, 4 * sizeof(double));
  k<<<1, 4>>>(r);
  cudaDeviceSynchronize();
  std::printf("%g %g %g %g (expected 1.5 3 4.5 6)\n", r[0], r[1], r[2], r[3]);
}
