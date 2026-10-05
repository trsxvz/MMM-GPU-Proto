#include <cstdio>
template <typename I>
__global__ void ne(double* o) { double s = 0.; for (I k = 0; k != 4; ++k) s += 1.; o[threadIdx.x] = s; }
template <typename I>
__global__ void ne_arr(double* o) { double a[6]; for (I k = 0; k != 6; ++k) a[k] = 1.; double s = 0.; for (I k = 0; k != 6; ++k) s += a[k] * 0.5; o[threadIdx.x] = s; }
template <typename I>
__global__ void ne_post(double* o) { double s = 0.; for (I k = 0; k != 4; k++) s += 1.; o[threadIdx.x] = s; }
template <typename I>
__global__ void ne_mul(double* o) { double s = 1.; for (I k = 0; k != 4; ++k) s *= 2.; o[threadIdx.x] = s; }
template <typename I>
void run(const char* n, double* o) {
  ne<I><<<1, 32>>>(o); cudaDeviceSynchronize(); std::printf("%-15s sum 4: %g", n, o[0]);
  ne_arr<I><<<1, 32>>>(o); cudaDeviceSynchronize(); std::printf(" | arr 3: %g", o[0]);
  ne_post<I><<<1, 32>>>(o); cudaDeviceSynchronize(); std::printf(" | post 4: %g", o[0]);
  ne_mul<I><<<1, 32>>>(o); cudaDeviceSynchronize(); std::printf(" | mul 16: %g\n", o[0]);
}
int main() {
  double* o; cudaMallocManaged(&o, 8 * 32);
  run<int>("int", o); run<unsigned>("unsigned", o); run<unsigned short>("unsigned short", o);
  run<long>("long", o); run<std::size_t>("size_t", o); run<short>("short", o);
}
