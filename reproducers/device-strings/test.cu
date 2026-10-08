#include <cstdio>
#include "model.hxx"
__global__ void kernel(double* r) { r[threadIdx.x] = Model<double>{}.execute(threadIdx.x); }
int main() {
  double* r;
  cudaMallocManaged(&r, 4 * sizeof(double));
  kernel<<<1, 4>>>(r);
  if (cudaDeviceSynchronize() != cudaSuccess) { std::printf("kernel failed\n"); return 2; }
  std::printf("%g %g %g %g\n", r[0], r[1], r[2], r[3]);
  return (r[3] == 6 * 3 + 1) ? 0 : 1;
}
