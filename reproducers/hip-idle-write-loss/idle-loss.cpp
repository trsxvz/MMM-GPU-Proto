// Writes of a finished kernel to device memory are lost if the GPU stays
// idle for a short time before the next operation.
//
// usage: idle-loss [idle time in microseconds] [sync]
// With "sync", the null stream is synchronised right after the launch.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <hip/hip_runtime.h>

#define CHECK(call)                                                  \
  do {                                                               \
    const hipError_t e = (call);                                     \
    if (e != hipSuccess) {                                           \
      std::fprintf(stderr, "%s: %s\n", #call, hipGetErrorString(e)); \
      std::exit(2);                                                  \
    }                                                                \
  } while (0)

__global__ void set_kernel(int* p, int n, int v) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    p[i] = v;
  }
}

int main(int argc, char** argv) {
  const int idle_us = argc > 1 ? std::atoi(argv[1]) : 1000;
  const bool sync = argc > 2 && std::strcmp(argv[2], "sync") == 0;
  constexpr int n = 256;
  constexpr int trials = 20;
  int* d;
  CHECK(hipMalloc(&d, n * sizeof(int)));
  int lost = 0;
  for (int t = 1; t <= trials; ++t) {
    hipLaunchKernelGGL(set_kernel, dim3(1), dim3(n), 0, 0, d, n, t);
    CHECK(hipGetLastError());
    if (sync) {
      CHECK(hipStreamSynchronize(nullptr));
    }
    std::this_thread::sleep_for(std::chrono::microseconds(idle_us));
    int h[n];
    CHECK(hipMemcpy(h, d, n * sizeof(int), hipMemcpyDeviceToHost));
    if (h[0] != t) {
      ++lost;
    }
  }
  std::printf("idle %d us%s: %d / %d trials lost the writes of the kernel\n", idle_us,
              sync ? ", synchronised" : "", lost, trials);
  return lost == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
