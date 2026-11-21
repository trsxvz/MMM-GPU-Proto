// A kernel parameter passed by value whose member functions take the
// address, the address escaping (as the function objects of the
// generic-parallel interface of MFront), is copied by nvcc in the stack
// (local memory) of each thread, unless it is declared __grid_constant__
// (compute capability 7.0 or later). The run script counts the local stores of the kernel in both
// cases and checks that the results are the same.
#include <cstdio>
#include <cuda_runtime.h>

#ifdef USE_GRID_CONSTANT
#define GRID_CONSTANT __grid_constant__
#else
#define GRID_CONSTANT
#endif

// function object holding data (as the fields and the parameters of a
// behaviour): its call operator reads the data through `this`
struct Function {
  float values[256];
  __device__ float operator()(const unsigned i) const {
    const float* const p = this->values;
    return p[i % 256] + p[(7 * i) % 256];
  }
};

// view of a point holding the address of the function object, given to a
// function which is not inlined (as the views of the points of MFront,
// given to deeply nested functions): the address escapes
struct View {
  const Function* f;
  unsigned i;
};

__device__ __noinline__ float evaluate(const View& v) { return (*v.f)(v.i); }

__global__ void apply(const GRID_CONSTANT Function f, float* const out,
                      const unsigned n) {
  const auto i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    out[i] = evaluate(View{&f, i});
  }
}

int main() {
  auto f = Function{};
  for (unsigned i = 0; i != 256; ++i) {
    f.values[i] = static_cast<float>(i);
  }
  constexpr unsigned n = 1u << 20;
  float* out = nullptr;
  cudaMalloc(&out, n * sizeof(float));
  apply<<<n / 256, 256>>>(f, out, n);
  float h[2];
  cudaMemcpy(h, out + 1000, 2 * sizeof(float), cudaMemcpyDeviceToHost);
  std::printf("%g %g\n", h[0], h[1]);
  cudaFree(out);
  return 0;
}
