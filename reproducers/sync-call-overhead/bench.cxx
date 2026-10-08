// Time of a synchronous call of the generic-parallel interface (CUDA
// backend) on a batch of n points of an elastic behaviour, data in managed
// memory. usage: bench library n calls
#include <dlfcn.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cuda_runtime_api.h>
#include "MFront/GenericParallelBehaviour/BatchData.h"

using Function = int (*)(mfront_gpb_OutputStatus*, const mfront_gpb_BatchData*);

static double* managed(const std::size_t n, const double v) {
  void* p = nullptr;
  if (cudaMallocManaged(&p, n * sizeof(double)) != cudaSuccess) {
    std::exit(EXIT_FAILURE);
  }
  std::fill_n(static_cast<double*>(p), n, v);
  return static_cast<double*>(p);
}

int main(const int argc, const char* const* const argv) {
  if (argc != 4) {
    return EXIT_FAILURE;
  }
  const auto n = static_cast<std::size_t>(std::atol(argv[2]));
  const auto calls = std::atoi(argv[3]);
  auto* const l = dlopen(argv[1], RTLD_NOW);
  const auto f = reinterpret_cast<Function>(dlsym(l, "Elasticity_Tridimensional"));
  if (f == nullptr) {
    std::printf("%s\n", dlerror());
    return EXIT_FAILURE;
  }
  auto* const e0 = managed(6 * n, 0);
  auto* const e1 = managed(6 * n, 1e-3);
  auto* const s0 = managed(6 * n, 0);
  auto* const s1 = managed(6 * n, 0);
  auto* const K = managed(36 * n, 0);
  auto* const mp = managed(2, 150e9);
  mp[1] = 0.3;
  auto* const T = managed(1, 293.15);
  auto d = mfront_gpb_BatchData{};
  d.number_of_points = n;
  d.dt = 1;
  d.options[0] = 4;  // consistent tangent operator
  d.rdt = 1;
  d.tangent_operators = {K, 36, 1};
  d.s0.gradients = {e0, 6, 1};
  d.s0.thermodynamic_forces = {s0, 6, 1};
  d.s0.material_properties = {mp, 0, 1};
  d.s0.external_state_variables = {T, 0, 1};
  d.s1.gradients = {e1, 6, 1};
  d.s1.thermodynamic_forces = {s1, 6, 1};
  d.s1.material_properties = {mp, 0, 1};
  d.s1.external_state_variables = {T, 0, 1};
  auto o = mfront_gpb_OutputStatus{};
  for (int i = 0; i != 20; ++i) {
    if (f(&o, &d) != 1) {
      std::printf("status %d: %s\n", o.status, o.msg);
      return EXIT_FAILURE;
    }
  }
  const auto start = std::chrono::steady_clock::now();
  for (int i = 0; i != calls; ++i) {
    f(&o, &d);
  }
  const auto stop = std::chrono::steady_clock::now();
  const auto us =
      std::chrono::duration<double, std::micro>(stop - start).count() / calls;
  std::printf("n=%zu %.2f us per call (status %d, s1[0]=%g)\n", n, us,
              o.status, s1[0]);
  return EXIT_SUCCESS;
}
