#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cstring>
#include "kernels.hxx"

#define CK(x)                                                         \
  do {                                                                \
    cudaError_t ck_e_ = (x);                                              \
    if (ck_e_ != cudaSuccess) {                                           \
      std::printf("CUDA error %s at %d\n", cudaGetErrorString(ck_e_), __LINE__); \
      std::exit(1);                                                   \
    }                                                                 \
  } while (0)

template <int G, typename A>
__global__ void per_point(const double* K, const double* g, double* y, A a,
                          idx_t n) {
  const idx_t o = blockIdx.x * (idx_t)blockDim.x + threadIdx.x;
  if (o < n) point<G>(K, g, y, a, o);
}

template <int G, typename A>
__global__ void per_element(const double* K, const double* g, double* y, A a,
                            idx_t ne, int nq) {
  const idx_t e = blockIdx.x * (idx_t)blockDim.x + threadIdx.x;
  if (e >= ne) return;
  for (int q = 0; q != nq; ++q) point<G>(K, g, y, a, e * nq + q);
}

template <int G, typename A>
float run(bool by_point, const double* K, const double* g, double* y, A a,
          idx_t n, int nq) {
  cudaEvent_t b, e;
  CK(cudaEventCreate(&b));
  CK(cudaEventCreate(&e));
  std::vector<float> t;
  const int bs = 256;
  for (int i = 0; i != 23; ++i) {
    CK(cudaEventRecord(b));
    if (by_point) {
      per_point<G><<<(n + bs - 1) / bs, bs>>>(K, g, y, a, n);
    } else {
      const idx_t ne = n / nq;
      per_element<G><<<(ne + bs - 1) / bs, bs>>>(K, g, y, a, ne, nq);
    }
    CK(cudaEventRecord(e));
    CK(cudaEventSynchronize(e));
    float ms;
    CK(cudaEventElapsedTime(&ms, b, e));
    if (i >= 3) t.push_back(ms);
  }
  std::sort(t.begin(), t.end());
  return t[t.size() / 2];
}

template <int G>
void bench(idx_t n, int nq, idx_t ld) {
  constexpr int NK = G * G;
  std::vector<double> Ka(n * NK), ga(n * G), Ks(ld * NK), gs(ld * G);
  for (idx_t o = 0; o != n; ++o) {
    for (int c = 0; c != NK; ++c) {
      const double v = 1.0 + 1e-3 * ((o * 7 + c * 13) % 101);
      Ka[o * NK + c] = v;
      Ks[c * ld + o] = v;
    }
    for (int c = 0; c != G; ++c) {
      const double v = 1.0 + 1e-2 * ((o * 3 + c * 5) % 37);
      ga[o * G + c] = v;
      gs[c * ld + o] = v;
    }
  }
  double *dKa, *dga, *dKs, *dgs, *dy;
  CK(cudaMalloc(&dKa, n * NK * 8));
  CK(cudaMalloc(&dKs, ld * NK * 8));
  CK(cudaMalloc(&dga, n * G * 8));
  CK(cudaMalloc(&dgs, ld * G * 8));
  CK(cudaMalloc(&dy, ld * G * 8));
  CK(cudaMemcpy(dKa, Ka.data(), n * NK * 8, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dKs, Ks.data(), ld * NK * 8, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dga, ga.data(), n * G * 8, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dgs, gs.data(), ld * G * 8, cudaMemcpyHostToDevice));
  const double bytes = double(n) * (NK + 2 * G) * 8;
  const RuntimeAccess ra{NK, 1, G, 1};
  const RuntimeAccess rs{1, ld, 1, ld};
  const StaticAccess<false, G> sa{ld};
  const StaticAccess<true, G> ss{ld};
  // check: same values in both layouts
  std::vector<double> ya(n * G), ys(ld * G);
  run<G>(true, dKa, dga, dy, sa, n, nq);
  CK(cudaMemcpy(ya.data(), dy, n * G * 8, cudaMemcpyDeviceToHost));
  run<G>(true, dKs, dgs, dy, ss, n, nq);
  CK(cudaMemcpy(ys.data(), dy, ld * G * 8, cudaMemcpyDeviceToHost));
  bool same = true;
  for (idx_t o = 0; o != n && same; ++o)
    for (int c = 0; c != G; ++c)
      if (std::memcmp(&ya[o * G + c], &ys[c * ld + o], 8) != 0) same = false;
  std::printf("G=%d n=%lld ld=%lld nq=%d bitwise identical AoS/SoA: %s\n", G, n, ld, nq,
              same ? "yes" : "NO");
  for (int bp = 1; bp >= 0; --bp) {
    const char* m = bp ? "thread/point  " : "thread/element";
    const float t1 = run<G>(bp, dKa, dga, dy, sa, n, nq);
    const float t2 = run<G>(bp, dKa, dga, dy, ra, n, nq);
    const float t3 = run<G>(bp, dKs, dgs, dy, ss, n, nq);
    const float t4 = run<G>(bp, dKs, dgs, dy, rs, n, nq);
    std::printf(
        "  %s AoS static %7.3f ms (%5.0f GB/s) runtime %7.3f ms (%5.0f GB/s)"
        " | SoA static %7.3f ms (%5.0f GB/s) runtime %7.3f ms (%5.0f GB/s)\n",
        m, t1, bytes / t1 * 1e-6, t2, bytes / t2 * 1e-6, t3,
        bytes / t3 * 1e-6, t4, bytes / t4 * 1e-6);
  }
  cudaFree(dKa);
  cudaFree(dKs);
  cudaFree(dga);
  cudaFree(dgs);
  cudaFree(dy);
}

int main() {
  const idx_t n2 = idx_t{1} << 21;
  const idx_t n = 2000000;  // multiple of 8 and of 27? no: 2000016 below
  bench<6>(2000016, 8, 2000016);
  bench<6>(2000016, 8, 2000016 + 16);
  bench<9>(2000016, 8, 2000016);
  bench<9>(2000016, 8, 2000016 + 16);
  bench<9>(2000016, 27, 2000016);
  // copy bandwidth reference
  double *a, *b;
  const idx_t m = idx_t{1} << 27;
  CK(cudaMalloc(&a, m * 8));
  CK(cudaMalloc(&b, m * 8));
  cudaEvent_t e0, e1;
  cudaEventCreate(&e0);
  cudaEventCreate(&e1);
  std::vector<float> t;
  for (int i = 0; i != 13; ++i) {
    cudaEventRecord(e0);
    cudaMemcpy(b, a, m * 8, cudaMemcpyDeviceToDevice);
    cudaEventRecord(e1);
    cudaEventSynchronize(e1);
    float ms;
    cudaEventElapsedTime(&ms, e0, e1);
    if (i >= 3) t.push_back(ms);
  }
  std::sort(t.begin(), t.end());
  std::printf("device memcpy: %.0f GB/s\n", 2.0 * m * 8 / t[t.size() / 2] * 1e-6);
}
