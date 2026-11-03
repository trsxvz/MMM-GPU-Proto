// Matrix free action of a small strain tangent (3D, hexahedra Q1: nd = 8,
// nq = 8), as the fused kernel of mfem-mgis: per element, per point,
// g = sym(grad u_e), y = K g, f_e += w B^T y.
// Mapping E: one thread per element (current kernel).
// Mapping Q: one thread per point, reduction over the points of the element
// with warp shuffles (nq = 8 divides the warp size).
// Layouts of K: AoS (o * 36 + c) or SoA (c * ld + o), compile time.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>

#define CK(x)                                                              \
  do {                                                                     \
    cudaError_t ck_e_ = (x);                                               \
    if (ck_e_ != cudaSuccess) {                                            \
      std::printf("CUDA error %s at %d\n", cudaGetErrorString(ck_e_),      \
                  __LINE__);                                               \
      std::exit(1);                                                        \
    }                                                                      \
  } while (0)

#ifndef NQV
#define NQV 8
#endif
#ifndef GV
#define GV 6
#endif
constexpr int ND = 8, NQ = NQV, V = 3, G = GV, NK = G * G;
using idx_t = long long;

template <bool soa>
__device__ inline idx_t kidx(idx_t o, int c, idx_t ld) {
  if constexpr (soa) {
    return c * ld + o;
  } else {
    return o * NK + c;
  }
}

// contribution of point q of element e: fl[node * V + i]
template <bool soa>
__device__ inline void point(double* fl, const double* __restrict__ K,
                             const double* __restrict__ dN,
                             const double* __restrict__ ue, double w, idx_t o,
                             int q, idx_t ld) {
  // gradient of u at q
  double gu[3][3] = {};
  for (int a = 0; a != ND; ++a) {
    for (int j = 0; j != 3; ++j) {
      const double d = dN[(q * ND + a) * 3 + j];
      for (int i = 0; i != 3; ++i) gu[i][j] += ue[a + ND * i] * d;
    }
  }
#if GV == 6
  const double g[G] = {gu[0][0], gu[1][1], gu[2][2], gu[0][1] + gu[1][0],
                       gu[0][2] + gu[2][0], gu[1][2] + gu[2][1]};
#else
  const double g[G] = {gu[0][0], gu[1][1], gu[2][2], gu[0][1], gu[1][0],
                       gu[0][2], gu[2][0], gu[1][2], gu[2][1]};
#endif
  double y[G];
  for (int r = 0; r != G; ++r) {
    double s = 0;
    for (int c = 0; c != G; ++c) s += K[kidx<soa>(o, r * G + c, ld)] * g[c];
    y[r] = w * s;
  }
  for (int a = 0; a != ND; ++a) {
    const double* d = dN + (q * ND + a) * 3;
#if GV == 6
    fl[a * V + 0] += y[0] * d[0] + y[3] * d[1] + y[4] * d[2];
    fl[a * V + 1] += y[3] * d[0] + y[1] * d[1] + y[5] * d[2];
    fl[a * V + 2] += y[4] * d[0] + y[5] * d[1] + y[2] * d[2];
#else
    fl[a * V + 0] += y[0] * d[0] + y[3] * d[1] + y[5] * d[2];
    fl[a * V + 1] += y[4] * d[0] + y[1] * d[1] + y[7] * d[2];
    fl[a * V + 2] += y[6] * d[0] + y[8] * d[1] + y[2] * d[2];
#endif
  }
}

template <bool soa>
__global__ void per_element(const double* __restrict__ K,
                            const double* __restrict__ dN,
                            const double* __restrict__ u, double* f, idx_t ne,
                            idx_t ld) {
  const idx_t e = blockIdx.x * (idx_t)blockDim.x + threadIdx.x;
  if (e >= ne) return;
  const double* ue = u + e * ND * V;
  double fl[ND * V] = {};
  for (int q = 0; q != NQ; ++q) point<soa>(fl, K, dN, ue, 0.125, e * NQ + q, q, ld);
  for (int k = 0; k != ND * V; ++k) f[e * ND * V + k] = fl[k];
}

template <bool soa>
__global__ void per_point(const double* __restrict__ K,
                          const double* __restrict__ dN,
                          const double* __restrict__ u, double* f, idx_t ne,
                          idx_t ld) {
  const idx_t t = blockIdx.x * (idx_t)blockDim.x + threadIdx.x;
  const idx_t e = t / NQ;
  const int q = static_cast<int>(t % NQ);
  double fl[ND * V] = {};
  if (e < ne) point<soa>(fl, K, dN, u + e * ND * V, 0.125, t, q, ld);
  // reduction over the NQ lanes of the element (fixed order: bitwise
  // reproducible, but different from the sequential sum of mapping E)
  for (int k = 0; k != ND * V; ++k) {
    double v = fl[k];
    for (int s = NQ / 2; s > 0; s /= 2) v += __shfl_down_sync(0xffffffff, v, s, NQ);
    fl[k] = v;
  }
  if ((e < ne) && (q == 0))
    for (int k = 0; k != ND * V; ++k) f[e * ND * V + k] = fl[k];
}


// one thread per point, deterministic reduction in shared memory: the
// contributions of the points are stored, then summed in the order of the
// points by one thread per (element, value)
constexpr int EPB = 256 / NQ;  // elements per block
template <bool soa>
__global__ void per_point_shared(const double* __restrict__ K,
                                 const double* __restrict__ dN,
                                 const double* __restrict__ u, double* f,
                                 idx_t ne, idx_t ld) {
  __shared__ double sh[EPB * NQ * ND * V];
  const int el = threadIdx.x / NQ, q = threadIdx.x % NQ;
  const idx_t e = blockIdx.x * (idx_t)EPB + el;
  double fl[ND * V] = {};
  if ((el < EPB) && (e < ne))
    point<soa>(fl, K, dN, u + e * ND * V, 0.125, e * NQ + q, q, ld);
  if (el < EPB)
    for (int k = 0; k != ND * V; ++k) sh[(el * NQ + q) * ND * V + k] = fl[k];
  __syncthreads();
  for (int t = threadIdx.x; t < EPB * ND * V; t += blockDim.x) {
    const int el2 = t / (ND * V), k = t % (ND * V);
    const idx_t e2 = blockIdx.x * (idx_t)EPB + el2;
    if (e2 >= ne) continue;
    double v = 0;
    for (int q2 = 0; q2 != NQ; ++q2) v += sh[(el2 * NQ + q2) * ND * V + k];
    f[e2 * ND * V + k] = v;
  }
}

template <typename F>
float timeit(F&& f) {
  cudaEvent_t b, e;
  cudaEventCreate(&b);
  cudaEventCreate(&e);
  std::vector<float> t;
  for (int i = 0; i != 23; ++i) {
    cudaEventRecord(b);
    f();
    cudaEventRecord(e);
    cudaEventSynchronize(e);
    float ms;
    cudaEventElapsedTime(&ms, b, e);
    if (i >= 3) t.push_back(ms);
  }
  CK(cudaGetLastError());
  std::sort(t.begin(), t.end());
  return t[t.size() / 2];
}

int main() {
  const idx_t ne = 2000016 / NQ;
  const idx_t n = ne * NQ, ld = n;
  std::vector<double> Ka(n * NK), Ks(ld * NK), u(ne * ND * V), dN(NQ * ND * 3);
  for (idx_t o = 0; o != n; ++o)
    for (int c = 0; c != NK; ++c) {
      const double v = 1.0 + 1e-3 * ((o * 7 + c * 13) % 101);
      Ka[o * NK + c] = v;
      Ks[c * ld + o] = v;
    }
  for (idx_t i = 0; i != ne * ND * V; ++i) u[i] = 1e-3 * (i % 97);
  for (int i = 0; i != NQ * ND * 3; ++i) dN[i] = 0.1 * ((i % 7) - 3);
  double *dKa, *dKs, *du, *ddN, *df;
  CK(cudaMalloc(&dKa, n * NK * 8));
  CK(cudaMalloc(&dKs, ld * NK * 8));
  CK(cudaMalloc(&du, ne * ND * V * 8));
  CK(cudaMalloc(&ddN, NQ * ND * 3 * 8));
  CK(cudaMalloc(&df, ne * ND * V * 8));
  CK(cudaMemcpy(dKa, Ka.data(), n * NK * 8, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dKs, Ks.data(), ld * NK * 8, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(du, u.data(), ne * ND * V * 8, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(ddN, dN.data(), NQ * ND * 3 * 8, cudaMemcpyHostToDevice));
  const double bytes = double(n) * NK * 8 + 2.0 * ne * ND * V * 8;
  const int bs = 256;
  const auto ge = (ne + bs - 1) / bs, gp = (n + bs - 1) / bs;
  std::vector<double> r[6];
  auto get = [&](int k) {
    r[k].resize(ne * ND * V);
    CK(cudaMemcpy(r[k].data(), df, ne * ND * V * 8, cudaMemcpyDeviceToHost));
  };
  const float t0 = timeit([&] { per_element<false><<<ge, bs>>>(dKa, ddN, du, df, ne, ld); });
  get(0);
  const float t1 = timeit([&] { per_element<true><<<ge, bs>>>(dKs, ddN, du, df, ne, ld); });
  get(1);
#if (32 % NQV) == 0
  const float t2 = timeit([&] { per_point<false><<<gp, bs>>>(dKa, ddN, du, df, ne, ld); });
  get(2);
  const float t3 = timeit([&] { per_point<true><<<gp, bs>>>(dKs, ddN, du, df, ne, ld); });
  get(3);
#else
  const float t2 = 0, t3 = 0;
  r[2] = r[0]; r[3] = r[0];
#endif
  const auto gs = (ne + EPB - 1) / EPB;
  const float t4 = timeit([&] { per_point_shared<false><<<gs, 256>>>(dKa, ddN, du, df, ne, ld); });
  get(4);
  const float t5 = timeit([&] { per_point_shared<true><<<gs, 256>>>(dKs, ddN, du, df, ne, ld); });
  get(5);
  auto same = [&](int a, int b) {
    return std::memcmp(r[a].data(), r[b].data(), r[a].size() * 8) == 0;
  };
  std::printf("matrix free action, %lld elements, %lld points, K %d comps\n", ne, n, NK);
  std::printf("  thread/element (current)  AoS %.3f ms (%4.0f GB/s)  SoA %.3f ms (%4.0f GB/s)  bitwise AoS==SoA: %s\n",
              t0, bytes / t0 * 1e-6, t1, bytes / t1 * 1e-6, same(0, 1) ? "yes" : "NO");
  std::printf("  thread/point + reduction  AoS %.3f ms (%4.0f GB/s)  SoA %.3f ms (%4.0f GB/s)  bitwise AoS==SoA: %s\n",
              t2, bytes / t2 * 1e-6, t3, bytes / t3 * 1e-6, same(2, 3) ? "yes" : "NO");
  std::printf("  thread/point + shared     AoS %.3f ms (%4.0f GB/s)  SoA %.3f ms (%4.0f GB/s)  bitwise AoS==SoA: %s\n",
              t4, bytes / t4 * 1e-6, t5, bytes / t5 * 1e-6, same(4, 5) ? "yes" : "NO");
  double md = 0;
  for (size_t i = 0; i != r[0].size(); ++i)
    md = std::max(md, std::abs(r[0][i] - r[4][i]) / (std::abs(r[0][i]) + 1e-300));
  std::printf("  max relative difference between the two mappings: %.2e\n", md);
}
