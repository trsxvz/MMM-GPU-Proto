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

// L: 0 AoS static, 1 SoA static, 2 runtime strides
struct Strides { idx_t ps, cs; };
__constant__ Strides rs;
template <int L>
__device__ inline idx_t kidx(idx_t o, int c, idx_t ld) {
  if constexpr (L == 1) {
    return c * ld + o;
  } else if constexpr (L == 0) {
    return o * NK + c;
  } else {
    return o * rs.ps + c * rs.cs;
  }
}

// contribution of point q of element e: fl[node * V + i]
template <int soa>
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

template <int soa>
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

template <int soa>
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
template <int soa>
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
#ifndef NEV
  const idx_t ne = 2000016 / NQ;
#else
  const idx_t ne = NEV;
#endif
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
  std::vector<double> r[8];
  auto get = [&](int k) {
    r[k].resize(ne * ND * V);
    CK(cudaMemcpy(r[k].data(), df, ne * ND * V * 8, cudaMemcpyDeviceToHost));
  };
  auto setrs = [](idx_t ps, idx_t cs) { Strides h{ps, cs}; CK(cudaMemcpyToSymbol(rs, &h, sizeof(h))); };
  const auto gs = (ne + EPB - 1) / EPB;
  float t[2][2][2];
  t[0][0][0] = timeit([&] { per_element<0><<<ge, bs>>>(dKa, ddN, du, df, ne, ld); }); get(0);
  t[0][1][0] = timeit([&] { per_element<1><<<ge, bs>>>(dKs, ddN, du, df, ne, ld); }); get(1);
  t[1][0][0] = timeit([&] { per_point_shared<0><<<gs, 256>>>(dKa, ddN, du, df, ne, ld); }); get(2);
  t[1][1][0] = timeit([&] { per_point_shared<1><<<gs, 256>>>(dKs, ddN, du, df, ne, ld); }); get(3);
  setrs(NK, 1);
  t[0][0][1] = timeit([&] { per_element<2><<<ge, bs>>>(dKa, ddN, du, df, ne, ld); }); get(4);
  t[1][0][1] = timeit([&] { per_point_shared<2><<<gs, 256>>>(dKa, ddN, du, df, ne, ld); }); get(5);
  setrs(1, ld);
  t[0][1][1] = timeit([&] { per_element<2><<<ge, bs>>>(dKs, ddN, du, df, ne, ld); }); get(6);
  t[1][1][1] = timeit([&] { per_point_shared<2><<<gs, 256>>>(dKs, ddN, du, df, ne, ld); }); get(7);
  bool all = true;
  for (int k = 1; k != 8; ++k) all = all && (std::memcmp(r[0].data(), r[k].data(), r[0].size() * 8) == 0);
  std::printf("matrix free action, %lld elements, %lld points, K %d comps, nq %d; 8 variants bitwise identical: %s\n", ne, n, NK, NQ, all ? "yes" : "NO");
  const char* mp[2] = {"thread/element     ", "thread/point+shared"};
  for (int m = 0; m != 2; ++m)
    std::printf("  %s AoS static %4.0f GB/s runtime %4.0f GB/s (%+5.1f %%) | SoA static %4.0f GB/s runtime %4.0f GB/s (%+5.1f %%)\n", mp[m],
                bytes / t[m][0][0] * 1e-6, bytes / t[m][0][1] * 1e-6, 100 * (t[m][0][1] / t[m][0][0] - 1),
                bytes / t[m][1][0] * 1e-6, bytes / t[m][1][1] * 1e-6, 100 * (t[m][1][1] / t[m][1][0] - 1));
  return 0;
  double md = 0;
  for (size_t i = 0; i != r[0].size(); ++i)
    md = std::max(md, std::abs(r[0][i] - r[4][i]) / (std::abs(r[0][i]) + 1e-300));
  std::printf("  max relative difference between the two mappings: %.2e\n", md);
}
