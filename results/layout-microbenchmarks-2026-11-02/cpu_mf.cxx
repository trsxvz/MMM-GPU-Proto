// CPU version of the matrix free action (one iteration per element, as
// mfem-mgis on the host), AoS and SoA layouts of K, compile time.
#include <chrono>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <omp.h>
#ifndef NQV
#define NQV 8
#endif
#ifndef GV
#define GV 6
#endif
constexpr int ND = 8, NQ = NQV, V = 3, G = GV, NK = G * G;
using idx_t = long long;
template <bool soa>
inline idx_t kidx(idx_t o, int c, idx_t ld) {
  if constexpr (soa) return c * ld + o; else return o * NK + c;
}
template <bool soa>
inline void point(double* fl, const double* K, const double* dN,
                  const double* ue, double w, idx_t o, int q, idx_t ld) {
  double gu[3][3] = {};
  for (int a = 0; a != ND; ++a)
    for (int j = 0; j != 3; ++j) {
      const double d = dN[(q * ND + a) * 3 + j];
      for (int i = 0; i != 3; ++i) gu[i][j] += ue[a + ND * i] * d;
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
double run(const double* K, const double* dN, const double* u, double* f,
           idx_t ne, idx_t ld) {
  std::vector<double> t;
  for (int i = 0; i != 11; ++i) {
    const auto b = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
    for (idx_t e = 0; e < ne; ++e) {
      double fl[ND * V] = {};
      for (int q = 0; q != NQ; ++q) point<soa>(fl, K, dN, u + e * ND * V, 0.125, e * NQ + q, q, ld);
      for (int k = 0; k != ND * V; ++k) f[e * ND * V + k] = fl[k];
    }
    const auto e = std::chrono::steady_clock::now();
    if (i >= 3) t.push_back(std::chrono::duration<double, std::milli>(e - b).count());
  }
  std::sort(t.begin(), t.end());
  return t[t.size() / 2];
}
void bench(idx_t n0) {
  const idx_t ne = n0 / NQ, n = ne * NQ, ld = n;
  std::vector<double> Ka(n * NK), Ks(ld * NK), u(ne * ND * V), dN(NQ * ND * 3), fa(ne * ND * V), fs(ne * ND * V);
#pragma omp parallel for schedule(static)
  for (idx_t e = 0; e < ne; ++e)
    for (int q = 0; q != NQ; ++q) {
      const idx_t o = e * NQ + q;
      for (int c = 0; c != NK; ++c) {
        const double v = 1.0 + 1e-3 * ((o * 7 + c * 13) % 101);
        Ka[o * NK + c] = v;
        Ks[c * ld + o] = v;
      }
    }
  for (idx_t i = 0; i != ne * ND * V; ++i) u[i] = 1e-3 * (i % 97);
  for (int i = 0; i != NQ * ND * 3; ++i) dN[i] = 0.1 * ((i % 7) - 3);
  const double bytes = double(n) * NK * 8 + 2.0 * ne * ND * V * 8;
  const double ta = run<false>(Ka.data(), dN.data(), u.data(), fa.data(), ne, ld);
  const double ts = run<true>(Ks.data(), dN.data(), u.data(), fs.data(), ne, ld);
  std::printf("  n=%lld K %d comps nq=%d threads=%d: AoS %.2f ms (%.1f GB/s) SoA %.2f ms (%.1f GB/s) ratio SoA/AoS %.2f bitwise %s\n",
              n, NK, NQ, omp_get_max_threads(), ta, bytes / ta * 1e-6, ts, bytes / ts * 1e-6, ts / ta,
              std::memcmp(fa.data(), fs.data(), fa.size() * 8) == 0 ? "yes" : "NO");
}
int main() {
  bench(2000016);
  bench(idx_t{1} << 21);
}
