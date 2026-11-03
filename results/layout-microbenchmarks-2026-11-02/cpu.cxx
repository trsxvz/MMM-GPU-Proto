#include <chrono>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "kernels.hxx"

template <int G, typename A>
double run(bool by_point, const double* K, const double* g, double* y,
           const A a, idx_t n, int nq) {
  std::vector<double> t;
  for (int i = 0; i != 13; ++i) {
    const auto b = std::chrono::steady_clock::now();
    if (by_point) {
#pragma omp parallel for schedule(static)
      for (idx_t o = 0; o < n; ++o) point<G>(K, g, y, a, o);
    } else {
      const idx_t ne = n / nq;
#pragma omp parallel for schedule(static)
      for (idx_t e = 0; e < ne; ++e)
        for (int q = 0; q != nq; ++q) point<G>(K, g, y, a, e * nq + q);
    }
    const auto e = std::chrono::steady_clock::now();
    if (i >= 3) t.push_back(std::chrono::duration<double, std::milli>(e - b).count());
  }
  std::sort(t.begin(), t.end());
  return t[t.size() / 2];
}

template <int G>
void bench(idx_t n, int nq, idx_t ld) {
  constexpr int NK = G * G;
  std::vector<double> Ka(n * NK), ga(n * G), Ks(ld * NK), gs(ld * G), y(ld * G);
  // first touch in parallel, as the kernels
#pragma omp parallel for schedule(static)
  for (idx_t o = 0; o < n; ++o) {
    for (int c = 0; c != NK; ++c) {
      const double v = 1.0 + 1e-3 * ((o * 7 + c * 13) % 101);
      Ka[o * NK + c] = v;
      Ks[c * ld + o] = v;
    }
    for (int c = 0; c != G; ++c) {
      const double v = 1.0 + 1e-2 * ((o * 3 + c * 5) % 37);
      ga[o * G + c] = v;
      gs[c * ld + o] = v;
      y[o * G + c] = 0;
    }
  }
  const double bytes = double(n) * (NK + 2 * G) * 8;
  const RuntimeAccess ra{NK, 1, G, 1};
  const RuntimeAccess rs{1, ld, 1, ld};
  const StaticAccess<false, G> sa{ld};
  const StaticAccess<true, G> ss{ld};
  std::printf("G=%d n=%lld ld=%lld nq=%d threads=%d\n", G, n, ld, nq, omp_get_max_threads());
  for (int bp = 1; bp >= 0; --bp) {
    const char* m = bp ? "loop/point  " : "loop/element";
    const double t1 = run<G>(bp, Ka.data(), ga.data(), y.data(), sa, n, nq);
    const double t2 = run<G>(bp, Ka.data(), ga.data(), y.data(), ra, n, nq);
    const double t3 = run<G>(bp, Ks.data(), gs.data(), y.data(), ss, n, nq);
    const double t4 = run<G>(bp, Ks.data(), gs.data(), y.data(), rs, n, nq);
    std::printf(
        "  %s AoS static %7.2f ms (%4.1f GB/s) runtime %7.2f ms (%4.1f GB/s)"
        " | SoA static %7.2f ms (%4.1f GB/s) runtime %7.2f ms (%4.1f GB/s)\n",
        m, t1, bytes / t1 * 1e-6, t2, bytes / t2 * 1e-6, t3,
        bytes / t3 * 1e-6, t4, bytes / t4 * 1e-6);
  }
}

int main() {
  const idx_t n2 = idx_t{1} << 21;
  const idx_t n = 2000003;  // not a power of two
  bench<6>(n2, 8, n2);
  bench<6>(n, 8, n);
  bench<6>(n, 8, (n + 7) / 8 * 8 + 8);
  bench<9>(n2, 8, n2);
  bench<9>(n, 8, n);
  bench<9>(n, 8, (n + 7) / 8 * 8 + 8);
}
