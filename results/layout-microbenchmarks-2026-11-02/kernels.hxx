// Micro benchmark: layout of the data at the points (AoS / SoA), runtime or
// compile-time strides, thread per point or per element.
// y(o, r) = sum_c K(o, r * G + c) * g(o, c)
#pragma once
#include <cstddef>
#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

using idx_t = long long;

// runtime strides, as the batch kernels of mfem-mgis
struct RuntimeAccess {
  idx_t psK, csK, psG, csG;
  HD idx_t k(idx_t o, int c) const { return o * psK + c * csK; }
  HD idx_t g(idx_t o, int c) const { return o * psG + c * csG; }
};

// compile time layout
template <bool soa, int G>
struct StaticAccess {
  idx_t n; // leading dimension in SoA
  HD idx_t k(idx_t o, int c) const {
    if constexpr (soa) {
      return c * n + o;
    } else {
      return o * (G * G) + c;
    }
  }
  HD idx_t g(idx_t o, int c) const {
    if constexpr (soa) {
      return c * n + o;
    } else {
      return o * G + c;
    }
  }
};

template <int G, typename A>
HD inline void point(const double* K, const double* g, double* y,
                     const A& a, idx_t o) {
  double gv[G];
  for (int c = 0; c != G; ++c) {
    gv[c] = g[a.g(o, c)];
  }
  for (int r = 0; r != G; ++r) {
    double s = 0;
    for (int c = 0; c != G; ++c) {
      s += K[a.k(o, r * G + c)] * gv[c];
    }
    y[a.g(o, r)] = s;
  }
}
