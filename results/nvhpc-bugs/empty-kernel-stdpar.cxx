// nvc++ -O2 -stdpar=gpu bug.cxx && ./a.out  -> expected 23.6
#include <algorithm>
#include <cstdio>
#include <execution>
#include <numeric>
#include <vector>

int main() {
  const std::size_t n = 1000;
  std::vector<std::size_t> idx(n);
  std::iota(idx.begin(), idx.end(), std::size_t{0});
  std::vector<double> out(n, 0.);
  double* const pout = out.data();
  std::for_each(std::execution::par_unseq, idx.begin(), idx.end(),
                [pout](const std::size_t i) {
                  double D[36];
                  for (int k = 0; k != 36; ++k) D[k] = 2.;
                  double s = D[0];
                  for (int k = 0; k != 36; ++k) s += D[k] * 0.3;
                  pout[i] = s;
                });
  std::printf("%g (expected 23.6)\n", out[0]);
}
