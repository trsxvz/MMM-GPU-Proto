#include <algorithm>
#include <execution>
#include <vector>
#include <cstdio>
#include "model.hxx"
int main() {
  auto r = std::vector<double>{0, 1, 2, 3};
  std::transform(std::execution::par_unseq, r.begin(), r.end(), r.begin(),
                 [](const double x) { return Model<double>{}.execute(x); });
  std::printf("%g %g %g %g\n", r[0], r[1], r[2], r[3]);
  return (r[3] == 6 * 3 + 1) ? 0 : 1;
}
