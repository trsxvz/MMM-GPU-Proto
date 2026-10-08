#include <cstdio>
#include "model.hxx"
int main() {
  double r[4] = {0, 1, 2, 3};
#pragma omp target teams distribute parallel for map(tofrom : r[0:4])
  for (int i = 0; i < 4; ++i) {
    r[i] = Model<double>{}.execute(r[i]);
  }
  std::printf("%g %g %g %g\n", r[0], r[1], r[2], r[3]);
  return (r[3] == 6 * 3 + 1) ? 0 : 1;
}
