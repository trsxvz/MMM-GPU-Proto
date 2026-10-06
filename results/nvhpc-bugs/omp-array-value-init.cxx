// nvc++ 25.11, -mp=gpu, -O2: the value-initialization of a member array
// (`v()`, `v{}` or `= {}`) of an object used in an OpenMP target region is
// implemented by the copy of an uninitialized global array (`.common`), and
// the optimizer removes the code of the kernel: the results are not
// written. Correct at -O0, -O1 and -O3, and when the values are set by a
// loop.
#include <cstdio>
struct Array;
struct Scaled {
  double a;
  const Array& x;
};
struct Array {
  double v[9];
  Array() : v() {}
  explicit Array(const Scaled& e) {
    for (int i = 0; i != 9; ++i) {
      v[i] = e.a * e.x.v[i];
    }
  }
};
int main() {
  double in[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  double r[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
#pragma omp target teams distribute parallel for map(to : in) map(tofrom : r)
  for (int k = 0; k < 1; ++k) {
    Array F;
    for (int i = 0; i != 9; ++i) {
      F.v[i] = in[i];
    }
    const Array P(Scaled{2, F});
    for (int i = 0; i != 9; ++i) {
      r[i] = P.v[i];
    }
  }
  for (const auto x : r) {
    std::printf("%g ", x);
  }
  std::printf("(expected 2 4 6 8 10 12 14 16 18)\n");
  return r[0] == 2 ? 0 : 1;
}
