// nvc++ 25.11, -stdpar=gpu: a member array of doubles with a default member
// initializer `= {}` makes the device compiler emit an invalid LLVM constant
// ("parse integer constant must have integer type") when the energies are
// written through a pointer.
#include <algorithm>
#include <cstdio>
#include <execution>
#include <vector>
struct Values {
  double values[2];
};
struct View {
  explicit View(const double x) : v0{{x, 2 * x}} {}
  double* pointer() { return this->v1.values; }
  Values v0;
  Values v1 = {};
};
int main() {
  auto v = std::vector<double>(1000, 1.);
  auto* const p = v.data();
  std::for_each(std::execution::par_unseq, p, p + v.size(), [](double& x) {
    auto view = View{x};
    double* e = view.pointer();
    e[1] += view.v0.values[1];
    x = e[0] + e[1];
  });
  std::printf("%g (expected 2)\n", v[0]);
  return v[0] == 2 ? 0 : 1;
}
