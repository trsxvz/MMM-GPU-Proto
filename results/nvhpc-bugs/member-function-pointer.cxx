// nvc++ 25.11, -stdpar=gpu: calling a member function through a pointer to
// member function given as a template argument leads to an illegal memory
// access on the GPU, although the pointer is a constant expression.
#include <algorithm>
#include <cstdio>
#include <execution>
#include <vector>
struct Behaviour {
  void twice(double* const v) { v[0] *= 2; }
};
template <void (Behaviour::*m)(double* const)>
void call(double& x) {
  auto b = Behaviour{};
  (b.*m)(&x);
}
int main() {
  auto v = std::vector<double>(1000, 1.);
  auto* const p = v.data();
  std::for_each(std::execution::par_unseq, p, p + v.size(),
                [](double& x) { call<&Behaviour::twice>(x); });
  std::printf("%g (expected 2)\n", v[0]);
  return v[0] == 2 ? 0 : 1;
}
