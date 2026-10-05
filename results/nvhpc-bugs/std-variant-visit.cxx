// nvc++ 25.11, -stdpar=gpu: visiting a `std::variant` (libstdc++) in a
// parallel algorithm. Expected output: 5.
#include <algorithm>
#include <array>
#include <cstdio>
#include <execution>
#include <type_traits>
#include <variant>
#include <vector>
using Value = std::variant<std::monostate, double, std::array<double, 4>>;
int main() {
  auto v = std::vector<double>(1000, 1.);
  auto* const p = v.data();
  std::for_each(std::execution::par_unseq, p, p + v.size(), [](double& x) {
    auto value = Value{};
    if (x > 0) {
      value = std::array<double, 4>{x, 2 * x, 3 * x, 4 * x};
    } else {
      value = x;
    }
    x = std::visit(
        [](const auto& a) -> double {
          using A = std::decay_t<decltype(a)>;
          if constexpr (std::is_same_v<A, std::array<double, 4>>) {
            return a[0] + a[3];
          } else if constexpr (std::is_same_v<A, double>) {
            return a;
          } else {
            return -1;
          }
        },
        value);
  });
  std::printf("%g (expected 5)\n", v[0]);
  return v[0] == 5 ? 0 : 1;
}
