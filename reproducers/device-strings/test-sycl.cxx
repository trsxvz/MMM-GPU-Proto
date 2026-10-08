#include <sycl/sycl.hpp>
#include <cstdio>
#include "model.hxx"
int main() {
  auto q = sycl::queue{};
  auto* r = sycl::malloc_shared<double>(4, q);
  q.parallel_for(sycl::range<1>{4}, [=](sycl::id<1> i) { r[i] = Model<double>{}.execute(i[0]); }).wait();
  std::printf("%g %g %g %g (%s)\n", r[0], r[1], r[2], r[3], q.get_device().get_info<sycl::info::device::name>().c_str());
  return (r[3] == 6 * 3 + 1) ? 0 : 1;
}
