// nvc++ 25.11, -mp=gpu, -O2 and -O3: internal compiler error "Computed goto
// control expression must be a 32-bit integer" on the product of a tensor by
// a symmetric tensor of TFEL in an OpenMP target region. The switch
// statement of the tensor view of the symmetric tensor is inlined in the
// switch statement of the product. Not minimal: requires headers of TFEL
// where TensorViewFromStensor::operator() is a host-device function using a
// switch statement (before the workaround).
// nvc++ -O2 -mp=gpu -std=c++20 -DDIM=2 -I <tfel>/include omp-tensor-view-switch.cxx
#include <cstdio>
#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
int main() {
  double r[5] = {-1, -1, -1, -1, -1};
#pragma omp target teams distribute parallel for map(tofrom : r)
  for (int k = 0; k < 1; ++k) {
    const auto F = tfel::math::tensor<DIM, double>::Id();
    auto S = tfel::math::stensor<DIM, double>{};
    S(0) = 1;
    S(1) = 2;
    S(2) = 3;
    tfel::math::tensor<DIM, double> P;
    P = F * S;
    for (unsigned short i = 0; i != 5; ++i) {
      r[i] = P(i);
    }
  }
  std::printf("%g %g %g %g %g (expected 1 2 3 0 0)\n", r[0], r[1], r[2], r[3],
              r[4]);
}
