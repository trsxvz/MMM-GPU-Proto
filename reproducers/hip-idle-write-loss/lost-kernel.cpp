// The writes of a kernel are lost when the next kernel is the first kernel
// of libmfem launched by the program, whose code object is then loaded
// (about 2 ms): the GPU stays idle long enough to drop them (idle-loss.cpp,
// ROCm 10.1 on the integrated GPU gfx1036). The pending kernel cancels the
// entries of a vector whose indices are given by an array, as
// mfem::ParNonlinearForm::Mult does for the essential degrees of freedom.
//
// usage: lost-kernel [warm-up]
// with any argument, a kernel of libmfem is launched once before: nothing is
// then lost. AMD_SERIALIZE_KERNEL=3 AMD_SERIALIZE_COPY=3 also avoids it.
#include <cstdio>
#include <cstdlib>
#include "mfem.hpp"

int main(int argc, char** argv) {
  mfem::Device device("hip");
  if (argc > 1) {
    mfem::Vector w(16);
    w.UseDevice(true);
    w = 1.;  // first kernel of libmfem, nothing pending
  }
  constexpr auto n = 81;
  mfem::Vector y(n);
  auto* const h = y.HostWrite();
  for (int i = 0; i != n; ++i) {
    h[i] = 1.;
  }
  mfem::Array<int> indices;
  for (int i = 0; i < n; i += 3) {
    indices.Append(i);
  }
  // counter of the threads of the cancellation kernel which have run
  mfem::Array<int> count(1);
  count.HostWrite()[0] = 0;
  auto* const c = count.ReadWrite();
  const auto idx = indices.Read();
  auto* const yd = y.ReadWrite();
  mfem::forall(indices.Size(), [=] MFEM_HOST_DEVICE(int i) {
    yd[idx[i]] = 0.;
    AtomicAdd(c[0], 1);
  });
  // first launch of a kernel of libmfem while the cancellation is pending
  mfem::Vector other(16);
  other.UseDevice(true);
  other = 3.;
  const auto* const hy = y.HostRead();
  const auto ran = count.HostRead()[0];
  std::printf("threads of the cancellation kernel which ran: %d / %d, "
              "y[0] = %g (expected 0)\n", ran, indices.Size(), hy[0]);
  return (hy[0] == 0.) ? EXIT_SUCCESS : EXIT_FAILURE;
}
