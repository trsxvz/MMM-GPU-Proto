/*!
 * \file   DeviceForallTest.cxx
 * \brief  Test of a kernel of MFEM (`mfem::forall`) in a source of mfem-mgis
 * compiled as CUDA code when MFEM is compiled with CUDA, on the device
 * selected by the `--device` option
 * \date   09/10/2026
 */

#include <cstdlib>
#include <cstring>
#include <iostream>
#include "mfem/general/device.hpp"
#include "mfem/general/forall.hpp"
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/vector.hpp"
#include "MGIS/Profiling.hxx"
#include "MFEMMGIS/Config.hxx"

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  // `cpu` if the kernel is expected to run on the host, `cuda` if it is
  // expected to run on a CUDA device
  const char* expected = "cpu";
  auto args = mfem::OptionsParser(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&expected, "-e", "--expected",
                 "expected backend of the kernel: cpu or cuda");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  const auto cuda = mfem::Device::Allows(mfem::Backend::CUDA_MASK);
  if ((std::strcmp(expected, "cuda") == 0) != cuda) {
    std::cerr << "unexpected backend\n";
    return EXIT_FAILURE;
  }
  // with a GPU, the profiling sections wait for the kernels
  if ((mgis::getProfilingSynchronization() != nullptr) != cuda) {
    std::cerr << "unexpected synchronization of the profiling sections\n";
    return EXIT_FAILURE;
  }
  constexpr int n = 1000;
  auto v = mfem::Vector(n);
  v.UseDevice(true);
  auto* const d = v.Write();
  mfem::forall(n, [=] MFEM_HOST_DEVICE(int i) { d[i] = 2 * i + 1; });
  // the values are written by the kernel where the vector lives (on the
  // device for a CUDA device), then copied to the host
  const auto* const h = v.HostRead();
  for (int i = 0; i != n; ++i) {
    if (h[i] != 2 * i + 1) {
      std::cerr << "invalid value at " << i << ": " << h[i] << '\n';
      return EXIT_FAILURE;
    }
  }
  // with a CUDA device, the values have been computed on the device
  if (cuda && !v.GetMemory().DeviceIsValid()) {
    std::cerr << "the vector is not valid on the device\n";
    return EXIT_FAILURE;
  }
  std::cout << "kernel run on the " << (cuda ? "CUDA device" : "host")
            << '\n';
  return EXIT_SUCCESS;
}  // end of main
