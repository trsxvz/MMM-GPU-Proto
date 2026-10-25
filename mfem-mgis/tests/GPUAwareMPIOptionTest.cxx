/*!
 * \file   GPUAwareMPIOptionTest.cxx
 * \brief  Test of GPU-aware MPI, requested from `mfem_mgis::initialize` by
 * the option `--gpu-aware-mpi` (or `-gam`) or by the environment variable
 * `MFEM_MGIS_GPU_AWARE_MPI` (overridden by the option `--no-gpu-aware-mpi`).
 * The option `--expected` gives the expected state: `enabled`, `disabled`,
 * or `error` if the initialization must fail (no device of MFEM using a GPU,
 * the variable being then ignored with a warning). The test is skipped
 * (exit status 77) if GPU-aware MPI is expected but the MPI library states
 * that it can't read the memory of the GPUs.
 * \date   25/10/2026
 */

#include <string>
#include <string_view>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <exception>
#include "mfem/general/device.hpp"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"

int main(int argc, char** argv) {
  // the expected result is read before the initialization, which may fail
  auto expected = std::string{};
  for (auto i = 1; i + 1 < argc; ++i) {
    if (std::strcmp(argv[i], "--expected") == 0) {
      expected = argv[i + 1];
    }
  }
  if ((expected != "enabled") && (expected != "disabled") &&
      (expected != "error")) {
    std::cerr << "invalid or missing expected result (--expected)\n";
    return EXIT_FAILURE;
  }
  try {
    mfem_mgis::initialize(argc, argv);
  } catch (std::exception& e) {
    std::cout << "initialize: " << e.what() << '\n';
    mfem_mgis::finalize();
    const auto message = std::string_view{e.what()};
    if ((expected == "enabled") &&
        (message.find("not supported by the MPI library") !=
         std::string_view::npos)) {
      std::cout << "skipped\n";
      return 77;
    }
    const auto refused =
        message.find("GPU-aware MPI") != std::string_view::npos;
    return ((expected == "error") && refused) ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  auto args = mfem::OptionsParser(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  const char* e = "";
  args.AddOption(&e, "-e", "--expected",
                 "expected state of GPU-aware MPI: enabled, disabled or error");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    mfem_mgis::finalize();
    return EXIT_FAILURE;
  }
  auto success = true;
  if (expected == "error") {
    std::cerr << "the initialization did not fail\n";
    success = false;
  } else if (mfem::Device::GetGPUAwareMPI() != (expected == "enabled")) {
    std::cerr << "unexpected state of GPU-aware MPI after the "
                 "initialization\n";
    success = false;
  }
  mfem_mgis::finalize();
  std::cout << (success ? "success" : "failure") << '\n';
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
