/*!
 * \file   HypreDevicePoolSizeTest.cxx
 * \brief  Test of the size of the pool of memory of hypre on the GPU, set by
 * `mfem_mgis::initialize` from the option `--hypre-device-pool-size` (or
 * `-hdps`) or from the environment variable
 * `MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE`. The option `--expected` gives the
 * expected size in GiB, or `error` if the initialization must fail. No
 * memory is allocated on the GPU.
 * \date   11/10/2026
 */

#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <exception>
#include "mfem/config/config.hpp"
#include "mfem/general/optparser.hpp"
#ifdef MFEM_USE_MPI
#include "_hypre_utilities.h"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Config.hxx"

int main(int argc, char** argv) {
  // the expected result is read before the initialization, which may fail
  auto expected = std::string{};
  for (auto i = 1; i + 1 < argc; ++i) {
    if (std::strcmp(argv[i], "--expected") == 0) {
      expected = argv[i + 1];
    }
  }
  if (expected.empty()) {
    std::cerr << "no expected result given (--expected)\n";
    return EXIT_FAILURE;
  }
  try {
    mfem_mgis::initialize(argc, argv);
  } catch (std::exception& e) {
    std::cout << "initialize: " << e.what() << '\n';
    mfem_mgis::finalize();
    if (expected == "error") {
      return EXIT_SUCCESS;
    }
#if defined(MFEM_USE_MPI) && defined(HYPRE_USING_UMPIRE_DEVICE)
    return EXIT_FAILURE;
#else
    // the option is refused if hypre does not use Umpire on the GPU
    return EXIT_SUCCESS;
#endif
  }
  auto args = mfem::OptionsParser(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  const char* e = "";
  args.AddOption(&e, "-e", "--expected",
                 "expected size in GiB, or error if the initialization must "
                 "fail");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    mfem_mgis::finalize();
    return EXIT_FAILURE;
  }
  auto success = expected != "error";
  if (!success) {
    std::cerr << "the initialization did not fail\n";
  }
#if defined(MFEM_USE_MPI) && defined(HYPRE_USING_UMPIRE_DEVICE)
  if (success) {
    const auto size = hypre_HandleUmpireDevicePoolSize(hypre_handle());
    const auto value =
        static_cast<std::size_t>(std::stod(expected) * 1073741824.);
    std::cout << "size of the pool of memory of hypre on the GPU: " << size
              << " bytes\n";
    if (size != value) {
      std::cerr << "expected " << value << " bytes\n";
      success = false;
    }
  }
#endif
  mfem_mgis::finalize();
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
