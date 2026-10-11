/*!
 * \file   KokkosOptionTest.cxx
 * \brief  Test of the initialization of Kokkos by `mfem_mgis::initialize`,
 * requested by the option `--kokkos` (or `-kk`) or by the environment
 * variable `MFEM_MGIS_KOKKOS` (overridden by the option `--no-kokkos`), and
 * of its
 * finalization by `mfem_mgis::finalize`. The option `--expected` gives the
 * expected state: `initialized`, `none`, or `error` if the initialization
 * must fail (option given to mfem-mgis built without Kokkos, the variable
 * being then ignored with a warning).
 * \date   11/10/2026
 */

#include <string>
#include <string_view>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <exception>
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
  if ((expected != "initialized") && (expected != "none") &&
      (expected != "error")) {
    std::cerr << "invalid or missing expected result (--expected)\n";
    return EXIT_FAILURE;
  }
  try {
    mfem_mgis::initialize(argc, argv);
  } catch (std::exception& e) {
    std::cout << "initialize: " << e.what() << '\n';
    mfem_mgis::finalize();
    const auto refused = std::string_view{e.what()}.find(
                             "built without Kokkos") != std::string_view::npos;
    return ((expected == "error") && refused) ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  auto args = mfem::OptionsParser(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  const char* e = "";
  args.AddOption(&e, "-e", "--expected",
                 "expected state of Kokkos: initialized, none or error");
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
  } else if (mfem_mgis::isKokkosInitialized() != (expected == "initialized")) {
    std::cerr << "unexpected state of Kokkos after the initialization\n";
    success = false;
  }
  mfem_mgis::finalize();
  if (mfem_mgis::isKokkosInitialized()) {
    std::cerr << "Kokkos is not finalized\n";
    success = false;
  }
  std::cout << (success ? "success" : "failure") << '\n';
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
