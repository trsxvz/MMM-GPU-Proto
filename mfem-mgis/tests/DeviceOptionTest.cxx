/*!
 * \file   DeviceOptionTest.cxx
 * \brief  Test of the `--device` option, which creates the MFEM device in
 * `mfem_mgis::initialize`
 * \date   09/10/2026
 */

#include <cstdlib>
#include <cstring>
#include <iostream>
#include "mfem/general/device.hpp"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  // `none` if no device is expected, `debug` if the debug device of MFEM is
  // expected
  const char* expected = "none";
  auto args = mfem::OptionsParser(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&expected, "-e", "--expected",
                 "expected device: none or debug");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  const auto configured = mfem::Device::IsConfigured();
  const auto debug = mfem::Device::Allows(mfem::Backend::DEBUG_DEVICE);
  auto success = false;
  if (std::strcmp(expected, "none") == 0) {
    success = !configured;
  } else if (std::strcmp(expected, "debug") == 0) {
    success = configured && debug;
  } else {
    std::cerr << "invalid expected device '" << expected << "'\n";
  }
  std::cout << "device configured: " << configured
            << ", debug device: " << debug << '\n';
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
