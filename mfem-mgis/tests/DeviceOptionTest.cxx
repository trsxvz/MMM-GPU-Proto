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
#include "mfem/general/mem_manager.hpp"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"

/*!
 * \return true if the memory types of MFEM are the expected ones
 * \param[in] expected_host_memory: expected memory of the host, `host`
 * (standard memory) or `umpire`
 * \param[in] gpu: if true, the device uses a GPU
 */
static bool checkMemoryTypes(const char* const expected_host_memory,
                             const bool gpu) {
  const auto host = mfem::Device::GetHostMemoryType();
  const auto device = mfem::Device::GetDeviceMemoryType();
  std::cout << "memory of the host: " << mfem::MemoryTypeName[int(host)]
            << ", memory of the device: " << mfem::MemoryTypeName[int(device)]
            << '\n';
  auto check = [](const char* const what, const mfem::MemoryType value,
                  const mfem::MemoryType expected) {
    if (value == expected) {
      return true;
    }
    std::cerr << "invalid memory of the " << what << ": expected "
              << mfem::MemoryTypeName[int(expected)] << '\n';
    return false;
  };
  if (std::strcmp(expected_host_memory, "umpire") == 0) {
    return check("host", host, mfem::MemoryType::HOST_UMPIRE);
  }
  if (std::strcmp(expected_host_memory, "host") != 0) {
    std::cerr << "invalid expected memory of the host '" << expected_host_memory
              << "'\n";
    return false;
  }
  if (!check("host", host, mfem::MemoryType::HOST)) {
    return false;
  }
  if (!gpu) {
    // no memory of the device without GPU
    return check("device", device, mfem::MemoryType::HOST);
  }
  // default memory of the device of MFEM
#ifdef MFEM_USE_UMPIRE
  return check("device", device, mfem::MemoryType::DEVICE_UMPIRE);
#else  /* MFEM_USE_UMPIRE */
  return check("device", device, mfem::MemoryType::DEVICE);
#endif /* MFEM_USE_UMPIRE */
}  // end of checkMemoryTypes

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  // `none` if no device is expected, `debug` if the debug device of MFEM is
  // expected, `cpu` for a device on the host, `gpu` for a device on a GPU
  const char* expected = "none";
  // expected memory of the host for the `cpu` and `gpu` devices: `host` or
  // `umpire`
  const char* expected_host_memory = "host";
  auto args = mfem::OptionsParser(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&expected, "-e", "--expected",
                 "expected device: none, debug, cpu or gpu");
  args.AddOption(&expected_host_memory, "-ehm", "--expected-host-memory",
                 "expected memory of the host: host or umpire");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  const auto configured = mfem::Device::IsConfigured();
  const auto debug = mfem::Device::Allows(mfem::Backend::DEBUG_DEVICE);
  const auto gpu =
      mfem::Device::Allows(mfem::Backend::CUDA_MASK | mfem::Backend::HIP_MASK);
  auto success = false;
  if (std::strcmp(expected, "none") == 0) {
    success = !configured;
  } else if (std::strcmp(expected, "debug") == 0) {
    success = configured && debug;
  } else if (std::strcmp(expected, "cpu") == 0) {
    success = configured && !debug && !gpu &&
              checkMemoryTypes(expected_host_memory, false);
  } else if (std::strcmp(expected, "gpu") == 0) {
    success = configured && gpu && checkMemoryTypes(expected_host_memory, true);
  } else {
    std::cerr << "invalid expected device '" << expected << "'\n";
  }
  std::cout << "device configured: " << configured
            << ", debug device: " << debug << ", GPU: " << gpu << '\n';
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
