/*!
 * \file   DeviceOptionTest.cxx
 * \brief  Test of the `--device` option, which creates the MFEM device in
 * `mfem_mgis::initialize`, and of the `--hypre-on-host` option
 * \date   09/10/2026
 */

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include "mfem/config/config.hpp"
#ifdef MFEM_USE_MPI
#include "HYPRE_config.h"
#endif /* MFEM_USE_MPI */
#include "mfem/general/device.hpp"
#include "mfem/general/mem_manager.hpp"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"

/*!
 * \return true if the memory types of MFEM are the expected ones
 * \param[in] expected_host_memory: expected memory of the host, `host`
 * (standard memory), `umpire` or `any` (not checked)
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
  if (std::strcmp(expected_host_memory, "any") == 0) {
    return true;
  }
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

//! \return true if hypre may compute on a GPU
static constexpr bool hypreMayUseGPU() {
#if defined(MFEM_USE_MPI) && defined(HYPRE_USING_GPU) && \
    (MFEM_HYPRE_VERSION >= 23100)
  return true;
#else
  return false;
#endif
}  // end of hypreMayUseGPU

/*!
 * \return the value of the option `--expected-error`, empty if the option
 * is not given. The arguments are read before the initialization of
 * mfem-mgis, which may fail.
 * \param[in] argc: number of arguments
 * \param[in] argv: arguments
 */
static std::string getExpectedError(const int argc, char** const argv) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (std::strcmp(argv[i], "--expected-error") == 0) {
      return argv[i + 1];
    }
  }
  return {};
}  // end of getExpectedError

int main(int argc, char** argv) {
  // an error of the initialization is expected if the option
  // `--expected-error` gives a part of its message
  const auto expected_error = getExpectedError(argc, argv);
  try {
    mfem_mgis::initialize(argc, argv);
  } catch (std::exception& e) {
    std::cout << "initialize: " << e.what() << '\n';
    mfem_mgis::finalize();
    const auto found = (!expected_error.empty()) &&
                       (std::string_view{e.what()}.find(expected_error) !=
                        std::string_view::npos);
    return found ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  if (!expected_error.empty()) {
    std::cerr << "the initialization did not fail\n";
    return EXIT_FAILURE;
  }
  // `none` if no device is expected, `debug` if the debug device of MFEM is
  // expected, `cpu` for a device on the host, `gpu` for a device on a GPU
  const char* expected = "none";
  // expected memory of the host for the `cpu` and `gpu` devices: `host` or
  // `umpire`
  const char* expected_host_memory = "host";
  // where hypre is expected to compute: `host`, `device` (if hypre may
  // compute on a GPU, `host` otherwise) or `any`
  const char* expected_hypre = "any";
  // configuration of a device created by the program after the
  // initialization, empty if none
  const char* created_device = "";
  auto args = mfem::OptionsParser(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&expected, "-e", "--expected",
                 "expected device: none, debug, cpu or gpu");
  args.AddOption(&expected_host_memory, "-ehm", "--expected-host-memory",
                 "expected memory of the host: host or umpire");
  args.AddOption(&expected_hypre, "-eh", "--expected-hypre",
                 "where hypre is expected to compute: host, device or any");
  args.AddOption(&created_device, "-cd", "--create-device",
                 "configuration of a device created by the program after the "
                 "initialization of mfem-mgis");
  // read before the initialization
  const char* expected_error_option = "";
  args.AddOption(&expected_error_option, "-ee", "--expected-error",
                 "part of the message of the expected error of the "
                 "initialization");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  auto program_device = std::unique_ptr<mfem::Device>{};
  if (*created_device != '\0') {
    program_device = std::make_unique<mfem::Device>(created_device);
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
  const auto hypre_on_gpu = mfem_mgis::isHypreUsingGPU();
  if (std::strcmp(expected_hypre, "host") == 0) {
    success = success && !hypre_on_gpu;
  } else if (std::strcmp(expected_hypre, "device") == 0) {
    success = success && (hypre_on_gpu == hypreMayUseGPU());
  } else if (std::strcmp(expected_hypre, "any") != 0) {
    std::cerr << "invalid expected location of hypre '" << expected_hypre
              << "'\n";
    success = false;
  }
  std::cout << "device configured: " << configured
            << ", debug device: " << debug << ", GPU: " << gpu
            << ", hypre on the GPU: " << hypre_on_gpu << '\n';
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
