/*!
 * \file   mfront/src/GenericParallel/Behaviour/SYCLBackend.cxx
 * \brief  This file implements the SYCL backend of the `generic-parallel`
 * behaviour interface
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <ostream>
#include "MFront/GenericParallel/BackendOptions.hxx"
#include "MFront/GenericParallel/Behaviour/SYCLBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  SYCLBackend::SYCLBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    declareDataTransferPolicyOption(validator);
    declareDeviceStackSizeOption(validator);
    validator.validate(opts);
    this->automatic_data_transfers = areDataTransfersAutomatic(opts);
    this->device_stack_size = getDeviceStackSize(opts);
  }  // end of SYCLBackend

  std::string SYCLBackend::getName() const { return "sycl"; }

  std::string SYCLBackend::getSourceFileExtension() const { return "cxx"; }

  void SYCLBackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include \"MFront/GenericParallelBehaviour/SYCL.hxx\"\n";
  }  // end of writeSpecificIncludes

  std::string SYCLBackend::getApply(const std::string_view status) const {
    auto apply = "mfront::gpb::sycl::makeApply(" + std::string{status};
    if (this->device_stack_size != 0) {
      apply += ", " + std::to_string(this->device_stack_size);
    }
    return apply + ")";
  }  // end of getApply

  std::string SYCLBackend::getMemorySpace() const {
    // with automatic transfers, the fields which are not accessible by the
    // device are copied
    if (this->automatic_data_transfers) {
      return "mfront::generic_parallel::any_memory_space";
    }
    return "mfront::generic_parallel::device_memory_space";
  }  // end of getMemorySpace

  std::string SYCLBackend::getExecutionPolicyType() const {
    return "::mfront::generic_parallel::SyclBackendExecutionPolicy";
  }  // end of getExecutionPolicyType

  std::string SYCLBackend::getSharedDataCopy(
      const std::string_view getter, const std::string_view execution) const {
    return "mfront::gpb::sycl::getSharedDataCopy<" + std::string{getter} +
           ">(" + std::string{execution} + ")";
  }  // end of getSharedDataCopy

  void SYCLBackend::writeCall(std::ostream& os,
                              const std::string_view f,
                              const std::string_view template_arguments,
                              const std::string_view status,
                              const std::string_view arguments) const {
    // the functions of the SYCL backend transfer the fields which are not
    // accessible by the device before calling the generic ones
    os << "return mfront::gpb::sycl::" << f << "<" << template_arguments << ">("
       << status << ", " << arguments
       << ", mfront::gpb::sycl::DataTransferPolicy::"
       << (this->automatic_data_transfers ? "AUTOMATIC" : "NONE") << ", "
       << this->getApply(status) << ");\n";
  }  // end of writeCall

  SYCLBackend::~SYCLBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
