/*!
 * \file   mfront/src/GenericParallel/OpenACCBackend.cxx
 * \brief  This file implements the backend of the `generic-parallel`
 * interfaces based on OpenACC
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <ostream>
#include "MFront/GenericParallel/BackendOptions.hxx"
#include "MFront/GenericParallel/OpenACCBackend.hxx"

namespace mfront::generic_parallel {

  OpenACCBackend::OpenACCBackend() = default;

  OpenACCBackend::OpenACCBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    declareDeviceStackSizeOption(validator);
    validator.validate(opts);
    this->device_stack_size = getDeviceStackSize(opts);
  }  // end of OpenACCBackend

  std::string OpenACCBackend::getName() const { return "openacc"; }

  std::string OpenACCBackend::getSourceFileExtension() const {
    return "cxx";
  }  // end of getSourceFileExtension

  void OpenACCBackend::writeSpecificIncludes(std::ostream& os,
                                             const Runtime r) const {
    os << "#include \"MFront/GenericParallel/OpenACC.hxx\"\n";
    if (r == Runtime::BEHAVIOURS) {
      // copies of the shared data
      os << "#include \"MFront/GenericParallelBehaviour/OpenACC.hxx\"\n";
    }
  }  // end of writeSpecificIncludes

  std::string OpenACCBackend::getSharedDataCopy(const std::string_view getter,
                                                const std::string_view) const {
    // copy in the memory of the current device of OpenACC
    return "mfront::gpb::openacc::getSharedDataCopy<" + std::string{getter} +
           ">()";
  }  // end of getSharedDataCopy

  std::string OpenACCBackend::getApply() const {
    if (this->device_stack_size == 0) {
      return "mfront::generic_parallel::openacc::makeApply()";
    }
    return "mfront::generic_parallel::openacc::makeApply(" +
           std::to_string(this->device_stack_size) + ")";
  }  // end of getApply

  std::string OpenACCBackend::getMemorySpace() const {
    return "mfront::generic_parallel::getOpenACCMemorySpace()";
  }  // end of getMemorySpace

  std::string OpenACCBackend::getExecutionPolicyType() const {
    return "::mfront::generic_parallel::OpenACCBackendExecutionPolicy";
  }  // end of getExecutionPolicyType

  OpenACCBackend::~OpenACCBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel
