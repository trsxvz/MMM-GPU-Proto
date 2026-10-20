/*!
 * \file   mfront/src/GenericParallel/Behaviour/OpenACCBackend.cxx
 * \brief  This file implements the backend of the `generic-parallel`
 * behaviour interface based on OpenACC
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
#include "MFront/GenericParallel/Behaviour/OpenACCBackend.hxx"

namespace mfront::generic_parallel::behaviour {

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

  void OpenACCBackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include \"MFront/GenericParallelBehaviour/OpenACC.hxx\"\n";
  }  // end of writeSpecificIncludes

  std::string OpenACCBackend::getSharedDataCopy(const std::string_view getter,
                                                const std::string_view) const {
    // copy in the memory of the current device of OpenACC
    return "mfront::gpb::openacc::getSharedDataCopy<" + std::string{getter} +
           ">()";
  }  // end of getSharedDataCopy

  std::string OpenACCBackend::getApply(const std::string_view status) const {
    auto apply = "mfront::gpb::openacc::makeApply(" + std::string{status};
    if (this->device_stack_size != 0) {
      apply += ", " + std::to_string(this->device_stack_size);
    }
    return apply + ")";
  }  // end of getApply

  std::string OpenACCBackend::getMemorySpace() const {
    return "mfront::generic_parallel::getOpenACCMemorySpace()";
  }  // end of getMemorySpace

  OpenACCBackend::~OpenACCBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
