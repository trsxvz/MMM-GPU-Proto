/*!
 * \file   mfront/src/GenericParallel/Behaviour/RAJABackend.cxx
 * \brief  This file implements the RAJA backend of the `generic-parallel`
 * behaviour interface
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
#include "MFront/GenericParallel/Behaviour/RAJABackend.hxx"

namespace mfront::generic_parallel::behaviour {

  RAJABackend::RAJABackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    declareRAJAExecutionPolicyOption(validator);
    validator.validate(opts);
    this->execution_policy = getRAJAExecutionPolicy(opts);
  }  // end of RAJABackend

  std::string RAJABackend::getName() const { return "raja"; }

  std::string RAJABackend::getSourceFileExtension() const { return "cxx"; }

  void RAJABackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include \"MFront/GenericParallelBehaviour/RAJA.hxx\"\n";
  }  // end of writeSpecificIncludes

  std::string RAJABackend::getApply(const std::string_view status) const {
    return "mfront::gpb::raja::makeApply<" + this->execution_policy + ">(" +
           std::string{status} + ")";
  }  // end of getApply

  std::string RAJABackend::getMemorySpace() const {
    // the policies of CUDA and HIP may not access the memory of the host
    return "(mfront::generic_parallel::raja::isDeviceExecutionPolicy<" +
           this->execution_policy +
           "> ? mfront::generic_parallel::device_memory_space : "
           "mfront::generic_parallel::host_memory_space)";
  }  // end of getMemorySpace

  std::string RAJABackend::getSharedDataCopy(const std::string_view getter,
                                             const std::string_view) const {
    return "mfront::gpb::raja::getSharedDataCopy<" + this->execution_policy +
           ", " + std::string{getter} + ">()";
  }  // end of getSharedDataCopy

  RAJABackend::~RAJABackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
