/*!
 * \file   mfront/src/GenericParallel/Behaviour/ParallelSTLBackend.cxx
 * \brief
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
#include "MFront/GenericParallel/Behaviour/ParallelSTLBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  ParallelSTLBackend::ParallelSTLBackend() = default;

  ParallelSTLBackend::ParallelSTLBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    declareExecutionPolicyOption(validator);
    declareDeviceStackSizeOption(validator);
    validator.validate(opts);
    this->execution_policy = getExecutionPolicy(opts);
    this->device_stack_size = getDeviceStackSize(opts);
  }  // end of ParallelSTLBackend

  std::string ParallelSTLBackend::getName() const { return "stlpar"; }

  std::string ParallelSTLBackend::getSourceFileExtension() const {
    return "cxx";
  }  // end of getSourceFileExtension

  void ParallelSTLBackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include <execution>\n"
       << "#include \"MFront/GenericParallelBehaviour/ParallelSTL.hxx\"\n";
  }  // end of writeSpecificIncludes

  std::string ParallelSTLBackend::getApply(
      const std::string_view status) const {
    auto apply = "mfront::gpb::stlpar::makeApply(" + std::string{status} +
                 ", " + this->execution_policy;
    if (this->device_stack_size != 0) {
      apply += ", " + std::to_string(this->device_stack_size);
    }
    return apply + ")";
  }  // end of getApply

  ParallelSTLBackend::~ParallelSTLBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
