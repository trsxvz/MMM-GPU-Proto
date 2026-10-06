/*!
 * \file   mfront/src/GenericParallel/Behaviour/OpenMPBackend.cxx
 * \brief  This file implements the backend of the `generic-parallel`
 * behaviour interface based on OpenMP
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
#include "MFront/GenericParallel/Behaviour/OpenMPBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  OpenMPBackend::OpenMPBackend() = default;

  OpenMPBackend::OpenMPBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addDataTypeValidator<bool>("offload");
    declareDeviceStackSizeOption(validator);
    validator.validate(opts);
    if (contains(opts, "offload")) {
      this->offload = get<bool>(opts, "offload");
    }
    this->device_stack_size = getDeviceStackSize(opts);
  }  // end of OpenMPBackend

  std::string OpenMPBackend::getName() const { return "openmp"; }

  std::string OpenMPBackend::getSourceFileExtension() const {
    return "cxx";
  }  // end of getSourceFileExtension

  void OpenMPBackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include \"MFront/GenericParallelBehaviour/OpenMP.hxx\"\n";
  }  // end of writeSpecificIncludes

  std::string OpenMPBackend::getApply(const std::string_view status) const {
    // the offload option is a template argument, so that the code
    // executed on the host does not depend on the support of offloading by
    // the OpenMP runtime
    auto apply = std::string{"mfront::gpb::openmp::makeApply<"} +
                 (this->offload ? "true" : "false") + ">(" +
                 std::string{status};
    if (this->device_stack_size != 0) {
      apply += ", " + std::to_string(this->device_stack_size);
    }
    return apply + ")";
  }  // end of getApply

  OpenMPBackend::~OpenMPBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
