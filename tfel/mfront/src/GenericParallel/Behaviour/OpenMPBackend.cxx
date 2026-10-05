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
#include "MFront/GenericParallel/Behaviour/OpenMPBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  OpenMPBackend::OpenMPBackend() = default;

  OpenMPBackend::OpenMPBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addDataTypeValidator<bool>("offload");
    validator.addStrictlyPositiveIntegerCheck("device_stack_size");
    validator.validate(opts);
    if (contains(opts, "offload")) {
      this->offload = get<bool>(opts, "offload");
    }
    if (contains(opts, "device_stack_size")) {
      this->device_stack_size = get<int>(opts, "device_stack_size");
    }
  }  // end of OpenMPBackend

  std::string OpenMPBackend::getName() const { return "openmp"; }

  std::string OpenMPBackend::getSourceFileExtension() const {
    return "cxx";
  }  // end of getSourceFileExtension

  void OpenMPBackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include \"MFront/GenericParallelBehaviour/OpenMP.hxx\"\n";
  }  // end of writeSpecificIncludes

  void OpenMPBackend::writeCall(std::ostream& os,
                                const std::string_view f,
                                const std::string_view template_arguments,
                                const std::string_view arguments) const {
    // the offload option is a template argument, so that the code
    // executed on the host does not depend on the support of offloading by
    // the OpenMP runtime
    os << "return mfront::gpb::openmp::" << f << "<" << template_arguments
       << ", " << (this->offload ? "true" : "false") << ">(" << arguments;
    if (this->device_stack_size != 0) {
      os << ", " << this->device_stack_size;
    }
    os << ");\n";
  }  // end of writeCall

  OpenMPBackend::~OpenMPBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
