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
#include "TFEL/Raise.hxx"
#include "MFront/GenericParallel/Behaviour/ParallelSTLBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  ParallelSTLBackend::ParallelSTLBackend() = default;

  ParallelSTLBackend::ParallelSTLBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addDataTypeValidator<std::string>("execution_policy");
    validator.addStrictlyPositiveIntegerCheck("device_stack_size");
    validator.validate(opts);
    if (contains(opts, "device_stack_size")) {
      this->device_stack_size = get<int>(opts, "device_stack_size");
    }
    if (!contains(opts, "execution_policy")) {
      return;
    }
    const auto p = get<std::string>(opts, "execution_policy");
    if ((p == "seq") || (p == "sequenced_policy")) {
      this->execution_policy = "std::execution::seq";
    } else if ((p == "unseq") || (p == "unsequenced_policy")) {
      this->execution_policy = "std::execution::unseq";
    } else if ((p == "par") || (p == "parallel_policy")) {
      this->execution_policy = "std::execution::par";
    } else if ((p == "par_unseq") || (p == "parallel_unsequenced_policy")) {
      this->execution_policy = "std::execution::par_unseq";
    } else {
      tfel::raise(
          "ParallelSTLBackend::ParallelSTLBackend: unsupported execution "
          "policy '" +
          p +
          "'. Valid execution policies are: 'sequenced_policy' (or "
          "'seq'), 'unsequenced_policy' (or 'unseq'), 'parallel_policy' (or "
          "'par') and 'parallel_unsequenced_policy' (or 'par_unseq')");
    }
  }  // end of ParallelSTLBackend

  std::string ParallelSTLBackend::getName() const { return "stlpar"; }

  std::string ParallelSTLBackend::getSourceFileExtension() const {
    return "cxx";
  }  // end of getSourceFileExtension

  void ParallelSTLBackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include <execution>\n"
       << "#include \"MFront/GenericParallelBehaviour/ParallelSTL.hxx\"\n";
  }  // end of writeSpecificIncludes

  void ParallelSTLBackend::writeCall(std::ostream& os,
                                     const std::string_view f,
                                     const std::string_view template_arguments,
                                     const std::string_view arguments) const {
    os << "return mfront::gpb::stlpar::" << f << "<" << template_arguments
       << ">(" << arguments << ", " << this->execution_policy;
    if (this->device_stack_size != 0) {
      os << ", " << this->device_stack_size;
    }
    os << ");\n";
  }  // end of writeCall

  ParallelSTLBackend::~ParallelSTLBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
