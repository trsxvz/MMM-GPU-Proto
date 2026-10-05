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
#include "TFEL/Raise.hxx"
#include "MFront/GenericParallel/Behaviour/SYCLBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  SYCLBackend::SYCLBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addDataTypeValidator<std::string>("data_transfer_policy");
    validator.addStrictlyPositiveIntegerCheck("device_stack_size");
    validator.validate(opts);
    if (contains(opts, "device_stack_size")) {
      this->device_stack_size = get<int>(opts, "device_stack_size");
    }
    if (contains(opts, "data_transfer_policy")) {
      const auto p = get<std::string>(opts, "data_transfer_policy");
      if ((p == "auto") || (p == "automatic")) {
        this->automatic_data_transfers = true;
      } else if (p == "none") {
        this->automatic_data_transfers = false;
      } else {
        tfel::raise("SYCLBackend::SYCLBackend: invalid data transfer policy '" +
                    p +
                    "'. Valid policies are 'automatic' (or 'auto') and "
                    "'none'");
      }
    }
  }  // end of SYCLBackend

  std::string SYCLBackend::getName() const { return "sycl"; }

  std::string SYCLBackend::getSourceFileExtension() const { return "cxx"; }

  void SYCLBackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include \"MFront/GenericParallelBehaviour/SYCL.hxx\"\n";
  }  // end of writeSpecificIncludes

  void SYCLBackend::writeCall(std::ostream& os,
                              const std::string_view f,
                              const std::string_view template_arguments,
                              const std::string_view arguments) const {
    os << "return mfront::gpb::sycl::" << f << "<" << template_arguments << ">("
       << arguments << ", mfront::gpb::sycl::DataTransferPolicy::"
       << (this->automatic_data_transfers ? "AUTOMATIC" : "NONE");
    if (this->device_stack_size != 0) {
      os << ", " << this->device_stack_size;
    }
    os << ");\n";
  }  // end of writeCall

  SYCLBackend::~SYCLBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
