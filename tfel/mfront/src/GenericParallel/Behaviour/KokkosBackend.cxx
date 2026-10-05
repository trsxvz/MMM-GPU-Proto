/*!
 * \file   mfront/src/GenericParallel/Behaviour/KokkosBackend.cxx
 * \brief  This file implements the Kokkos backend of the `generic-parallel`
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
#include "MFront/GenericParallel/Behaviour/KokkosBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  KokkosBackend::KokkosBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addDataTypeValidator<std::string>("execution_space");
    validator.validate(opts);
    if (contains(opts, "execution_space")) {
      this->execution_space = get<std::string>(opts, "execution_space");
    }
  }  // end of KokkosBackend

  std::string KokkosBackend::getName() const { return "kokkos"; }

  std::string KokkosBackend::getSourceFileExtension() const { return "cxx"; }

  void KokkosBackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include \"MFront/GenericParallelBehaviour/Kokkos.hxx\"\n";
  }  // end of writeSpecificIncludes

  void KokkosBackend::writeCall(std::ostream& os,
                                const std::string_view f,
                                const std::string_view template_arguments,
                                const std::string_view arguments) const {
    os << "return mfront::gpb::kokkos::" << f << "<" << template_arguments
       << ", " << this->execution_space << ">(" << arguments << ");\n";
  }  // end of writeCall

  KokkosBackend::~KokkosBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
