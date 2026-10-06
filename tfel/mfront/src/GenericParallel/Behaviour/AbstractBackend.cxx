/*!
 * \file   mfront/src/GenericParallel/Behaviour/AbstractBackend.cxx
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
#include "MFront/CodeGeneratorUtilities.hxx"
#include "MFront/GenericParallel/Behaviour/AbstractBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  void AbstractBackend::writeSymbols(std::ostream& os,
                                     const std::string_view f) const {
    const auto n = std::string{f} + "_generic_parallel_";
    exportStringSymbol(os, n + "backend", this->getName());
    // the version is defined in `MFront/GenericParallelBehaviour/BatchData.h`
    exportSymbol(os, "unsigned short", n + "api_version",
                 "MFRONT_GPB_API_VERSION");
  }  // end of writeSymbols

  std::string AbstractBackend::getSharedDataCopy(const std::string_view value,
                                                 const std::string_view) const {
    return "mfront::gpb::copyToHeap(" + std::string{value} + ")";
  }  // end of getSharedDataCopy

  void AbstractBackend::writeCall(std::ostream& os,
                                  const std::string_view f,
                                  const std::string_view template_arguments,
                                  const std::string_view status,
                                  const std::string_view arguments) const {
    os << "return mfront::gpb::" << f << "<" << template_arguments << ">("
       << status << ", " << arguments << ", " << this->getApply(status)
       << ");\n";
  }  // end of writeCall

  AbstractBackend::~AbstractBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
