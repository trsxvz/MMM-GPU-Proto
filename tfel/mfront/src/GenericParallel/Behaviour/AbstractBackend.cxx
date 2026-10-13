/*!
 * \file   mfront/src/GenericParallel/Behaviour/AbstractBackend.cxx
 * \brief  This file implements the `AbstractBackend` class of the
 * `generic-parallel` behaviour interface
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

  std::string AbstractBackend::getLibraryName(
      const std::string_view library,
      const std::string_view material,
      const std::string_view default_name) const {
    const auto base = !library.empty()    ? library
                      : !material.empty() ? material
                                          : default_name;
    return "GenericParallel" + std::string{base} + '-' + this->getName();
  }  // end of getLibraryName

  std::string AbstractBackend::getHeaderFileName(
      const std::string_view name) const {
    return std::string{name} + "-generic-parallel-" + this->getName() + ".hxx";
  }  // end of getHeaderFileName

  std::string AbstractBackend::getSourceFileName(
      const std::string_view name) const {
    return std::string{name} + "-generic-parallel-" + this->getName() + '.' +
           this->getSourceFileExtension();
  }  // end of getSourceFileName

  void AbstractBackend::writeSymbols(std::ostream& os,
                                     const std::string_view f) const {
    const auto n = std::string{f} + "_generic_parallel_";
    exportStringSymbol(os, n + "backend", this->getName());
    exportSymbol(os, "const char*", n + "memory_space", this->getMemorySpace());
    // the version is defined in `MFront/GenericParallelBehaviour/BatchData.h`
    exportSymbol(os, "unsigned short", n + "api_version",
                 "MFRONT_GPB_API_VERSION");
  }  // end of writeSymbols

  std::string AbstractBackend::getSharedDataCopy(const std::string_view getter,
                                                 const std::string_view) const {
    return "mfront::gpb::getSharedDataCopy<" + std::string{getter} + ">()";
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
