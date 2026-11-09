/*!
 * \file   mfront/src/GenericParallel/AbstractBackend.cxx
 * \brief  This file implements the `AbstractBackend` class of the
 * `generic-parallel` interfaces
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
#include "MFront/GenericParallel/AbstractBackend.hxx"

namespace mfront::generic_parallel {

  std::string_view AbstractBackend::getNamespace(const Runtime r) noexcept {
    return r == Runtime::BEHAVIOURS ? "mfront::gpb" : "mfront::gpmp";
  }  // end of getNamespace

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
                                     const std::string_view f,
                                     const Runtime r) const {
    const auto n = std::string{f} + "_generic_parallel_";
    exportStringSymbol(os, n + "backend", this->getName());
    exportSymbol(os, "const char*", n + "memory_space", this->getMemorySpace());
    // identifier of the execution policy expected by the generated functions.
    // As for the memory space, the expression is enclosed in parentheses:
    // the type of the policy may contain commas (for instance, a policy of
    // RAJA with several template arguments), which would split the arguments
    // of the macro exporting the symbol.
    exportSymbol(os, "const char*", n + "execution_policy",
                 "(" + this->getExecutionPolicyType() + "::identifier())");
    // the version is defined in `MFront/GenericParallelBehaviour/BatchData.h`
    // or in `MFront/GenericParallelMaterialProperty/MaterialProperty.h`, and
    // shifted for a numerical type other than double (see
    // `MFront/GenericBehaviour/Types.h` and
    // `MFront/GenericMaterialProperty/Types.h`)
    exportSymbol(os, "unsigned short", n + "api_version",
                 r == Runtime::BEHAVIOURS
                     ? "MFRONT_GPB_API_VERSION + "
                       "MFRONT_GB_REAL_TYPE_API_VERSION_OFFSET"
                     : "MFRONT_GPMP_API_VERSION + "
                       "MFRONT_GMP_REAL_TYPE_API_VERSION_OFFSET");
  }  // end of writeSymbols

  std::string AbstractBackend::getSharedDataCopy(const std::string_view getter,
                                                 const std::string_view) const {
    return "mfront::gpb::getSharedDataCopy<" + std::string{getter} + ">()";
  }  // end of getSharedDataCopy

  void AbstractBackend::writeCall(std::ostream& os,
                                  const Runtime r,
                                  const std::string_view f,
                                  const std::string_view template_arguments,
                                  const std::string_view status,
                                  const std::string_view arguments) const {
    os << "return " << getNamespace(r) << "::" << f;
    if (!template_arguments.empty()) {
      os << "<" << template_arguments << ">";
    }
    os << "(" << status << ", " << arguments << ", " << this->getApply()
       << ");\n";
  }  // end of writeCall

  AbstractBackend::~AbstractBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel
