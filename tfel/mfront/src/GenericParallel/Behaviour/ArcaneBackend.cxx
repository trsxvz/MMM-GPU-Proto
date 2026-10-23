/*!
 * \file   mfront/src/GenericParallel/Behaviour/ArcaneBackend.cxx
 * \brief  This file implements the arcane backend of the `generic-parallel`
 * behaviour interface
 * \date   23/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <ostream>
#include "MFront/GenericParallel/Behaviour/ArcaneBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  ArcaneBackend::ArcaneBackend(const tfel::utilities::DataMap& opts) {
    // no option
    auto validator = tfel::utilities::DataMapValidator{};
    validator.validate(opts);
  }  // end of ArcaneBackend

  std::string ArcaneBackend::getName() const { return "arcane"; }

  std::string ArcaneBackend::getSourceFileExtension() const {
    // the source is compiled by the CUDA compiler, as with the `cuda`
    // backend, and the sources of the `generic` interface by the C++
    // compiler
    return "cu";
  }  // end of getSourceFileExtension

  void ArcaneBackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include \"MFront/GenericParallelBehaviour/Arcane.hxx\"\n";
  }  // end of writeSpecificIncludes

  std::string ArcaneBackend::getApply(const std::string_view status) const {
    return "mfront::gpb::arcane::makeApply(" + std::string{status} + ")";
  }  // end of getApply

  std::string ArcaneBackend::getMemorySpace() const {
    // the device of the queue is only known at runtime
    return "mfront::generic_parallel::device_memory_space";
  }  // end of getMemorySpace

  std::string ArcaneBackend::getExecutionPolicyType() const {
    return "::mfront::generic_parallel::ArcaneBackendExecutionPolicy";
  }  // end of getExecutionPolicyType

  std::string ArcaneBackend::getSharedDataCopy(
      const std::string_view getter, const std::string_view execution) const {
    return "mfront::gpb::arcane::getSharedDataCopy<" + std::string{getter} +
           ">(" + std::string{execution} + ")";
  }  // end of getSharedDataCopy

  ArcaneBackend::~ArcaneBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
