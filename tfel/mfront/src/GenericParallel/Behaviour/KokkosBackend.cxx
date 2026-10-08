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
#include "MFront/GenericParallel/BackendOptions.hxx"
#include "MFront/GenericParallel/Behaviour/KokkosBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  KokkosBackend::KokkosBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    declareKokkosExecutionSpaceOption(validator);
    validator.validate(opts);
    this->execution_space = getKokkosExecutionSpace(opts);
  }  // end of KokkosBackend

  std::string KokkosBackend::getName() const { return "kokkos"; }

  std::string KokkosBackend::getSourceFileExtension() const { return "cxx"; }

  void KokkosBackend::writeSpecificIncludes(std::ostream& os) const {
    os << "#include \"MFront/GenericParallelBehaviour/Kokkos.hxx\"\n";
  }  // end of writeSpecificIncludes

  std::string KokkosBackend::getApply(const std::string_view status) const {
    return "mfront::gpb::kokkos::makeApply<" + this->execution_space + ">(" +
           std::string{status} + ")";
  }  // end of getApply

  std::string KokkosBackend::getMemorySpace() const {
    // the execution space may access the memory of the host or not
    return "(Kokkos::SpaceAccessibility<" + this->execution_space +
           ", Kokkos::HostSpace>::accessible ? "
           "mfront::generic_parallel::host_memory_space : "
           "mfront::generic_parallel::device_memory_space)";
  }  // end of getMemorySpace

  std::string KokkosBackend::getSharedDataCopy(const std::string_view getter,
                                               const std::string_view) const {
    return "mfront::gpb::kokkos::getSharedDataCopy<" + this->execution_space +
           ", " + std::string{getter} + ">()";
  }  // end of getSharedDataCopy

  KokkosBackend::~KokkosBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::behaviour
