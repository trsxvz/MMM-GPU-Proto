/*!
 * \file   mfront/src/GenericParallel/MaterialProperty/KokkosBackend.cxx
 * \brief  This file implements the backend of the `generic-parallel`
 * material property interface based on Kokkos
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <ostream>
#include "MFront/MaterialPropertyDescription.hxx"
#include "MFront/GenericParallelMaterialPropertyInterface.hxx"
#include "MFront/GenericParallel/MaterialProperty/KokkosBackend.hxx"

namespace mfront::generic_parallel::material_property {

  KokkosBackend::KokkosBackend() = default;

  KokkosBackend::KokkosBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addDataTypeValidator<std::string>("execution_space");
    validator.validate(opts);
    if (contains(opts, "execution_space")) {
      this->execution_space = get<std::string>(opts, "execution_space");
    }
  }  // end of KokkosBackend

  std::string KokkosBackend::getName() const { return "kokkos"; }

  void KokkosBackend::writeSpecificIncludesInSourceFile(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {
    os << "#include<algorithm>\n"
       << "#include<Kokkos_Core.hpp>\n"
       << "#include\"TFEL/FSAlgorithm/copy.hxx\"\n\n";
  }  // end of writeSpecificIncludesInSourceFile

  std::string KokkosBackend::getKernelSpecifiers() const {
    // the kernels are called in the execution space
    return " KOKKOS_FUNCTION";
  }  // end of getKernelSpecifiers

  void KokkosBackend::writeKernelCall(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd,
      const bool treatStrides) const {
    const auto types = i.getTypesDescription();
    os << "if (!Kokkos::is_initialized()) {\n"
       << "mfront_output_status->status = -6;\n"
       << "mfront_report(\"Kokkos is not initialized\");\n"
       << "return;\n"
       << "}\n"
       << "using mfront_execution_space = " << this->execution_space << ";\n";
    if (this->isBoundsCheckingRequired(mpd)) {
      // the statuses are written by the kernels and read on the host
      os << "auto mfront_bounds_statuses_storage = "
         << "Kokkos::View<int*, Kokkos::SharedSpace>("
         << "\"mfront_bounds_statuses\", " << getStatusesArraySize(mpd)
         << ");\n";
    }
    this->writeKernels(os, i, mpd, fd, treatStrides);
    auto write_loop = [&os, &types](const std::string_view npoints,
                                    const bool handleStrides) {
      os << "// loop over the points\n"
         << "Kokkos::parallel_for(\"mfront_kernel\", "
         << "Kokkos::RangePolicy<mfront_execution_space>(0, "
         << "static_cast<" << types.integer_type << ">(" << npoints << ")), "
         << "mfront_kernel";
      if (!handleStrides) {
        os << "_without_strides";
      }
      os << ");\n"
         << "Kokkos::fence();\n";
    };
    this->writeKernelSelection(os, i, mpd, treatStrides, write_loop);
  }  // end of writeKernelCall

  KokkosBackend::~KokkosBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::material_property
