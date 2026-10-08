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
#include "MFront/GenericParallel/BackendOptions.hxx"
#include "MFront/GenericParallel/MaterialProperty/KokkosBackend.hxx"

namespace mfront::generic_parallel::material_property {

  KokkosBackend::KokkosBackend() = default;

  KokkosBackend::KokkosBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    declareKokkosExecutionSpaceOption(validator);
    validator.validate(opts);
    this->execution_space = getKokkosExecutionSpace(opts);
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
       << "using mfront_execution_space = " << this->execution_space
       << ";\n"
       // instance of the execution space given by the execution options
       << "const auto mfront_space = (mfront_options.queue == nullptr) ? "
       << "mfront_execution_space{} : "
       << "*static_cast<const mfront_execution_space*>("
       << "mfront_options.queue);\n";
    if (this->isBoundsCheckingRequired(mpd)) {
      // the statuses are written by the kernels and read on the host. Kokkos
      // uses one device per process. The blocks are freed when Kokkos is
      // finalized.
      writeBoundsStatusesDeclaration(
          os, mpd, "int", "0",
          "auto* const mfront_ptr = static_cast<int*>("
          "Kokkos::kokkos_malloc<Kokkos::SharedSpace>("
          "\"mfront_bounds_statuses\", mfront_n * sizeof(int)));\n"
          "Kokkos::push_finalize_hook([mfront_ptr] {\n"
          "Kokkos::kokkos_free<Kokkos::SharedSpace>(mfront_ptr);\n"
          "});\n"
          "return mfront_ptr;");
    }
    this->writeKernels(os, i, mpd, fd, treatStrides);
    auto write_loop = [&os, &types](const std::string_view npoints,
                                    const bool handleStrides) {
      os << "// loop over the points\n"
         << "Kokkos::parallel_for(\"mfront_kernel\", "
         << "Kokkos::RangePolicy<mfront_execution_space>(mfront_space, 0, "
         << "static_cast<" << types.integer_type << ">(" << npoints << ")), "
         << "mfront_kernel";
      if (!handleStrides) {
        os << "_without_strides";
      }
      // an asynchronous evaluation returns once the kernel is submitted
      os << ");\n"
         << "if (!mfront_asynchronous) {\n"
         << "mfront_space.fence();\n"
         << "}\n";
    };
    this->writeKernelSelection(os, i, mpd, treatStrides, write_loop);
  }  // end of writeKernelCall

  std::string KokkosBackend::getMemorySpace() const {
    // the execution space may access the memory of the host or not
    return "(Kokkos::SpaceAccessibility<" + this->execution_space +
           ", Kokkos::HostSpace>::accessible ? "
           "mfront::generic_parallel::host_memory_space : "
           "mfront::generic_parallel::device_memory_space)";
  }  // end of getMemorySpace

  KokkosBackend::~KokkosBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::material_property
