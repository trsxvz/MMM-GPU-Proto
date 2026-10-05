/*!
 * \file   mfront/src/GenericParallel/MaterialProperty/OpenMPBackend.cxx
 * \brief  This file implements the backend of the `generic-parallel`
 * material property interface based on OpenMP
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
#include "MFront/GenericParallel/MaterialProperty/OpenMPBackend.hxx"

namespace mfront::generic_parallel::material_property {

  OpenMPBackend::OpenMPBackend() = default;

  OpenMPBackend::OpenMPBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addDataTypeValidator<bool>("offload");
    validator.addStrictlyPositiveIntegerCheck("device_stack_size");
    validator.validate(opts);
    if (contains(opts, "offload")) {
      this->offload = get<bool>(opts, "offload");
    }
    if (contains(opts, "device_stack_size")) {
      this->device_stack_size = get<int>(opts, "device_stack_size");
    }
  }  // end of OpenMPBackend

  std::string OpenMPBackend::getName() const { return "openmp"; }

  void OpenMPBackend::writeSpecificIncludesInSourceFile(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {
    os << "#include<vector>\n"
       << "#include<algorithm>\n"
       << "#include\"TFEL/FSAlgorithm/copy.hxx\"\n";
    if ((this->offload) && (this->device_stack_size != 0)) {
      os << "#include\"MFront/GenericParallel/DeviceStackSize.hxx\"\n";
    }
    os << "\n";
  }  // end of writeSpecificIncludesInSourceFile

  void OpenMPBackend::writeKernelCall(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd,
      const bool treatStrides) const {
    const auto types = i.getTypesDescription();
    if ((this->offload) && (this->device_stack_size != 0)) {
      os << "if (!mfront::generic_parallel::setMinimalCUDAStackSize("
         << this->device_stack_size << ")) {\n"
         << "mfront_output_status->status = -6;\n"
         << "mfront_report(\"unable to set the stack size of the threads\");\n"
         << "return;\n"
         << "}\n";
    }
    if (this->isBoundsCheckingRequired(mpd)) {
      os << "auto mfront_bounds_statuses_storage = std::vector<int>{};\n"
         << "mfront_bounds_statuses_storage.resize("
         << getStatusesArraySize(mpd) << ", 0);\n";
    }
    this->writeKernels(os, i, mpd, fd, treatStrides);
    auto write_loop = [this, &os, &types](const std::string_view npoints,
                                          const bool handleStrides) {
      const auto kernel = std::string{
          handleStrides ? "mfront_kernel" : "mfront_kernel_without_strides"};
      os << "{\n"
         << "// loop over the points\n"
         << "const auto mfront_n = static_cast<" << types.integer_type << ">("
         << npoints << ");\n";
      if (this->offload) {
        os << "#pragma omp target teams distribute parallel for "
           << "map(to : " << kernel << ")\n";
      } else {
        os << "#pragma omp parallel for\n";
      }
      os << "for (" << types.integer_type
         << " mfront_i = 0; mfront_i < mfront_n; ++mfront_i) {\n"
         << kernel << "(mfront_i);\n"
         << "}\n"
         << "}\n";
    };
    this->writeKernelSelection(os, i, mpd, treatStrides, write_loop);
  }  // end of writeKernelCall

  OpenMPBackend::~OpenMPBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::material_property
