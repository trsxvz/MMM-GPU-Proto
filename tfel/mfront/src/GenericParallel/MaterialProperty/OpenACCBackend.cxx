/*!
 * \file   mfront/src/GenericParallel/MaterialProperty/OpenACCBackend.cxx
 * \brief  This file implements the backend of the `generic-parallel`
 * material property interface based on OpenACC
 * \date   20/10/2026
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
#include "MFront/GenericParallel/MaterialProperty/OpenACCBackend.hxx"

namespace mfront::generic_parallel::material_property {

  OpenACCBackend::OpenACCBackend() = default;

  OpenACCBackend::OpenACCBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    declareDeviceStackSizeOption(validator);
    validator.validate(opts);
    this->device_stack_size = getDeviceStackSize(opts);
  }  // end of OpenACCBackend

  std::string OpenACCBackend::getName() const { return "openacc"; }

  void OpenACCBackend::writeSpecificIncludesInSourceFile(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {
    os << "#include<vector>\n"
       << "#include<algorithm>\n"
       << "#include<type_traits>\n"
       << "#include\"TFEL/FSAlgorithm/copy.hxx\"\n";
    if (this->device_stack_size != 0) {
      os << "#include\"MFront/GenericParallel/DeviceStackSize.hxx\"\n"
         << "#include\"MFront/GenericParallel/OpenACCStackSize.hxx\"\n";
    }
    os << "\n";
  }  // end of writeSpecificIncludesInSourceFile

  void OpenACCBackend::writeKernelCall(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd,
      const bool treatStrides) const {
    const auto types = i.getTypesDescription();
    // the queue of the execution options is not used: the computations are
    // performed on the current device of OpenACC
    if (this->device_stack_size != 0) {
      os << "if (!mfront::generic_parallel::openacc::setMinimalStackSize("
         << this->device_stack_size
         << ", mfront_options.preserve_stack_size != 0)) {\n"
         << "mfront_output_status->status = -6;\n"
         << "mfront_report("
         << "mfront::generic_parallel::stack_size_error_message);\n"
         << "return;\n"
         << "}\n";
    }
    const auto checks = this->isBoundsCheckingRequired(mpd);
    if (checks) {
      // the statuses of the points are merged by a reduction, as for the
      // behaviours
      writeLocalStatisticsDeclaration(os);
    }
    this->writeKernels(os, i, mpd, fd, treatStrides);
    auto write_loop = [this, &os, &types, checks](
                          const std::string_view npoints,
                          const bool handleStrides) {
      const auto kernel = std::string{
          handleStrides ? "mfront_kernel" : "mfront_kernel_without_strides"};
      os << "{\n"
         << "// loop over the points\n"
         << "const auto mfront_n = static_cast<" << types.integer_type << ">("
         << npoints << ");\n";
      if (checks) {
        // the minimum reductions are performed on `unsigned long` if this
        // type has the size of `unsigned long long`: the minimum reductions
        // of OpenACC on `unsigned long long` are wrong on NVIDIA GPUs with
        // `nvc++` 25.11
        os << "using mfront_index_type = std::conditional_t<"
           << "sizeof(unsigned long) == sizeof(unsigned long long), "
           << "unsigned long, unsigned long long>;\n"
           << "auto mfront_number_of_failures = 0ull;\n"
           << "auto mfront_first_failure_index = "
           << "static_cast<mfront_index_type>(::mfront::gpmp::no_event);\n"
           << "auto mfront_number_of_unreliable_results = 0ull;\n"
           << "auto mfront_first_violation = "
           << "static_cast<mfront_index_type>(::mfront::gpmp::no_event);\n";
      }
      const auto reductions =
          std::string{checks ? " reduction(+ : mfront_number_of_failures, "
                               "mfront_number_of_unreliable_results) "
                               "reduction(min : mfront_first_failure_index, "
                               "mfront_first_violation)"
                             : ""};
      os << "#pragma acc parallel loop copyin(" << kernel << ")" << reductions
         << "\n";
      os << "for (" << types.integer_type
         << " mfront_i = 0; mfront_i < mfront_n; ++mfront_i) {\n";
      if (checks) {
        os << "const auto mfront_point_status = " << kernel << "(mfront_i);\n"
           << "if (mfront_point_status.hasEvent()) {\n"
           << "if (mfront_point_status.failure) {\n"
           << "mfront_number_of_failures += 1;\n"
           << "const auto mfront_ui = static_cast<mfront_index_type>("
           << "mfront_i);\n"
           << "mfront_first_failure_index = "
           << "(mfront_ui < mfront_first_failure_index) ? mfront_ui : "
           << "mfront_first_failure_index;\n"
           << "} else {\n"
           << "mfront_number_of_unreliable_results += 1;\n"
           << "}\n"
           << "const auto mfront_violation = static_cast<mfront_index_type>("
           << "mfront_point_status.violation);\n"
           << "mfront_first_violation = "
           << "(mfront_violation < mfront_first_violation) ? "
           << "mfront_violation : mfront_first_violation;\n"
           << "}\n";
      } else {
        os << kernel << "(mfront_i);\n";
      }
      os << "}\n";
      if (checks) {
        os << "*mfront_statistics = ::mfront::gpmp::merge(*mfront_statistics, "
           << "mfront_gpmp_Statistics{"
           << ".number_of_failures = mfront_number_of_failures, "
           << ".first_failure_index = "
           << "static_cast<unsigned long long>(mfront_first_failure_index), "
           << ".number_of_unreliable_results = "
           << "mfront_number_of_unreliable_results, "
           << ".first_violation = "
           << "static_cast<unsigned long long>(mfront_first_violation)});\n";
      }
      os << "}\n";
    };
    this->writeKernelSelection(os, i, mpd, treatStrides, write_loop);
  }  // end of writeKernelCall

  std::string OpenACCBackend::getMemorySpace() const {
    return "mfront::generic_parallel::getOpenACCMemorySpace()";
  }  // end of getMemorySpace

  OpenACCBackend::~OpenACCBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::material_property
