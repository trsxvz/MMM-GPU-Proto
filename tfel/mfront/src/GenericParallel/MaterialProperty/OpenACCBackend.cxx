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
      const FileDescription& fd) const {
    const auto types = i.getTypesDescription();
    // only the common options of the execution policy are used: the
    // computations are performed on the current device of OpenACC
    if (this->device_stack_size != 0) {
      os << "if (!mfront::generic_parallel::openacc::setMinimalStackSize("
         << this->device_stack_size << ", mfront_preserve_stack_size)) {\n";
      writeBatchError(os, -7,
                      "mfront::generic_parallel::stack_size_error_message");
      os << "}\n";
    }
    const auto checks = this->isBoundsCheckingRequired(mpd);
    if (checks) {
      // the statuses of the points are merged by a reduction, as for the
      // behaviours
      writeLocalStatisticsDeclaration(os);
    }
    this->writeKernels(os, i, mpd, fd);
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
        os << "using mfront_key_type = std::conditional_t<"
           << "sizeof(unsigned long) == sizeof(unsigned long long), "
           << "unsigned long, unsigned long long>;\n"
           << "auto mfront_number_of_failures = 0ull;\n"
           << "auto mfront_first_failure = "
           << "static_cast<mfront_key_type>(::mfront::gpmp::no_event);\n"
           << "auto mfront_number_of_unreliable_results = 0ull;\n"
           << "auto mfront_first_unreliable_result = "
           << "static_cast<mfront_key_type>(::mfront::gpmp::no_event);\n";
      }
      const auto reductions =
          std::string{checks ? " reduction(+ : mfront_number_of_failures, "
                               "mfront_number_of_unreliable_results) "
                               "reduction(min : mfront_first_failure, "
                               "mfront_first_unreliable_result)"
                             : ""};
      os << "#pragma acc parallel loop copyin(" << kernel << ")" << reductions
         << "\n";
      os << "for (" << types.integer_type
         << " mfront_i = 0; mfront_i < mfront_n; ++mfront_i) {\n";
      if (checks) {
        os << "const auto mfront_point_status = " << kernel << "(mfront_i);\n"
           << "if (mfront_point_status.hasEvent()) {\n"
           << "const auto mfront_key = static_cast<mfront_key_type>("
           << "::mfront::gpmp::makeKey("
           << "static_cast<unsigned long long>(mfront_i), "
           << "mfront_point_status.check));\n"
           << "if (mfront_point_status.isFailure()) {\n"
           << "mfront_number_of_failures += 1;\n"
           << "mfront_first_failure = (mfront_key < mfront_first_failure) ? "
           << "mfront_key : mfront_first_failure;\n"
           << "} else {\n"
           << "mfront_number_of_unreliable_results += 1;\n"
           << "mfront_first_unreliable_result = "
           << "(mfront_key < mfront_first_unreliable_result) ? "
           << "mfront_key : mfront_first_unreliable_result;\n"
           << "}\n"
           << "}\n";
      } else {
        os << kernel << "(mfront_i);\n";
      }
      os << "}\n";
      if (checks) {
        os << "*mfront_statistics = ::mfront::gpmp::merge(*mfront_statistics, "
           << "mfront_gpmp_Statistics{"
           << ".number_of_failures = mfront_number_of_failures, "
           << ".first_failure = "
           << "static_cast<unsigned long long>(mfront_first_failure), "
           << ".number_of_unreliable_results = "
           << "mfront_number_of_unreliable_results, "
           << ".first_unreliable_result = "
           << "static_cast<unsigned long long>("
           << "mfront_first_unreliable_result)});\n";
      }
      os << "}\n";
    };
    this->writeKernelSelection(os, i, mpd, write_loop);
  }  // end of writeKernelCall

  std::string OpenACCBackend::getMemorySpace() const {
    return "mfront::generic_parallel::getOpenACCMemorySpace()";
  }  // end of getMemorySpace

  std::string OpenACCBackend::getExecutionPolicyType() const {
    return "::mfront::generic_parallel::OpenACCBackendExecutionPolicy";
  }  // end of getExecutionPolicyType

  std::string OpenACCBackend::getExecutionPolicyHeader() const {
    return "MFront/GenericParallel/OpenACCBackendExecutionPolicy.hxx";
  }  // end of getExecutionPolicyHeader

  OpenACCBackend::~OpenACCBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::material_property
