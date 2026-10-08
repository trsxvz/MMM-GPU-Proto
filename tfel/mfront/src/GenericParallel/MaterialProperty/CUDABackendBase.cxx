/*!
 * \file   CUDABackendBase.cxx
 * \brief  This file implements the `CUDABackendBase` class
 * \author Thomas Helfer
 * \date   07/04/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <ostream>
#include "MFront/VariableDescription.hxx"
#include "MFront/MaterialPropertyDescription.hxx"
#include "MFront/MaterialPropertyParametersHandler.hxx"
#include "MFront/GenericParallelMaterialPropertyInterface.hxx"
#include "MFront/GenericParallel/MaterialProperty/CUDABackendBase.hxx"

namespace mfront::generic_parallel::material_property {

  CUDABackendBase::CUDABackendBase() = default;

  CUDABackendBase::CUDABackendBase(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addStrictlyPositiveIntegerCheck("number_of_threads_per_block");
    validator.validate(opts);
    if (contains(opts, "number_of_threads_per_block")) {
      this->number_of_threads_per_block =
          get<int>(opts, "number_of_threads_per_block");
    }
  }  // end of CUDABackendBase

  void CUDABackendBase::writeSpecificIncludesInHeaderFile(
      std::ostream&,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {
  }  // end of writeSpecificIncludesInSourceFile

  void CUDABackendBase::writeSpecificIncludesInSourceFile(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {
    os << "#include <memory>\n"
       << "#include <algorithm>\n"
       << "#include<type_traits>\n"
       << "#include\"TFEL/FSAlgorithm/copy.hxx\"\n"
       << "#include\"TFEL/FSAlgorithm/fill.hxx\"\n\n";
  }  // end of writeSpecificIncludesInSourceFile

  void CUDABackendBase::writeGlobalFunctions(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd) const {
    for (const auto b : {true, false}) {
      this->writeGlobalFunctions2(os, i, mpd, fd, b);
    }
  }  // end of writeGlobalFunctions

  void CUDABackendBase::writeGlobalFunctions2(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd,
      const bool treatStrides) const {
    const auto name = i.getFunctionName(mpd);
    const auto types = i.getTypesDescription();
    const auto kernel =
        "mfront_" + name + "_kernel" + (treatStrides ? "" : "_without_strides");
    os << "__global__ void " << kernel << "(" << types.real_type
       << "* const mfront_output";
    if (treatStrides) {
      os << ",\nconst " << types.integer_type << " mfront_output_stride";
    }
    if (isBoundsCheckingRequired(mpd)) {
      os << ",\nint* const mfront_bounds_statuses"
         << ",\nconst " << types.integer_type << " mfront_out_of_bounds_policy";
    }
    if ((!areParametersTreatedAsStaticVariables(mpd)) &&
        (!mpd.parameters.empty())) {
      const auto pcn = getMaterialPropertyParametersClassName(name);
      os << ",\nconst " << i.getInterfaceInternalNamespace() << "::" << pcn
         << " mfront_parameters";
    }
    if ((treatStrides) && (!mpd.inputs.empty())) {
      os << ",\nconst std::array<" << types.integer_type << ", "
         << mpd.inputs.size() << "> mfront_strides";
    }
    for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
      os << ",\nconst " << types.real_type << "* const mfront_arg" << idx;
    }
    os << ",\nconst int mfront_npoints){\n";
    // the body of the kernel is only executed in device code
    writeBeginningOfMaterialPropertyBody(os, mpd, fd, "double", true, true);
    os << "auto mfront_idx = blockIdx.x * blockDim.x + threadIdx.x;\n"
       << "if(mfront_idx >= mfront_npoints){\n"
       << "return;\n"
       << "}\n";
    this->writeKernelBody(os, i, mpd, treatStrides);
    os << "} // end of " << kernel << "\n\n";
  }  // end of writeGlobalFunctions2

  void CUDABackendBase::writeCxxDeclarations(
      std::ostream&,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {}

  void CUDABackendBase::writeCxxImplementations(
      std::ostream&,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&,
      const FileDescription&) const {}

  void CUDABackendBase::writeCImplementations(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd) const {
    for (const auto b : {true, false}) {
      this->writeCImplementations2(os, i, mpd, fd, b);
    }
  }  // end of writeCImplementations

  void CUDABackendBase::writeCImplementations2(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription&,
      const bool treatStrides) const {
    const auto pm = this->getProgrammingModel();
    this->writeCFunctionBeginning(os, i, mpd, treatStrides);
    if (isBoundsCheckingRequired(mpd)) {
      os << "int *mfront_bounds_statuses;\n"
         << "const auto mfront_" << pm << "_alloc_managed__error = " << pm
         << "MallocManaged(&mfront_bounds_statuses,\n"
         << getStatusesArraySize(mpd) << " * sizeof(int));\n"
         << "if (mfront_" << pm << "_alloc_managed__error != " << pm
         << "Success) {\n"
         << "mfront_output_status->status = -6;\n"
         << "mfront_report(" << pm << "GetErrorString(mfront_" << pm
         << "_alloc_managed__error)"
         << ");\n"
         << "return;\n"
         << "}\n"
         << "std::unique_ptr<int, void (*)(int *const) noexcept> "
         << "mfront_bounds_statuses_ptr(mfront_bounds_statuses, "
         << "+[](int *const mfront_ptr) noexcept {"
         << "std::ignore = " << pm << "Free(mfront_ptr);});"
         << "tfel::fsalgo::fill<" << getStatusesArraySize(mpd)
         << ">::exe(mfront_bounds_statuses, 0);\n";
    }
    this->writeKernelCall(os, i, mpd, treatStrides);
    this->writeResultsAnalysis(os, i, mpd);
    os << "} // end of " << i.getFunctionName(mpd) << "\n\n";
  }  // end of writeCImplementations2

  void CUDABackendBase::writeKernelCall(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const bool treatStrides) const {
    const auto pm = this->getProgrammingModel();
    const auto name = i.getFunctionName(mpd);
    const auto block_size = this->number_of_threads_per_block.value_or(64);
    auto write_kernel_call = [&](const std::string_view npoints,
                                 const bool handleStrides) {
      os << "mfront_" << name << "_kernel"
         << (handleStrides ? "" : "_without_strides") << "<<<((" << npoints
         << " + " << block_size - 1 << ") / " << block_size << "), "
         << block_size << ">>>(mfront_output, ";
      if (handleStrides) {
        os << "mfront_output_stride, ";
      }
      if (isBoundsCheckingRequired(mpd)) {
        os << "mfront_bounds_statuses, mfront_out_of_bounds_policy, ";
      }
      if ((!areParametersTreatedAsStaticVariables(mpd)) &&
          (!mpd.parameters.empty())) {
        const auto hn = getMaterialPropertyParametersHandlerClassName(name);
        os << i.getInterfaceInternalNamespace() << "::"  //
           << hn << "::get" << hn << "(), ";
      }
      if ((handleStrides) && (!mpd.inputs.empty())) {
        os << "mfront_strides, ";
      }
      for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
        os << "mfront_args[" << idx << "], ";
      }
      os << "static_cast<int>(" << npoints << "));\n"
         << "const auto mfront_" << pm << "_launch_error = " << pm
         << "GetLastError();\n"
         << "if (mfront_" << pm << "_launch_error != " << pm << "Success) {\n"
         << "mfront_output_status->status = -6;\n"
         << "mfront_report(" << pm << "GetErrorString(mfront_" << pm
         << "_launch_error));\n"
         << "return;\n"
         << "}\n"
         << "const auto mfront_" << pm << "_kernel_error = " << pm
         << "DeviceSynchronize();\n"
         << "if (mfront_" << pm << "_kernel_error != " << pm << "Success) {\n"
         << "mfront_output_status->status = -6;\n"
         << "mfront_report(" << pm << "GetErrorString(mfront_" << pm
         << "_kernel_error));\n"
         << "return;\n"
         << "}\n";
    };
    this->writeKernelSelection(os, i, mpd, treatStrides, write_kernel_call);
  }  // end of writeKernelCall

  CUDABackendBase::~CUDABackendBase() = default;

}  // end of namespace mfront::generic_parallel::material_property
