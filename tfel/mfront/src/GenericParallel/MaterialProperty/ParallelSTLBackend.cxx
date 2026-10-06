/*!
 * \file   ParallelSTLBackend.cxx
 * \brief  This file implements the `ParallelSTLBackend` class
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
#include "MFront/GenericParallel/BackendOptions.hxx"
#include "MFront/GenericParallel/MaterialProperty/ParallelSTLBackend.hxx"

namespace mfront::generic_parallel::material_property {

  ParallelSTLBackend::ParallelSTLBackend()
      : execution_policy("std::execution::par_unseq") {
  }  // end of ParallelSTLBackend

  ParallelSTLBackend::ParallelSTLBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    declareExecutionPolicyOption(validator);
    declareDeviceStackSizeOption(validator);
    validator.validate(opts);
    this->execution_policy = getExecutionPolicy(opts);
    this->device_stack_size = getDeviceStackSize(opts);
  }  // end of ParallelSTLBackend

  std::string ParallelSTLBackend::getName() const {
    return "stlpar";
  }  // end of getName

  void ParallelSTLBackend::writeSpecificIncludesInSourceFile(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {
    os << "#include<ranges>\n"
       << "#include<concepts>\n"
       << "#include<execution>\n"
       << "#include<type_traits>\n"
       << "#include\"TFEL/FSAlgorithm/copy.hxx\"\n"
       << "#include\"MFront/GenericParallel/CountingIterator.hxx\"\n";
    if (this->device_stack_size != 0) {
      os << "#include\"MFront/GenericParallel/DeviceStackSize.hxx\"\n";
    }
    os << "\n";
  }  // end of writeSpecificIncludesInSourceFile

  void ParallelSTLBackend::writeKernelCall(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd,
      const bool treatStrides) const {
    const auto types = i.getTypesDescription();
    const auto requiresBoundsCheck = this->isBoundsCheckingRequired(mpd);
    if (this->device_stack_size != 0) {
      os << "if (!mfront::generic_parallel::setMinimalCUDAStackSize("
         << this->device_stack_size << ")) {\n"
         << "mfront_output_status->status = -6;\n"
         << "mfront_report(\"unable to set the stack size of the threads\");\n"
         << "return;\n"
         << "}\n";
    }
    if (requiresBoundsCheck) {
      // we do need to use heap allocations as nvc++ would automatically use
      // managed memory
      os << "auto mfront_bounds_statuses_storage = std::vector<int>{};\n"
         << "mfront_bounds_statuses_storage.resize("
         << getStatusesArraySize(mpd) << ", 0);\n";
    }
    this->writeKernels(os, i, mpd, fd, treatStrides);
    auto write_kernel_call = [&os, this, &types](const std::string_view npoints,
                                                 const bool handleStrides) {
      // the iterators of std::views::iota can't be used: with libstdc++,
      // the parallel algorithms run sequentially with them (see the
      // documentation of mfront::generic_parallel::CountingIterator)
      os << "// loop over the points\n"
         << "#ifdef __ADAPTIVECPP__\n"
         << "using mfront_index_type = int;\n"
         << "#else /* __ADAPTIVECPP__ */ \n"
         << "using mfront_index_type = " << types.integer_type << ";\n"
         << "#endif /* __ADAPTIVECPP__ */ \n"
         << "using mfront_index_iterator = "
         << "mfront::generic_parallel::CountingIterator<mfront_index_type>;\n"
         << "std::for_each(" << this->execution_policy
         << ", mfront_index_iterator{0}, "
         << "mfront_index_iterator{"
         << "static_cast<mfront_index_type>(" << npoints << ")}, "
         << "mfront_kernel";
      if (!handleStrides) {
        os << "_without_strides";
      }
      os << ");\n";
    };
    this->writeKernelSelection(os, i, mpd, treatStrides, write_kernel_call);
  }  // end of writeKernelCall

  ParallelSTLBackend::~ParallelSTLBackend() = default;

}  // end of namespace mfront::generic_parallel::material_property
