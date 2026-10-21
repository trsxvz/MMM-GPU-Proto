/*!
 * \file   mfront/src/GenericParallel/MaterialProperty/RAJABackend.cxx
 * \brief  This file implements the RAJA backend of the `generic-parallel`
 * material property interface
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
#include "MFront/GenericParallel/MaterialProperty/RAJABackend.hxx"

namespace mfront::generic_parallel::material_property {

  RAJABackend::RAJABackend() = default;

  RAJABackend::RAJABackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    declareRAJAExecutionPolicyOption(validator);
    validator.validate(opts);
    this->execution_policy = getRAJAExecutionPolicy(opts);
  }  // end of RAJABackend

  std::string RAJABackend::getName() const { return "raja"; }

  void RAJABackend::writeSpecificIncludesInSourceFile(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {
    os << "#include<algorithm>\n"
       << "#include\"RAJA/RAJA.hpp\"\n"
       << "#include\"TFEL/FSAlgorithm/copy.hxx\"\n"
       << "#include\"MFront/GenericParallel/RAJAExecutionPolicy.hxx\"\n\n";
  }  // end of writeSpecificIncludesInSourceFile

  std::string RAJABackend::getKernelSpecifiers() const {
    // the kernels are called by the loops of RAJA, on the host or on the
    // device
    return " RAJA_HOST_DEVICE";
  }  // end of getKernelSpecifiers

  void RAJABackend::writeKernelCall(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd) const {
    const auto types = i.getTypesDescription();
    os << "using mfront_execution_policy = " << this->execution_policy
       << ";\n"
       // resource given by the execution options
       << "auto mfront_resource = "
       << "::mfront::generic_parallel::raja::getResource<"
       << "mfront_execution_policy>(mfront_options.queue);\n";
    if (this->isBoundsCheckingRequired(mpd)) {
      // the statistics are written by the kernels and read on the host. The
      // blocks of managed memory, by device, are kept until the end of the
      // process.
      writeStatisticsDeclaration(
          os, "int",
          "::mfront::generic_parallel::raja::getDevice(mfront_resource)",
          "return mfront_resource.template allocate<mfront_gpmp_Statistics>("
          "1, camp::resources::MemoryAccess::Managed);");
    }
    this->writeKernels(os, i, mpd, fd);
    auto write_loop = [&os, &types](const std::string_view npoints,
                                    const bool handleStrides) {
      const auto kernel = std::string{
          handleStrides ? "mfront_kernel" : "mfront_kernel_without_strides"};
      // stack size of the threads, as for the behaviours
      os << "const auto mfront_stack_size_error = "
         << "::mfront::generic_parallel::raja::ensureStackSize<"
         << "mfront_execution_policy>(" << kernel
         << ", mfront_options.preserve_stack_size != 0);\n"
         << "if (!mfront_stack_size_error.empty()) {\n";
      writeBatchError(os, -7, "mfront_stack_size_error");
      os << "}\n";
      os << "// loop over the points\n"
         << "RAJA::forall<mfront_execution_policy>(mfront_resource, "
         << "RAJA::TypedRangeSegment<" << types.integer_type << ">("
         << types.integer_type << "{0}, static_cast<" << types.integer_type
         << ">(" << npoints << ")), " << kernel;
      // an asynchronous evaluation returns once the kernel is submitted
      os << ");\n"
         << "if (!mfront_asynchronous) {\n"
         << "mfront_resource.wait();\n"
         << "}\n";
    };
    this->writeKernelSelection(os, i, mpd, write_loop);
  }  // end of writeKernelCall

  bool RAJABackend::mergesPointStatusesInKernels() const {
    return true;
  }  // end of mergesPointStatusesInKernels

  std::string RAJABackend::getAtomicIncrement(const std::string_view c) const {
    return "RAJA::atomicAdd<RAJA::auto_atomic>(&(" + std::string{c} +
           "), 1ull)";
  }  // end of getAtomicIncrement

  std::string RAJABackend::getAtomicMinimum(const std::string_view c,
                                            const std::string_view v) const {
    return "RAJA::atomicMin<RAJA::auto_atomic>(&(" + std::string{c} + "), " +
           std::string{v} + ")";
  }  // end of getAtomicMinimum

  std::string RAJABackend::getMemorySpace() const {
    // the policies of CUDA and HIP may not access the memory of the host
    return "(mfront::generic_parallel::raja::isDeviceExecutionPolicy<" +
           this->execution_policy +
           "> ? mfront::generic_parallel::device_memory_space : "
           "mfront::generic_parallel::host_memory_space)";
  }  // end of getMemorySpace

  RAJABackend::~RAJABackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::material_property
