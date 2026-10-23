/*!
 * \file   mfront/src/GenericParallel/MaterialProperty/ArcaneBackend.cxx
 * \brief  This file implements the arcane backend of the `generic-parallel`
 * material property interface
 * \date   23/10/2026
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
#include "MFront/GenericParallel/MaterialProperty/ArcaneBackend.hxx"

namespace mfront::generic_parallel::material_property {

  ArcaneBackend::ArcaneBackend() = default;

  ArcaneBackend::ArcaneBackend(const tfel::utilities::DataMap& opts) {
    // no option
    auto validator = tfel::utilities::DataMapValidator{};
    validator.validate(opts);
  }  // end of ArcaneBackend

  std::string ArcaneBackend::getName() const { return "arcane"; }

  std::string ArcaneBackend::getSourceFileExtension() const {
    // the source is compiled by the CUDA compiler, as with the `cuda`
    // backend
    return "cu";
  }  // end of getSourceFileExtension

  void ArcaneBackend::writeSpecificIncludesInSourceFile(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd) const {
    CxxProgrammingModelBackendBase::writeSpecificIncludesInSourceFile(os, i,
                                                                      mpd);
    os << "#include\"MFront/GenericParallel/ArcaneUtilities.hxx\"\n\n";
  }  // end of writeSpecificIncludesInSourceFile

  std::string ArcaneBackend::getKernelSpecifiers() const {
    // the kernels are called by the commands of Arccore, on the host or on
    // the device
    return " ARCCORE_HOST_DEVICE";
  }  // end of getKernelSpecifiers

  void ArcaneBackend::writeExecutionPreparation(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {
    // the queue of the execution policy gives the execution policy of
    // Arccore, which is only known at runtime: it is required
    os << "if (mfront_backend_policy.policy == nullptr) {\n";
    writeBatchError(os, -5,
                    "\"no execution policy given: the arcane backend "
                    "requires a queue\"");
    os << "}\n"
       << "auto mfront_arcane_queue = mfront_backend_policy.policy->queue;\n"
       << "if (const auto* const mfront_queue_error = "
       << "::mfront::generic_parallel::arcane::checkQueue("
       << "mfront_arcane_queue); mfront_queue_error != nullptr) {\n";
    writeBatchError(os, -5, "mfront_queue_error");
    os << "}\n";
  }  // end of writeExecutionPreparation

  void ArcaneBackend::writeKernelCall(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd) const {
    const auto types = i.getTypesDescription();
    if (this->isBoundsCheckingRequired(mpd)) {
      // the statistics are written by the kernels and read on the host. The
      // blocks of unified memory, by allocator of Arccore, are kept until
      // the end of the process.
      writeStatisticsDeclaration(
          os, "::Arcane::IMemoryAllocator*",
          "::mfront::generic_parallel::arcane::getUnifiedMemoryAllocator()",
          "return ::mfront::generic_parallel::arcane::allocate<"
          "mfront_gpmp_Statistics>(::mfront::generic_parallel::arcane::"
          "getUnifiedMemoryAllocator());");
    }
    this->writeKernels(os, i, mpd, fd);
    auto write_loop = [&os, &types](const std::string_view npoints,
                                    const bool handleStrides) {
      const auto kernel = std::string{
          handleStrides ? "mfront_kernel" : "mfront_kernel_without_strides"};
      // stack size of the threads, as for the behaviours
      os << "const auto mfront_stack_size_error = "
         << "::mfront::generic_parallel::arcane::ensureStackSize("
         << "mfront_arcane_queue, " << kernel
         << ", mfront_preserve_stack_size);\n"
         << "if (!mfront_stack_size_error.empty()) {\n";
      writeBatchError(os, -7, "mfront_stack_size_error");
      os << "}\n";
      os << "// loop over the points\n"
         << "::mfront::generic_parallel::arcane::launch(mfront_arcane_queue, "
         << "static_cast<" << types.integer_type << ">(" << npoints << "), "
         << kernel
         << ");\n"
         // an asynchronous evaluation returns once the commands are launched
         << "if (!mfront_asynchronous) {\n"
         << "mfront_arcane_queue.barrier();\n"
         << "}\n";
    };
    this->writeKernelSelection(os, i, mpd, write_loop);
  }  // end of writeKernelCall

  bool ArcaneBackend::mergesPointStatusesInKernels() const {
    return true;
  }  // end of mergesPointStatusesInKernels

  std::string ArcaneBackend::getAtomicIncrement(
      const std::string_view c) const {
    return "::mfront::generic_parallel::arcane::atomicAdd(&(" + std::string{c} +
           "), 1ull)";
  }  // end of getAtomicIncrement

  std::string ArcaneBackend::getAtomicMinimum(const std::string_view c,
                                              const std::string_view v) const {
    return "::mfront::generic_parallel::arcane::atomicMin(&(" + std::string{c} +
           "), " + std::string{v} + ")";
  }  // end of getAtomicMinimum

  std::string ArcaneBackend::getMemorySpace() const {
    // the device of the queue is only known at runtime
    return "mfront::generic_parallel::device_memory_space";
  }  // end of getMemorySpace

  std::string ArcaneBackend::getExecutionPolicyType() const {
    return "::mfront::generic_parallel::ArcaneBackendExecutionPolicy";
  }  // end of getExecutionPolicyType

  std::string ArcaneBackend::getExecutionPolicyHeader() const {
    return "MFront/GenericParallel/ArcaneBackendExecutionPolicy.hxx";
  }  // end of getExecutionPolicyHeader

  ArcaneBackend::~ArcaneBackend() noexcept = default;

}  // end of namespace mfront::generic_parallel::material_property
