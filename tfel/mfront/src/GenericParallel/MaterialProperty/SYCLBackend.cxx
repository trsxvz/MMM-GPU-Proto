/*!
 * \file   SYCLBackend.cxx
 * \brief  This file implements the `SYCLBackend` class
 * \author Thomas Helfer
 * \date   18/05/2026
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
#include "MFront/GenericParallel/MaterialProperty/SYCLBackend.hxx"

namespace mfront::generic_parallel::material_property {

  SYCLBackend::SYCLBackend() {}  // end of SYCLBackend

  SYCLBackend::SYCLBackend(const tfel::utilities::DataMap& opts) {
    auto validator = tfel::utilities::DataMapValidator{};
    declareDataTransferPolicyOption(validator);
    declareDeviceStackSizeOption(validator);
    validator.validate(opts);
    this->data_transfer_policy = areDataTransfersAutomatic(opts)
                                     ? DataTransferPolicy::AUTOMATIC
                                     : DataTransferPolicy::NONE;
    this->device_stack_size = getDeviceStackSize(opts);
  }  // end of SYCLBackend

  std::string SYCLBackend::getName() const { return "sycl"; }  // end of getName

  void SYCLBackend::writeSpecificIncludesInSourceFile(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd) const {
    CxxProgrammingModelBackendBase::writeSpecificIncludesInSourceFile(os, i,
                                                                      mpd);
    os << "#include <sycl/sycl.hpp>\n";
    if (this->device_stack_size != 0) {
      os << "#include\"MFront/GenericParallel/DeviceStackSize.hxx\"\n";
    }
    os << "\n";

    if (this->handlesDataTransfer()) {
      os << "namespace mfront::gpmp::sycl{\n\n"
         << "template<typename T>\n"
         << "struct [[nodiscard]] DataBuffer{\n"
         << "DataBuffer(::sycl::queue& q)\n"
         << " : queue(q)\n"
         << "{}\n"
         << "DataBuffer(DataBuffer&&) = delete;\n"
         << "DataBuffer(const DataBuffer&) = delete;\n"
         << "DataBuffer& operator=(DataBuffer&&) = delete;\n"
         << "DataBuffer& operator=(const DataBuffer&) = delete;\n"
         << "[[nodiscard]] ::sycl::event acquire(T* const src, "
         << "const std::size_t s, const bool do_copy){\n"
         << "const bool is_queue_device = !this->queue.get_device().is_cpu();\n"
         << "const auto ptr_type = ::sycl::get_pointer_type(src, "
         << "this->queue.get_context());\n"
         << "bool may_require_device_allocation = \n"
         << "                     (ptr_type == ::sycl::usm::alloc::host) ||\n"
         << "                     (ptr_type == ::sycl::usm::alloc::unknown);\n"
         << "this->cleanup();\n"
         << "this->host_ptr = src;\n"
         << "if (may_require_device_allocation && is_queue_device) {\n"
         << "  this->device_ptr = "
         << "  ::sycl::malloc_device<std::remove_const_t<T>>(s, this->queue);\n"
         << "  this->ptr = this->device_ptr;\n"
         << "  if (do_copy) {\n"
         << "    return this->queue.copy(src, this->device_ptr, s);\n"
         << "  }\n"
         << "  return ::sycl::event{};\n"
         << "}\n"
         << "this->device_ptr = nullptr;\n"
         << "this->ptr = src;\n"
         << "return ::sycl::event{};\n"
         << "}\n"
         << "[[nodiscard]] T* data() noexcept{\n"
         << "return this->ptr;\n"
         << "}\n"
         << "[[nodiscard]] bool requiresSynchronization() const "
         << "requires(!std::is_const_v<T>){ "
         << "  return this->device_ptr != nullptr;\n"
         << "}\n"
         << "[[nodiscard]] bool holdsCopy() const noexcept {\n"
         << "  return this->device_ptr != nullptr;\n"
         << "}\n"
         << "~DataBuffer() noexcept{\n"
         << "try{\n"
         << "this->cleanup();\n"
         << "} catch(...){}\n"
         << "}\n"
         << "void cleanup(){\n"
         << "if(this->device_ptr != nullptr) {\n"
         << "::sycl::free(this->device_ptr, this->queue);\n"
         << "}\n"
         << "this->host_ptr = nullptr;\n"
         << "this->device_ptr = nullptr;\n"
         << "this->ptr = nullptr;\n"
         << "}\n"
         << "private:\n"
         << "::sycl::queue queue;\n"
         << "T* ptr = nullptr;\n"
         << "T* host_ptr = nullptr;\n"
         << "std::remove_const_t<T>* device_ptr = nullptr;\n"
         << "};\n\n"
         << "} // end of namespace mfront::gpmp::sycl\n";
    }

  }  // end of writeSpecificIncludesInSourceFile

  void SYCLBackend::writeExecutionPreparation(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {
    // SYCL has no default queue: the queue gives the device and the context
    // of the unified shared memory of the calling code
    os << "if (mfront_options.queue == nullptr) {\n"
       << "mfront_output_status->status = -5;\n"
       << "mfront_report(\"no SYCL queue given\");\n"
       << "return;\n"
       << "}\n"
       << "auto* const mfront_queue = "
       << "static_cast<::sycl::queue*>(mfront_options.queue);\n";
  }  // end of writeExecutionPreparation

  void SYCLBackend::writeKernelCall(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd,
      const bool treatStrides) const {
    const auto types = i.getTypesDescription();
    const auto requiresBoundsCheck = this->isBoundsCheckingRequired(mpd);
    writeDeviceStackSizeSetting(os, this->device_stack_size);
    if (requiresBoundsCheck) {
      // the statuses are held by a block of shared memory, identified by
      // the device and the context of the queue
      writeBoundsStatusesDeclaration(
          os, mpd, "std::pair<sycl::device, sycl::context>",
          "{mfront_queue->get_device(), mfront_queue->get_context()}",
          "return sycl::malloc_shared<int>(mfront_n, *mfront_queue);");
    }
    this->writeKernels(os, i, mpd, fd, treatStrides);
    auto write_kernel_call = [&os, this, &types](const std::string_view npoints,
                                                 const bool handleStrides) {
      os << "// loop over the points\n";
      os << "auto mfront_sycl_kernel_event = ";
      os << "mfront_queue->parallel_for(" << npoints;
      if (this->handlesDataTransfer()) {
        os << ", mfront_sycl_copy_events";
      }
      os << ", mfront_kernel";
      if (!handleStrides) {
        os << "_without_strides";
      }
      os << ");\n";
      if (this->handlesDataTransfer()) {
        os << "if (mfront_output_ptr.requiresSynchronization()) {\n";
        if (handleStrides) {
          if (npoints == "1") {
            os << "mfront_queue->copy<" << types.real_type
               << ">(mfront_output_ptr.data(), mfront_output_host, " << 1
               << ", mfront_sycl_kernel_event).wait();\n";
          } else {
            os << "mfront_queue->copy<" << types.real_type
               << ">(mfront_output_ptr.data(), mfront_output_host, "  //
               << "(" << npoints << " - 1) * mfront_output_stride + 1, "
               << "mfront_sycl_kernel_event).wait();\n";
          }
        } else {
          os << "mfront_queue->copy<" << types.real_type
             << ">(mfront_output_ptr.data(), mfront_output_host, " << npoints
             << ", mfront_sycl_kernel_event).wait();\n";
        }
        os << "} else if ((!mfront_asynchronous) || mfront_sycl_copies) {\n"
           << "  mfront_sycl_kernel_event.wait();\n"
           << "}\n";
      } else {
        // an asynchronous evaluation returns once the kernel is submitted
        os << "if (!mfront_asynchronous) {\n"
           << "mfront_sycl_kernel_event.wait();\n"
           << "}\n";
      }
    };
    this->writeKernelSelection(os, i, mpd, treatStrides, write_kernel_call);
  }  // end of writeKernelCall

  bool SYCLBackend::handlesDataTransfer() const {
    return this->data_transfer_policy == DataTransferPolicy::AUTOMATIC;
  }  // end of handlesDataTransfer

  void SYCLBackend::writeDataTransfersToDevice(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const bool treatStrides) const {
    const auto types = i.getTypesDescription();
    if (!this->handlesDataTransfer()) {
      return;
    }
    os << "auto mfront_sycl_copy_events = std::vector<::sycl::event>{};\n";
    if (mpd.inputs.empty()) {
      os << "mfront_sycl_copy_events.reserve(1);\n";
      os << "const " << types.real_type
         << " * const * mfront_args = nullptr;\n";
    } else {
      os << "mfront_sycl_copy_events.reserve(" << mpd.inputs.size() + 1
         << ");\n"
         << "auto mfront_args_storage = std::array<const " << types.real_type
         << "*, " << mpd.inputs.size() << ">{};\n";
      for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
        os << "auto mfront_arg" << idx << "_ptr = "
           << "::mfront::gpmp::sycl::DataBuffer<const " << types.real_type
           << ">{*mfront_queue};\n";
        if (treatStrides) {
          os << "if(mfront_args_strides[" << idx << "] ==0){\n"
             << ";\n"
             << "mfront_sycl_copy_events.push_back("
             << "mfront_arg" << idx << "_ptr.acquire("
             << "mfront_args_host[" << idx << "], 1, true));\n"
             << "} else {\n"
             << "mfront_sycl_copy_events.push_back("
             << "mfront_arg" << idx << "_ptr.acquire("
             << "mfront_args_host[" << idx << "], "
             << "(mfront_npoints - 1) * mfront_args_strides[" << idx
             << "] + 1, true));\n"
             << "}\n";
        } else {
          os << "mfront_sycl_copy_events.push_back("
             << "mfront_arg" << idx << "_ptr.acquire("
             << "mfront_args_host[" << idx << "], mfront_npoints, true));\n";
        }
        os << "mfront_args_storage[" << idx << "] = "  //
           << "mfront_arg" << idx << "_ptr.data();\n";
      }
      os << "const " << types.real_type << " * * mfront_args = "
         << "mfront_args_storage.data();\n";
    }
    os << "auto mfront_output_ptr = "
       << "::mfront::gpmp::sycl::DataBuffer<" << types.real_type
       << ">{*mfront_queue};\n";
    if (treatStrides) {
      // with a stride greater than one, the values between the values of
      // the output are copied, so that copying the output back preserves
      // them
      os << "if (mfront_output_stride == 0) {\n"
         << "mfront_sycl_copy_events.push_back("
         << "mfront_output_ptr.acquire(mfront_output_host,  1, false));\n"
         << "} else {\n"
         << "mfront_sycl_copy_events.push_back("
         << "mfront_output_ptr.acquire(mfront_output_host, "
         << "(mfront_npoints - 1) * mfront_output_stride + 1, "
         << "mfront_output_stride > 1));\n"
         << "}\n";
    } else {
      os << "mfront_sycl_copy_events.push_back("
         << "mfront_output_ptr.acquire(mfront_output_host, "
         << "mfront_npoints, false));\n";
    }
    os << "auto *mfront_output = mfront_output_ptr.data();\n";
    // the copies must be completed before the buffers are destroyed: the
    // evaluation is then synchronous
    os << "const bool mfront_sycl_copies = mfront_output_ptr.holdsCopy()";
    for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
      os << " || mfront_arg" << idx << "_ptr.holdsCopy()";
    }
    os << ";\n";
  }

  void SYCLBackend::writeDataTransfersToHost(
      std::ostream&,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&,
      const bool) const {}

  std::string SYCLBackend::getMemorySpace() const {
    // with automatic transfers, the arguments and the output which are not
    // accessible by the device are copied
    if (this->handlesDataTransfer()) {
      return "mfront::generic_parallel::any_memory_space";
    }
    return "mfront::generic_parallel::device_memory_space";
  }  // end of getMemorySpace

  SYCLBackend::~SYCLBackend() = default;

}  // end of namespace mfront::generic_parallel::material_property
