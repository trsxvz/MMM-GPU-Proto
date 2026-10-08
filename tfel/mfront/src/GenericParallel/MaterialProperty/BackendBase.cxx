/*!
 * \file   BackendBase.cxx
 * \brief  This file implements the `Backend` class
 * \author Thomas Helfer
 * \date   07/04/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <string>
#include <ostream>
#include "MFront/VariableDescription.hxx"
#include "MFront/MaterialPropertyDescription.hxx"
#include "MFront/MaterialPropertyParametersHandler.hxx"
#include "MFront/GenericParallelMaterialPropertyInterface.hxx"
#include "MFront/GenericParallel/MaterialProperty/BackendBase.hxx"

namespace mfront::generic_parallel::material_property {

  BackendBase::BackendBase() = default;

  bool BackendBase::handlesDataTransfer() const {
    return false;
  }  // end of handlesDataTransfer

  std::size_t BackendBase::getStatusesArraySize(
      const MaterialPropertyDescription& mpd) {
    return getInvalidResultStatusIndex(mpd) + 1;
  }  // end of getStatusesArraySize

  std::size_t BackendBase::getInvalidResultStatusIndex(
      const MaterialPropertyDescription& mpd) {
    return 2 * (mpd.inputs.size() + 1);
  }  // end of getInvalidResultStatusIndex

  void BackendBase::writeBoundsStatusesDeclaration(
      std::ostream& os,
      const MaterialPropertyDescription& mpd,
      const std::string_view device_type,
      const std::string_view device,
      const std::string_view allocation) {
    const auto n = getStatusesArraySize(mpd);
    const auto pool = "mfront::generic_parallel::StatusBlockPool<int, " +
                      std::string{device_type} + ">";
    // statuses given by the calling code (asynchronous evaluation), or block
    // of the pool, initialized here
    os << "static auto mfront_bounds_statuses_blocks = " << pool << "{};\n"
       << "auto mfront_bounds_statuses_block = std::optional<" << pool
       << "::Block>{};\n"
       << "auto* mfront_bounds_statuses = mfront_options.statuses;\n"
       << "if (mfront_bounds_statuses == nullptr) {\n"
       << "mfront_bounds_statuses_block.emplace("
       << "mfront_bounds_statuses_blocks.acquire(" << device << ", [&]{\n"
       << "constexpr auto mfront_n = std::size_t{" << n << "};\n"
       << allocation << "\n"
       << "}));\n"
       << "mfront_bounds_statuses = mfront_bounds_statuses_block->get();\n"
       << "if (mfront_bounds_statuses == nullptr) {\n"
       << "mfront_output_status->status = -6;\n"
       << "mfront_report(\"unable to allocate the statuses of the bounds "
       << "checks\");\n"
       << "return;\n"
       << "}\n"
       << "std::fill(mfront_bounds_statuses, mfront_bounds_statuses + " << n
       << ", 0);\n"
       << "}\n";
  }  // end of writeBoundsStatusesDeclaration

  void BackendBase::writeErrorReporting(std::ostream& os) {
    os << "[[maybe_unused]] auto mfront_report = "
       << "[&mfront_output_status](const std::string& "
       << "mfront_error_message){\n"
       << "if(mfront_error_message.empty()){\n"
       << "return;\n"
       << "}\n"
       << "std::strncpy(mfront_output_status->msg,"
       << "mfront_error_message.c_str(),511);\n"
       << "mfront_output_status->msg[511]='\\0';\n"
       << "};\n";
  }  // end of writeErrorReporting

  std::size_t BackendBase::getNumberOfStatuses(
      const MaterialPropertyDescription& mpd) const {
    return isBoundsCheckingRequired(mpd) ? getStatusesArraySize(mpd) : 0;
  }  // end of getNumberOfStatuses

  void BackendBase::writeStatusesReportFunction(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd) const {
    const auto name = i.getFunctionName(mpd);
    const auto types = i.getTypesDescription();
    os << "MFRONT_SHAREDOBJ void " << name << "_generic_parallel_report("
       << types.output_status_type << "* const mfront_output_status,\n"
       << "[[maybe_unused]] const int* const mfront_bounds_statuses,\n"
       << "[[maybe_unused]] const " << types.out_of_bounds_policy_type
       << " mfront_out_of_bounds_policy){\n";
    writeErrorReporting(os);
    os << "mfront_output_status->status = 0;\n"
       << "mfront_output_status->bounds_status = 0;\n"
       << "mfront_output_status->c_error_number = 0;\n";
    this->writeBoundsStatusesAnalysis(os, i, mpd);
    os << "} // end of " << name << "_generic_parallel_report\n\n";
  }  // end of writeStatusesReportFunction

  void BackendBase::writeDeviceStackSizeSetting(std::ostream& os,
                                                const std::size_t stack_size) {
    if (stack_size == 0) {
      return;
    }
    os << "if (!mfront::generic_parallel::setMinimalCUDAStackSize("
       << stack_size << ", mfront_options.preserve_stack_size != 0)) {\n"
       << "mfront_output_status->status = -6;\n"
       << "mfront_report("
       << "mfront::generic_parallel::stack_size_error_message);\n"
       << "return;\n"
       << "}\n";
  }  // end of writeDeviceStackSizeSetting

  void BackendBase::writeInvalidResultCheck(
      std::ostream& os, const MaterialPropertyDescription& mpd) const {
    const auto v =
        useQuantities(mpd) ? mpd.output.name + ".getValue()" : mpd.output.name;
    os << "if(!tfel::math::ieee754::isfinite(" << v << ")){\n"
       << "mfront_bounds_statuses[" << getInvalidResultStatusIndex(mpd)
       << "] = 1;\n"
       << "}\n";
  }  // end of writeInvalidResultCheck

  void BackendBase::writeCDeclarations(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd) const {
    const auto name = i.getFunctionName(mpd);
    const auto types = i.getTypesDescription();
    const auto output =
        this->handlesDataTransfer() ? "mfront_output_host" : "mfront_output";
    const auto args =
        this->handlesDataTransfer() ? "mfront_args_host" : "mfront_args";
    auto document_arguments = [&os, args] {
      os << " * \\param[in] " << args
         << ": array of pointers to the values of the arguments\n";
    };
    auto document_last_arguments = [&os] {
      os << " * \\param[in] mfront_nargs: number of arguments\n"
         << " * \\param[in] mfront_npoints: number of points in which the "
         << "material property is computed\n"
         << " * \\param[in] mfront_out_of_bounds_policy: out of bounds "
         << "policy\n"
         << " * \\param[in] mfront_execution_options: options of the "
         << "execution (queue, statuses of an asynchronous evaluation, stack "
         << "size), or null for the default options\n";
    };
    os << "/*!\n"
       << " * \\brief compute the value of '" << mpd.className
       << "' material property on a set of points\n"
       << " *\n"
       << " * \\param[out] mfront_output_status: output status\n"
       << " * \\param[out] " << output << ": output values\n"
       << " * \\param[in] mfront_output_stride: output stride.\n"
       << " * If zero, only the first element of the output is computed.\n"
       << " * A zero value is only meaningful if all inputs are uniform\n";
    document_arguments();
    os << " * \\param[in] mfront_args_strides: array of integers giving the "
       << "stride of each argument. If zero, the argument is uniform\n";
    document_last_arguments();
    os << " *\n"
       << " * \\note `mfront_npoints` is also the size of the arrays giving\n"
       << " *   the output and the inputs when non uniform\n"
       << " */\n"
       << "MFRONT_SHAREDOBJ void " << name << "(\n"
       << types.output_status_type << "* const,\n"
       << types.real_type << "* const,\n"
       << "const " << types.integer_type << ",\n"
       << "const " << types.real_type << "* const* const,\n"
       << "const " << types.integer_type << "*,\n"
       << "const " << types.integer_type << ",\n"
       << "const " << types.integer_type << ",\n"
       << "const " << types.out_of_bounds_policy_type << ",\n"
       << "const mfront_gpmp_ExecutionOptions* const);\n\n"
       << "/*!\n"
       << " * \\brief compute the value of '" << mpd.className
       << "' material property on a set of points, the output and the "
       << "arguments being contiguous arrays\n"
       << " *\n"
       << " * \\param[out] mfront_output_status: output status\n"
       << " * \\param[out] " << output << ": output values\n";
    document_arguments();
    document_last_arguments();
    os << " *\n"
       << " * \\note `mfront_npoints` is also the size of the arrays giving\n"
       << " *   the output and the inputs\n"
       << " */\n"
       << "MFRONT_SHAREDOBJ void " << name << "2(\n"
       << types.output_status_type << "* const,\n"
       << types.real_type << "* const,\n"
       << "const " << types.real_type << "* const* const,\n"
       << "const " << types.integer_type << ",\n"
       << "const " << types.integer_type << ",\n"
       << "const " << types.out_of_bounds_policy_type << ",\n"
       << "const mfront_gpmp_ExecutionOptions* const);\n\n"
       << "/*!\n"
       << " * \\brief analyse the statuses written by an asynchronous "
       << "evaluation of '" << mpd.className << "' material property, once "
       << "the computations are completed, and report the result as a "
       << "synchronous evaluation does\n"
       << " *\n"
       << " * \\param[out] mfront_output_status: output status\n"
       << " * \\param[in] mfront_bounds_statuses: statuses given by the "
       << "execution options, read on the host\n"
       << " * \\param[in] mfront_out_of_bounds_policy: out of bounds policy "
       << "of the evaluation\n"
       << " */\n"
       << "MFRONT_SHAREDOBJ void " << name << "_generic_parallel_report(\n"
       << types.output_status_type << "* const,\n"
       << "const int* const,\n"
       << "const " << types.out_of_bounds_policy_type << ");\n\n";
  }  // end of writeDeclarations

  void BackendBase::writeProxyObjects(
      std::ostream&,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {}  // end of writeProxyObject

  void BackendBase::writePhysicalBounds(std::ostream& os,
                                        const VariableDescription& v,
                                        const size_t i) const {
    if (!v.hasPhysicalBounds()) {
      return;
    }
    const auto& b = v.getPhysicalBounds();
    if ((b.boundsType == VariableBoundsDescription::LOWER) ||
        (b.boundsType == VariableBoundsDescription::LOWERANDUPPER)) {
      os << "if(" << v.name << " < " << v.type << "(" << b.lowerBound << ")){\n"
         << "mfront_bounds_statuses[" << 2 * i << "] = -1;\n"
         << "return;\n"
         << "}\n";
    }
    if ((b.boundsType == VariableBoundsDescription::UPPER) ||
        (b.boundsType == VariableBoundsDescription::LOWERANDUPPER)) {
      os << "if(" << v.name << " > " << v.type << "(" << b.upperBound << ")){\n"
         << "mfront_bounds_statuses[" << 2 * i << "] = 1;\n"
         << "return;\n"
         << "}\n";
    }
  }  // end of writePhysicalBounds

  void BackendBase::writeBounds(std::ostream& os,
                                const std::string& prefix,
                                const VariableDescription& v,
                                const size_t i) const {
    if (!v.hasBounds()) {
      return;
    }
    const auto& b = v.getBounds();
    if ((b.boundsType == VariableBoundsDescription::LOWER) ||
        (b.boundsType == VariableBoundsDescription::LOWERANDUPPER)) {
      os << "if(" << v.name << " < " << v.type << "(" << b.lowerBound << ")){\n"
         << "mfront_bounds_statuses[" << 2 * i + 1 << "] = -1;\n"
         << "if (mfront_out_of_bounds_policy==" << prefix
         << "_STRICT_POLICY){\n"
         << "return;\n"
         << "}\n"
         << "}\n";
    }
    if ((b.boundsType == VariableBoundsDescription::UPPER) ||
        (b.boundsType == VariableBoundsDescription::LOWERANDUPPER)) {
      os << "if(" << v.name << " > " << v.type << "(" << b.upperBound << ")){\n"
         << "mfront_bounds_statuses[" << 2 * i + 1 << "] = 1;\n"
         << "if (mfront_out_of_bounds_policy==" << prefix
         << "_STRICT_POLICY){\n"
         << "return;\n"
         << "}\n"
         << "}\n";
    }
  }  // end of writeBounds

  void BackendBase::writeBoundsStatusesAnalysis(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd) const {
    if (areRuntimeChecksDisabled(mpd)) {
      return;
    }
    os << "if(mfront_bounds_statuses[" << getInvalidResultStatusIndex(mpd)
       << "] != 0){\n"
       << "mfront_report(\"invalid result at at least one point\");\n"
       << "mfront_output_status->status = -4;\n"
       << "return;\n"
       << "}\n";
    if (!((hasPhysicalBounds(mpd.inputs)) || (hasBounds(mpd.inputs)) ||
          (mpd.output.hasPhysicalBounds()) || (mpd.output.hasBounds()))) {
      return;
    }
    const auto iucname = i.getInterfaceNameInUpperCase();
    const auto prefix = i.getOutOfBoundsPolicyEnumerationPrefix();
    os << "#ifndef NO_" << iucname << "_BOUNDS_CHECK\n";
    auto write_physical_bounds_violation = [&os](const VariableDescription& v,
                                                 std::size_t idx) {
      os << "if(mfront_bounds_statuses[" << 2 * idx << "] > 0){\n"
         << "mfront_report(\"" << v.name
         << " is over its physical upper bound at at least one "
         << "point\");\n"
         << "mfront_output_status->status = -1;\n"
         << "mfront_output_status->bounds_status = -" << idx + 1 << ";\n"
         << "return;\n"
         << "}\n";
      os << "if(mfront_bounds_statuses[" << 2 * idx << "] < 0){\n"
         << "mfront_report(\"" << v.name
         << " is below its physical lower bound at at least one "
         << "point\");\n"
         << "mfront_output_status->status = -1;\n"
         << "mfront_output_status->bounds_status = -" << idx + 1 << ";\n"
         << "return;\n"
         << "}\n";
    };
    for (decltype(mpd.inputs.size()) idx = 0; idx != mpd.inputs.size(); ++idx) {
      const auto& v = mpd.inputs[idx];
      if (!v.hasPhysicalBounds()) {
        continue;
      }
      write_physical_bounds_violation(v, idx);
    }
    if (mpd.output.hasPhysicalBounds()) {
      write_physical_bounds_violation(mpd.output, mpd.inputs.size());
    }
    if ((hasBounds(mpd.inputs)) || (mpd.output.hasBounds())) {
      auto write_bounds_violation = [&os, &prefix](const VariableDescription& v,
                                                   std::size_t idx) {
        os << "if(mfront_bounds_statuses[" << 2 * idx + 1 << "] > 0){\n"
           << "mfront_report(\"" << v.name
           << " is over its upper bound at at least one point\");\n"
           << "mfront_output_status->bounds_status = " << idx + 1 << ";\n"
           << "if(mfront_out_of_bounds_policy==" << prefix
           << "_STRICT_POLICY){\n"
           << "mfront_output_status->status = -1;\n"
           << "return;\n"
           << "}\n"
           << "mfront_output_status->status = 1;\n"
           << "}\n";
        os << "if(mfront_bounds_statuses[" << 2 * idx + 1 << "] < 0){\n"
           << "mfront_report(\"" << v.name
           << " is below its lower bound at at least one point\");\n"
           << "mfront_output_status->bounds_status = " << idx + 1 << ";\n"
           << "if(mfront_out_of_bounds_policy==" << prefix
           << "_STRICT_POLICY){\n"
           << "mfront_output_status->status = -1;\n"
           << "return;\n"
           << "}\n"
           << "mfront_output_status->status = 1;\n"
           << "}\n";
      };
      // as with the `generic` interface, the result must be used with
      // caution (status 1) if a bound is violated with the warning policy
      os << "if((mfront_out_of_bounds_policy==" << prefix
         << "_STRICT_POLICY)||\n"
         << "(mfront_out_of_bounds_policy==" << prefix
         << "_WARNING_POLICY)){\n";
      for (decltype(mpd.inputs.size()) idx = 0; idx != mpd.inputs.size();
           ++idx) {
        const auto& v = mpd.inputs[idx];
        if (!v.hasBounds()) {
          continue;
        }
        write_bounds_violation(v, idx);
      }
      if (mpd.output.hasBounds()) {
        write_bounds_violation(mpd.output, mpd.inputs.size());
      }
      os << "}\n";
    }
    os << "#endif /* NO_" << iucname << "_BOUNDS_CHECK */\n";
  }

  bool BackendBase::isBoundsCheckingRequired(
      const MaterialPropertyDescription& mpd) {
    // the array of statuses is also used to report invalid results
    return !areRuntimeChecksDisabled(mpd);
  }  // end of isBoundsCheckingRequired

  void BackendBase::writeCFunctionBeginning(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const bool treatStrides) const {
    const auto name = i.getFunctionName(mpd);
    const auto nname = i.getInterfaceInternalNamespace();
    const auto types = i.getTypesDescription();
    const auto& params = mpd.parameters;
    // signature
    os << "MFRONT_SHAREDOBJ void " << name;
    if (!treatStrides) {
      os << '2';
    }
    os << "(" << types.output_status_type << "* const mfront_output_status,\n";
    if (this->handlesDataTransfer()) {
      os << types.real_type << "* const mfront_output_host";
    } else {
      os << types.real_type << "* const mfront_output";
    }
    if (treatStrides) {
      os << ",\nconst " << types.integer_type << " mfront_output_stride";
    }
    os << ",\nconst " << types.real_type << "* const * const";
    if (!mpd.inputs.empty()) {
      os << (this->handlesDataTransfer() ? " mfront_args_host"
                                         : " mfront_args");
    }
    if (treatStrides) {
      os << ",\nconst " << types.integer_type << "* const";
      if (!mpd.inputs.empty()) {
        os << " mfront_args_strides";
      }
    }
    // the policy is passed to the kernels with the array of statuses, even
    // if the material property has no bounds
    os << ",\nconst " << types.integer_type << " mfront_nargs"
       << ",\nconst " << types.integer_type << " mfront_npoints"
       << ",\n[[maybe_unused]] const " << types.out_of_bounds_policy_type
       << " mfront_out_of_bounds_policy"
       << ",\nconst mfront_gpmp_ExecutionOptions* const "
       << "mfront_execution_options";
    os << ")\n{\n";
    // execution options, default options if none is given
    os << "[[maybe_unused]] const auto mfront_options = "
       << "(mfront_execution_options != nullptr) ? *mfront_execution_options "
       << ": mfront_gpmp_ExecutionOptions{};\n"
       << "[[maybe_unused]] const bool mfront_asynchronous = "
       << "mfront_options.statuses != nullptr;\n";
    writeErrorReporting(os);
    if (!areRuntimeChecksDisabled(mpd)) {
      os << "auto mfront_check_positive_or_null = "
         << "[]<std::integral mfront_IntegerType>"
         << "(const mfront_IntegerType mfront_integral_value){\n"
         << "if constexpr (std::is_signed_v<mfront_IntegerType>){\n"
         << "if(mfront_integral_value < 0){\n"
         << "return false;\n"
         << "}\n"
         << "}\n"
         << "return true;\n"
         << "};\n";
    }
    // errno is restored before returning, even if the runtime checks are
    // disabled (checks of the strides)
    os << "[[maybe_unused]] const int mfront_errno_old = errno;\n";
    os << "mfront_output_status->status = 0;\n"
       << "mfront_output_status->bounds_status = 0;\n"
       << "mfront_output_status->c_error_number = 0;\n";
    // writes an error and returns
    auto write_error = [&os](const int status, const std::string& message) {
      os << "mfront_output_status->status = " << status << ";\n"
         << "mfront_report(" << message << ");\n"
         << "errno = mfront_errno_old;\n"
         << "return;\n"
         << "}\n";
    };
    if (!areRuntimeChecksDisabled(mpd)) {
      os << "errno = 0;\n";
      os << "if(!mfront_check_positive_or_null(mfront_npoints)){\n";
      write_error(
          -5,
          "\"invalid number of points (\"+std::to_string(mfront_npoints)"
          "+\" given)\"");
      os << "if(mfront_nargs != " << mpd.inputs.size() << "){\n";
      write_error(
          -5,
          "\"invalid number of arguments (\"+std::to_string(mfront_nargs)"
          "+\" given, " +
              std::to_string(mpd.inputs.size()) + " expected)\"");
      if (treatStrides) {
        os << "if(!mfront_check_positive_or_null(mfront_output_stride)){\n";
        write_error(-5,
                    "\"negative output stride (\"+std::to_string("
                    "mfront_output_stride)+\" given)\"");
        for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
          const auto stride =
              "mfront_args_strides[" + std::to_string(idx) + "]";
          os << "if(!mfront_check_positive_or_null(" << stride << ")){\n";
          write_error(-5, "\"negative stride given for variable '" +
                              mpd.inputs[idx].getExternalName() +
                              "' (\"+std::to_string(" + stride +
                              ")+\" given)\"");
        }
      }
      if ((!areParametersTreatedAsStaticVariables(mpd)) && (!params.empty())) {
        const auto hn = getMaterialPropertyParametersHandlerClassName(name);
        const auto handler = nname + "::" + hn + "::get" + hn + "()";
        os << "if(!" << handler << ".ok){\n";
        write_error(-6, handler + ".msg");
      }
    }
    if (treatStrides) {
      os << "if((mfront_output_stride == 0) && (mfront_npoints != 1)){\n";
      write_error(-7,
                  "\"if the output is uniform, the number of points must be "
                  "equal to one (\" + std::to_string(mfront_npoints) + \" "
                  "given)\"");
    }
    os << "if(mfront_npoints <= 0){\n"
       << "return;\n"
       << "}\n";
    if ((treatStrides) && (!mpd.inputs.empty())) {
      os << "auto mfront_strides = std::array<" << types.integer_type << ", "
         << mpd.inputs.size() << ">{};\n"
         << "tfel::fsalgo::copy<" << mpd.inputs.size() << ">::exe("
         << "mfront_args_strides, mfront_strides.begin());\n"
         << "const bool mfront_areAllArgumentsUniform = std::all_of("
         << "mfront_strides.begin(), mfront_strides.end(), "
         << "[](const " << types.integer_type
         << " mfront_stride){return mfront_stride == 0;});\n"
         << "if((mfront_output_stride == 0) && "
         << "(!mfront_areAllArgumentsUniform)){\n";
      write_error(-7,
                  "\"if the output is uniform, all the arguments must be "
                  "uniform (i.e. their strides must be null)\"");
    }
  }  // end of writeCFunctionBeginning

  void BackendBase::writeKernelBody(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const bool handleStrides) const {
    const auto iucname = i.getInterfaceNameInUpperCase();
    const auto prefix = i.getOutOfBoundsPolicyEnumerationPrefix();
    const auto qt = useQuantities(mpd);
    // parameters
    for (const auto& p : mpd.parameters) {
      if (areParametersTreatedAsStaticVariables(mpd)) {
        const auto pv =
            p.getAttribute<double>(VariableDescription::defaultValue);
        os << "constexpr auto " << p.name << " = "
           << (qt ? p.type : std::string{"double"}) << "(" << pv << ");\n";
      } else {
        os << "const auto " << p.name << " = "
           << (qt ? p.type : std::string{"real"}) << "(mfront_parameters."
           << p.name << ");\n";
      }
    }
    // arguments
    for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
      const auto& v = mpd.inputs[idx];
      const auto arg = "mfront_arg" + std::to_string(idx);
      const auto value = handleStrides
                             ? arg + "[mfront_strides[" + std::to_string(idx) +
                                   "] * mfront_idx]"
                             : arg + "[mfront_idx]";
      os << "const auto " << v.name << " = "
         << (qt ? v.type + "(" + value + ")" : value) << ";\n";
    }
    // output
    os << "auto " << mpd.output.name << " = " << mpd.output.type << "{};\n";
    const auto checks = !areRuntimeChecksDisabled(mpd);
    // bounds of the arguments: physical bounds first
    if ((checks) &&
        ((hasPhysicalBounds(mpd.inputs)) || (hasBounds(mpd.inputs)))) {
      os << "#ifndef NO_" << iucname << "_BOUNDS_CHECK\n";
      for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
        this->writePhysicalBounds(os, mpd.inputs[idx], idx);
      }
      for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
        this->writeBounds(os, prefix, mpd.inputs[idx], idx);
      }
      os << "#endif /* NO_" << iucname << "_BOUNDS_CHECK */\n";
    }
    os << mpd.f.body;
    // bounds of the output and invalid results
    if (checks) {
      if ((mpd.output.hasPhysicalBounds()) || (mpd.output.hasBounds())) {
        os << "#ifndef NO_" << iucname << "_BOUNDS_CHECK\n";
        this->writePhysicalBounds(os, mpd.output, mpd.inputs.size());
        this->writeBounds(os, prefix, mpd.output, mpd.inputs.size());
        os << "#endif /* NO_" << iucname << "_BOUNDS_CHECK */\n";
      }
      this->writeInvalidResultCheck(os, mpd);
    }
    os << (handleStrides ? "mfront_output[mfront_output_stride * mfront_idx] = "
                         : "mfront_output[mfront_idx] = ")
       << mpd.output.name << (qt ? ".getValue();\n" : ";\n");
  }  // end of writeKernelBody

  void BackendBase::writeKernelSelection(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd,
      const bool treatStrides,
      const std::function<void(std::string_view, bool)>& write_call) const {
    if (!treatStrides) {
      write_call("mfront_npoints", false);
      return;
    }
    if (mpd.inputs.empty()) {
      os << "if(mfront_output_stride == 0){\n";
      write_call("1", false);
      os << "} else if(mfront_output_stride == 1){\n";
      write_call("mfront_npoints", false);
      os << "} else {\n";
      write_call("mfront_npoints", true);
      os << "}\n";
      return;
    }
    const auto types = i.getTypesDescription();
    os << "if(mfront_areAllArgumentsUniform && "
       << "(mfront_output_stride == 0)){\n";
    write_call("1", false);
    os << "} else {\n"
       << "const bool mfront_areAllArgumentsStrideOne = std::all_of("
       << "mfront_strides.begin(), mfront_strides.end(), "
       << "[](const " << types.integer_type << " mfront_stride){"
       << "return mfront_stride == 1;});\n"
       << "if(mfront_areAllArgumentsStrideOne && (mfront_output_stride == "
          "1)){\n";
    write_call("mfront_npoints", false);
    os << "} else {\n";
    write_call("mfront_npoints", true);
    os << "}\n"
       << "}\n";
  }  // end of writeKernelSelection

  void BackendBase::writeResultsAnalysis(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd) const {
    if (areRuntimeChecksDisabled(mpd)) {
      return;
    }
    // errno can't be used to detect errors in the kernels: it is
    // thread-local, does not exist on GPUs, and may be modified by the
    // runtimes (for instance at the initialisation of TBB). Invalid results
    // are detected by the kernels.
    os << "errno = mfront_errno_old;\n"
       << "if (!mfront_asynchronous) {\n"
       << i.getFunctionName(mpd) << "_generic_parallel_report("
       << "mfront_output_status, mfront_bounds_statuses, "
       << "mfront_out_of_bounds_policy);\n"
       << "}\n";
  }  // end of writeResultsAnalysis

  BackendBase::~BackendBase() noexcept = default;

}  // end of namespace mfront::generic_parallel::material_property
