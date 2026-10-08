/*!
 * \file   BackendBase.cxx
 * \brief  This file implements the `BackendBase` class
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
#include "TFEL/Raise.hxx"
#include "MFront/VariableDescription.hxx"
#include "MFront/MaterialPropertyDescription.hxx"
#include "MFront/MaterialPropertyParametersHandler.hxx"
#include "MFront/GenericParallelMaterialPropertyInterface.hxx"
#include "MFront/GenericParallel/MaterialProperty/BackendBase.hxx"

namespace mfront::generic_parallel::material_property {

  /*
   * Identifiers of the checks performed by the kernels, in the order of the
   * analysis of the statistics (see `writeStatisticsAnalysis`): the first
   * check which failed at some point is reported. The physical bounds and
   * the bounds are numbered variable by variable, the output being after
   * the arguments, the upper bound before the lower bound.
   */

  //! \brief identifier of the check of the finiteness of the result
  static constexpr unsigned long long invalid_result_check = 0;

  /*!
   * \return the identifier of the check of a physical bound
   * \param[in] i: index of the variable
   * \param[in] upper: true for the upper bound
   */
  static unsigned long long getPhysicalBoundCheck(const std::size_t i,
                                                  const bool upper) {
    return 1 + 2 * i + (upper ? 0 : 1);
  }  // end of getPhysicalBoundCheck

  /*!
   * \return the identifier of the check of a bound
   * \param[in] mpd: material property description
   * \param[in] i: index of the variable
   * \param[in] upper: true for the upper bound
   */
  static unsigned long long getBoundCheck(
      const MaterialPropertyDescription& mpd,
      const std::size_t i,
      const bool upper) {
    return 1 + 2 * (mpd.inputs.size() + 1) + 2 * i + (upper ? 0 : 1);
  }  // end of getBoundCheck

  /*!
   * \brief write the check of the physical bounds of a variable in a
   * kernel: the evaluation fails at this point
   * \param[in] os: output stream
   * \param[in] v: variable
   * \param[in] i: index of the variable (the output is after the inputs)
   */
  static void writePhysicalBounds(std::ostream& os,
                                  const VariableDescription& v,
                                  const std::size_t i) {
    if (!v.hasPhysicalBounds()) {
      return;
    }
    const auto& b = v.getPhysicalBounds();
    auto write = [&os, &v, i](const char* const op, const auto value,
                              const bool upper) {
      os << "if(" << v.name << " " << op << " " << v.type << "(" << value
         << ")){\n"
         << "mfront_status.report(" << getPhysicalBoundCheck(i, upper)
         << "ull, true);\n"
         << "return mfront_status;\n"
         << "}\n";
    };
    if ((b.boundsType == VariableBoundsDescription::LOWER) ||
        (b.boundsType == VariableBoundsDescription::LOWERANDUPPER)) {
      write("<", b.lowerBound, false);
    }
    if ((b.boundsType == VariableBoundsDescription::UPPER) ||
        (b.boundsType == VariableBoundsDescription::LOWERANDUPPER)) {
      write(">", b.upperBound, true);
    }
  }  // end of writePhysicalBounds

  /*!
   * \brief write the check of the bounds of a variable in a kernel: the
   * evaluation fails at this point with the strict policy, and the result
   * must be used with caution with the warning policy
   * \param[in] os: output stream
   * \param[in] prefix: prefix of the out of bounds policies
   * \param[in] mpd: material property description
   * \param[in] v: variable
   * \param[in] i: index of the variable (the output is after the inputs)
   */
  static void writeBounds(std::ostream& os,
                          const std::string& prefix,
                          const MaterialPropertyDescription& mpd,
                          const VariableDescription& v,
                          const std::size_t i) {
    if (!v.hasBounds()) {
      return;
    }
    const auto& b = v.getBounds();
    auto write = [&os, &prefix, &mpd, &v, i](const char* const op,
                                             const auto value,
                                             const bool upper) {
      const auto c = getBoundCheck(mpd, i, upper);
      os << "if(" << v.name << " " << op << " " << v.type << "(" << value
         << ")){\n"
         << "if (mfront_out_of_bounds_policy == " << prefix
         << "_STRICT_POLICY){\n"
         << "mfront_status.report(" << c << "ull, true);\n"
         << "return mfront_status;\n"
         << "}\n"
         << "if (mfront_out_of_bounds_policy == " << prefix
         << "_WARNING_POLICY){\n"
         << "mfront_status.report(" << c << "ull, false);\n"
         << "}\n"
         << "}\n";
    };
    if ((b.boundsType == VariableBoundsDescription::LOWER) ||
        (b.boundsType == VariableBoundsDescription::LOWERANDUPPER)) {
      write("<", b.lowerBound, false);
    }
    if ((b.boundsType == VariableBoundsDescription::UPPER) ||
        (b.boundsType == VariableBoundsDescription::LOWERANDUPPER)) {
      write(">", b.upperBound, true);
    }
  }  // end of writeBounds

  /*!
   * \brief write the code reporting an invalid result (not finite) in a
   * kernel: the evaluation fails at this point
   * \param[in] os: output stream
   * \param[in] mpd: material property description
   */
  static void writeInvalidResultCheck(std::ostream& os,
                                      const MaterialPropertyDescription& mpd) {
    const auto v =
        useQuantities(mpd) ? mpd.output.name + ".getValue()" : mpd.output.name;
    os << "if(!tfel::math::ieee754::isfinite(" << v << ")){\n"
       << "mfront_status.report(" << invalid_result_check << "ull, true);\n"
       << "}\n";
  }  // end of writeInvalidResultCheck

  /*!
   * \brief write the code analysing the statistics `mfront_statistics` in
   * the function reporting the result: the first check which failed is
   * reported, with the number of points concerned
   * \param[in] os: output stream
   * \param[in] i: interface
   * \param[in] mpd: material property description
   */
  static void writeStatisticsAnalysis(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd) {
    const auto types = i.getTypesDescription();
    const auto prefix = i.getOutOfBoundsPolicyEnumerationPrefix();
    os << "const auto mfront_number_of_failures = static_cast<"
       << types.integer_type << ">(mfront_statistics->number_of_failures);\n"
       << "const auto mfront_number_of_unreliable_results = static_cast<"
       << types.integer_type
       << ">(mfront_statistics->number_of_unreliable_results);\n"
       << "mfront_output_status->number_of_failures = "
       << "mfront_number_of_failures;\n"
       << "mfront_output_status->number_of_unreliable_results = "
       << "mfront_number_of_unreliable_results;\n"
       << "if (mfront_number_of_failures != 0) {\n"
       << "mfront_output_status->first_failure_index = static_cast<"
       << types.integer_type << ">(mfront_statistics->first_failure_index);\n"
       << "}\n"
       // description of the points concerned by the reported check
       << "const auto mfront_points = (mfront_number_of_failures != 0) ? "
       << "\" (the evaluation failed at \" + "
       << "std::to_string(mfront_number_of_failures) + "
       << "\" points, the first one being the point \" + "
       << "std::to_string(mfront_output_status->first_failure_index) + "
       << "\")\" : "
       << "\" (\" + std::to_string(mfront_number_of_unreliable_results) + "
       << "\" points)\";\n"
       << "const auto mfront_check = mfront_statistics->first_violation;\n"
       << "if (mfront_check == " << invalid_result_check << "ull) {\n"
       << "mfront_output_status->status = -4;\n"
       << "mfront_report(\"invalid result\" + mfront_points);\n"
       << "return;\n"
       << "}\n";
    auto write_check = [&os](const unsigned long long c,
                             const std::string& msg,
                             const std::string& bounds_status,
                             const std::string& status) {
      os << "if (mfront_check == " << c << "ull) {\n"
         << "mfront_output_status->status = " << status << ";\n"
         << "mfront_output_status->bounds_status = " << bounds_status << ";\n"
         << "mfront_report(\"" << msg << "\" + mfront_points);\n"
         << "return;\n"
         << "}\n";
    };
    auto for_each_variable = [&mpd](const auto& f) {
      for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
        f(mpd.inputs[idx], idx);
      }
      f(mpd.output, mpd.inputs.size());
    };
    auto has_upper_bound = [](const VariableBoundsDescription& b) {
      return (b.boundsType == VariableBoundsDescription::UPPER) ||
             (b.boundsType == VariableBoundsDescription::LOWERANDUPPER);
    };
    auto has_lower_bound = [](const VariableBoundsDescription& b) {
      return (b.boundsType == VariableBoundsDescription::LOWER) ||
             (b.boundsType == VariableBoundsDescription::LOWERANDUPPER);
    };
    for_each_variable([&write_check, &has_upper_bound, &has_lower_bound](
                          const VariableDescription& v,
                          const std::size_t idx) {
      if (!v.hasPhysicalBounds()) {
        return;
      }
      const auto& b = v.getPhysicalBounds();
      const auto rank = "-" + std::to_string(idx + 1);
      if (has_upper_bound(b)) {
        write_check(getPhysicalBoundCheck(idx, true),
                    v.name + " is over its physical upper bound", rank, "-1");
      }
      if (has_lower_bound(b)) {
        write_check(getPhysicalBoundCheck(idx, false),
                    v.name + " is below its physical lower bound", rank,
                    "-1");
      }
    });
    // as with the `generic` interface, the status and the bounds status are
    // negative if the evaluation failed (strict policy), positive otherwise
    // (warning policy)
    const auto status = "((mfront_out_of_bounds_policy == " + prefix +
                        "_STRICT_POLICY) ? -1 : 1)";
    for_each_variable([&write_check, &has_upper_bound, &has_lower_bound, &mpd,
                       &status](const VariableDescription& v,
                                const std::size_t idx) {
      if (!v.hasBounds()) {
        return;
      }
      const auto& b = v.getBounds();
      const auto rank = status + " * " + std::to_string(idx + 1);
      if (has_upper_bound(b)) {
        write_check(getBoundCheck(mpd, idx, true),
                    v.name + " is over its upper bound", rank, status);
      }
      if (has_lower_bound(b)) {
        write_check(getBoundCheck(mpd, idx, false),
                    v.name + " is below its lower bound", rank, status);
      }
    });
  }  // end of writeStatisticsAnalysis

  BackendBase::BackendBase() = default;

  bool BackendBase::handlesDataTransfer() const {
    return false;
  }  // end of handlesDataTransfer

  bool BackendBase::mergesPointStatusesInKernels() const {
    return false;
  }  // end of mergesPointStatusesInKernels

  std::string BackendBase::getAtomicIncrement(const std::string_view) const {
    tfel::raise(
        "BackendBase::getAtomicIncrement: "
        "the statuses of the points are merged by a reduction");
  }  // end of getAtomicIncrement

  std::string BackendBase::getAtomicMinimum(const std::string_view,
                                            const std::string_view) const {
    tfel::raise(
        "BackendBase::getAtomicMinimum: "
        "the statuses of the points are merged by a reduction");
  }  // end of getAtomicMinimum

  void BackendBase::writeStatisticsDeclaration(
      std::ostream& os,
      const std::string_view device_type,
      const std::string_view device,
      const std::string_view allocation) {
    const auto pool = "mfront::generic_parallel::StatusBlockPool<"
                      "mfront_gpmp_Statistics, " +
                      std::string{device_type} + ">";
    // statistics given by the calling code (asynchronous evaluation), or
    // block of the pool, initialized here
    os << "static auto mfront_statistics_blocks = " << pool << "{};\n"
       << "auto mfront_statistics_block = std::optional<" << pool
       << "::Block>{};\n"
       << "auto* mfront_statistics = mfront_options.statistics;\n"
       << "if (mfront_statistics == nullptr) {\n"
       << "mfront_statistics_block.emplace("
       << "mfront_statistics_blocks.acquire(" << device << ", [&]{\n"
       << allocation << "\n"
       << "}));\n"
       << "mfront_statistics = mfront_statistics_block->get();\n"
       << "if (mfront_statistics == nullptr) {\n"
       << "mfront_output_status->status = -6;\n"
       << "mfront_report(\"unable to allocate the statistics of the "
       << "evaluation\");\n"
       << "return;\n"
       << "}\n"
       << "mfront_gpmp_initializeStatistics(mfront_statistics);\n"
       << "}\n";
  }  // end of writeStatisticsDeclaration

  void BackendBase::writeLocalStatisticsDeclaration(std::ostream& os) {
    os << "auto mfront_local_statistics = mfront_gpmp_Statistics{};\n"
       << "mfront_gpmp_initializeStatistics(&mfront_local_statistics);\n"
       << "auto* const mfront_statistics = "
       << "(mfront_options.statistics != nullptr) ? "
       << "mfront_options.statistics : &mfront_local_statistics;\n";
  }  // end of writeLocalStatisticsDeclaration


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

  void BackendBase::writeReportFunction(
      std::ostream& os,
      const GenericParallelMaterialPropertyInterface& i,
      const MaterialPropertyDescription& mpd) const {
    const auto name = i.getFunctionName(mpd);
    const auto types = i.getTypesDescription();
    os << "MFRONT_SHAREDOBJ void " << name << "_generic_parallel_report("
       << types.output_status_type << "* const mfront_output_status,\n"
       << "[[maybe_unused]] const mfront_gpmp_Statistics* const "
       << "mfront_statistics,\n"
       << "[[maybe_unused]] const " << types.out_of_bounds_policy_type
       << " mfront_out_of_bounds_policy){\n";
    writeErrorReporting(os);
    os << "mfront_output_status->status = 0;\n"
       << "mfront_output_status->bounds_status = 0;\n"
       << "mfront_output_status->number_of_failures = 0;\n"
       << "mfront_output_status->first_failure_index = "
       << "std::numeric_limits<" << types.integer_type << ">::max();\n"
       << "mfront_output_status->number_of_unreliable_results = 0;\n";
    if (!areRuntimeChecksDisabled(mpd)) {
      writeStatisticsAnalysis(os, i, mpd);
    }
    os << "} // end of " << name << "_generic_parallel_report\n\n";
  }  // end of writeReportFunction


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
         << "execution (queue, statistics of an asynchronous evaluation, "
         << "stack size), or null for the default options\n";
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
       << " * \\note the array giving the output (resp. an argument) holds at\n"
       << " *   least `(mfront_npoints - 1) * stride + 1` values, `stride`\n"
       << " *   being the output stride (resp. the stride of the argument)\n"
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
       << " * \\brief analyse the statistics of an asynchronous "
       << "evaluation of '" << mpd.className << "' material property, once "
       << "the computations are completed, and report the result as a "
       << "synchronous evaluation does\n"
       << " *\n"
       << " * \\param[out] mfront_output_status: output status\n"
       << " * \\param[in] mfront_statistics: statistics given by the "
       << "execution options, read on the host\n"
       << " * \\param[in] mfront_out_of_bounds_policy: out of bounds policy "
       << "of the evaluation\n"
       << " */\n"
       << "MFRONT_SHAREDOBJ void " << name << "_generic_parallel_report(\n"
       << types.output_status_type << "* const,\n"
       << "const mfront_gpmp_Statistics* const,\n"
       << "const " << types.out_of_bounds_policy_type << ");\n\n";
  }  // end of writeDeclarations

  void BackendBase::writeProxyObjects(
      std::ostream&,
      const GenericParallelMaterialPropertyInterface&,
      const MaterialPropertyDescription&) const {}  // end of writeProxyObject



  bool BackendBase::isBoundsCheckingRequired(
      const MaterialPropertyDescription& mpd) {
    // the statistics also report the invalid results
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
    // the policy is passed to the kernels with the statistics, even if the
    // material property has no bounds
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
       << "mfront_options.statistics != nullptr;\n";
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
    // errno is restored on every return. It is not used to detect the
    // errors: it is thread-local, does not exist on GPUs, and may be modified
    // by the runtimes (for instance at the initialisation of TBB). Invalid
    // results are detected by the kernels.
    os << "struct mfront_ErrnoGuard {\n"
       << "const int value = errno;\n"
       << "~mfront_ErrnoGuard() { errno = this->value; }\n"
       << "};\n"
       << "[[maybe_unused]] const auto mfront_errno_guard = "
       << "mfront_ErrnoGuard{};\n";
    os << "mfront_output_status->status = 0;\n"
       << "mfront_output_status->bounds_status = 0;\n"
       << "mfront_output_status->number_of_failures = 0;\n"
       << "mfront_output_status->first_failure_index = "
       << "std::numeric_limits<" << types.integer_type << ">::max();\n"
       << "mfront_output_status->number_of_unreliable_results = 0;\n";
    // an empty set of points is accepted without any check, as an empty
    // batch by the behaviours: nothing is read, written or submitted
    os << "if(mfront_npoints == 0){\n"
       << "return;\n"
       << "}\n";
    // writes an error and returns
    auto write_error = [&os](const int status, const std::string& message) {
      os << "mfront_output_status->status = " << status << ";\n"
         << "mfront_report(" << message << ");\n"
         << "return;\n"
         << "}\n";
    };
    if (!areRuntimeChecksDisabled(mpd)) {
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
    const auto checks = isBoundsCheckingRequired(mpd);
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
    // with the runtime checks, the evaluation is written in a lambda
    // function returning the status of the point, so that a failed check
    // ends the evaluation at this point. The variables are captured by copy:
    // nvc++ 25.11 (-mp=gpu) does not write through a pointer captured by
    // reference by such a lambda in an OpenMP target region
    if (checks) {
      os << "const auto mfront_point_status = [=]() -> "
         << "::mfront::gpmp::PointStatus {\n"
         << "auto mfront_status = ::mfront::gpmp::PointStatus{};\n";
    }
    // output
    os << "auto " << mpd.output.name << " = " << mpd.output.type << "{};\n";
    // bounds of the arguments: physical bounds first
    if ((checks) &&
        ((hasPhysicalBounds(mpd.inputs)) || (hasBounds(mpd.inputs)))) {
      os << "#ifndef NO_" << iucname << "_BOUNDS_CHECK\n";
      for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
        writePhysicalBounds(os, mpd.inputs[idx], idx);
      }
      for (std::size_t idx = 0; idx != mpd.inputs.size(); ++idx) {
        writeBounds(os, prefix, mpd, mpd.inputs[idx], idx);
      }
      os << "#endif /* NO_" << iucname << "_BOUNDS_CHECK */\n";
    }
    os << mpd.f.body;
    // bounds of the output and invalid results
    if (checks) {
      if ((mpd.output.hasPhysicalBounds()) || (mpd.output.hasBounds())) {
        os << "#ifndef NO_" << iucname << "_BOUNDS_CHECK\n";
        writePhysicalBounds(os, mpd.output, mpd.inputs.size());
        writeBounds(os, prefix, mpd, mpd.output, mpd.inputs.size());
        os << "#endif /* NO_" << iucname << "_BOUNDS_CHECK */\n";
      }
      writeInvalidResultCheck(os, mpd);
    }
    os << (handleStrides ? "mfront_output[mfront_output_stride * mfront_idx] = "
                         : "mfront_output[mfront_idx] = ")
       << mpd.output.name << (qt ? ".getValue();\n" : ";\n");
    if (!checks) {
      return;
    }
    os << "return mfront_status;\n"
       << "}();\n";
    if (!this->mergesPointStatusesInKernels()) {
      // merged by a reduction
      os << "return mfront_point_status;\n";
      return;
    }
    // the statistics are only accessed if an event occurred
    const auto increment = [this](const std::string_view c) {
      return this->getAtomicIncrement("mfront_statistics->" + std::string{c});
    };
    const auto minimum = [this](const std::string_view c,
                                const std::string_view v) {
      return this->getAtomicMinimum("mfront_statistics->" + std::string{c},
                                    v);
    };
    os << "if (mfront_point_status.hasEvent()) {\n"
       << "if (mfront_point_status.failure) {\n"
       << increment("number_of_failures") << ";\n"
       << minimum("first_failure_index",
                  "static_cast<unsigned long long>(mfront_idx)")
       << ";\n"
       << "} else {\n"
       << increment("number_of_unreliable_results") << ";\n"
       << "}\n"
       << minimum("first_violation", "mfront_point_status.violation") << ";\n"
       << "}\n";
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
    os << "if (!mfront_asynchronous) {\n"
       << i.getFunctionName(mpd) << "_generic_parallel_report("
       << "mfront_output_status, mfront_statistics, "
       << "mfront_out_of_bounds_policy);\n"
       << "}\n";
  }  // end of writeResultsAnalysis

  BackendBase::~BackendBase() noexcept = default;

}  // end of namespace mfront::generic_parallel::material_property
