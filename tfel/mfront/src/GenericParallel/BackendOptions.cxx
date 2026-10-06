/*!
 * \file   mfront/src/GenericParallel/BackendOptions.cxx
 * \brief  This file implements the options shared by the backends of the
 * `generic-parallel` interfaces for material properties and behaviours.
 * \date   06/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include "TFEL/Raise.hxx"
#include "TFEL/Utilities/CxxTokenizer.hxx"
#include "MFront/GenericParallel/BackendOptions.hxx"

namespace mfront::generic_parallel {

  void declareExecutionPolicyOption(tfel::utilities::DataMapValidator& v) {
    v.addDataTypeValidator<std::string>("execution_policy");
  }  // end of declareExecutionPolicyOption

  std::string getExecutionPolicy(const tfel::utilities::DataMap& opts) {
    if (!contains(opts, "execution_policy")) {
      return "std::execution::par_unseq";
    }
    const auto p = get<std::string>(opts, "execution_policy");
    if ((p == "seq") || (p == "sequenced_policy")) {
      return "std::execution::seq";
    } else if ((p == "unseq") || (p == "unsequenced_policy")) {
      return "std::execution::unseq";
    } else if ((p == "par") || (p == "parallel_policy")) {
      return "std::execution::par";
    } else if ((p != "par_unseq") && (p != "parallel_unsequenced_policy")) {
      tfel::raise(
          "getExecutionPolicy: unsupported execution policy '" + p +
          "'. Valid execution policies are: 'sequenced_policy' (or "
          "'seq'), 'unsequenced_policy' (or 'unseq'), 'parallel_policy' (or "
          "'par') and 'parallel_unsequenced_policy' (or 'par_unseq')");
    }
    return "std::execution::par_unseq";
  }  // end of getExecutionPolicy

  void declareDeviceStackSizeOption(tfel::utilities::DataMapValidator& v) {
    v.addStrictlyPositiveIntegerCheck("device_stack_size");
  }  // end of declareDeviceStackSizeOption

  int getDeviceStackSize(const tfel::utilities::DataMap& opts) {
    if (!contains(opts, "device_stack_size")) {
      return 0;
    }
    return get<int>(opts, "device_stack_size");
  }  // end of getDeviceStackSize

  void declareKokkosExecutionSpaceOption(tfel::utilities::DataMapValidator& v) {
    v.addDataTypeValidator<std::string>("execution_space");
  }  // end of declareKokkosExecutionSpaceOption

  std::string getKokkosExecutionSpace(const tfel::utilities::DataMap& opts) {
    if (!contains(opts, "execution_space")) {
      return "Kokkos::DefaultExecutionSpace";
    }
    // the execution space is written in the generated code
    const auto s = get<std::string>(opts, "execution_space");
    auto valid = !s.empty();
    for (auto p = std::string::size_type{}; valid;) {
      const auto e = s.find("::", p);
      const auto n = s.substr(p, e == std::string::npos ? e : e - p);
      valid = tfel::utilities::CxxTokenizer::isValidIdentifier(n, false);
      if (e == std::string::npos) {
        break;
      }
      p = e + 2;
    }
    if (!valid) {
      tfel::raise("getKokkosExecutionSpace: invalid execution space '" + s +
                  "'");
    }
    return s;
  }  // end of getKokkosExecutionSpace

  void declareDataTransferPolicyOption(tfel::utilities::DataMapValidator& v) {
    v.addDataTypeValidator<std::string>("data_transfer_policy");
  }  // end of declareDataTransferPolicyOption

  bool areDataTransfersAutomatic(const tfel::utilities::DataMap& opts) {
    if (!contains(opts, "data_transfer_policy")) {
      return true;
    }
    const auto p = get<std::string>(opts, "data_transfer_policy");
    if ((p != "auto") && (p != "automatic") && (p != "none")) {
      tfel::raise("areDataTransfersAutomatic: invalid data transfer policy '" +
                  p +
                  "'. Valid policies are 'automatic' (or 'auto') and 'none'");
    }
    return p != "none";
  }  // end of areDataTransfersAutomatic

}  // end of namespace mfront::generic_parallel
