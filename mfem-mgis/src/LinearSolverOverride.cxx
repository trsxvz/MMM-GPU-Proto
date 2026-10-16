/*!
 * \file   src/LinearSolverOverride.cxx
 * \brief  This file implements the reading of the linear solver which
 * replaces the linear solvers of the non linear evolution problems (see
 * `LinearSolverOverride.hxx`)
 * \date   16/10/2026
 */

#include <map>
#include <atomic>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#ifdef MFEM_MGIS_HAVE_TFEL_UTILITIES
#include "TFEL/Utilities/CxxTokenizer.hxx"
#include "TFEL/Utilities/Data.hxx"
#endif /* MFEM_MGIS_HAVE_TFEL_UTILITIES */
#include "mfem/general/globals.hpp"
#include "LinearSolverOverride.hxx"

namespace mfem_mgis::internals {

  //! \return the storage of the linear solver override
  static std::optional<LinearSolverOverride>& getLinearSolverOverrideStorage() {
    static auto o = std::optional<LinearSolverOverride>{};
    return o;
  }  // end of getLinearSolverOverrideStorage

#ifdef MFEM_MGIS_HAVE_TFEL_UTILITIES

  /*!
   * \return the parameter converted from a value read by TFEL
   * \param[in] d: value
   * \param[in] n: name of the value, for the error messages
   */
  static Parameter convertToParameter(const tfel::utilities::Data& d,
                                      const std::string& n) {
    using namespace tfel::utilities;
    // empty braces are read as an empty value
    if (d.empty()) {
      return Parameter{Parameters{}};
    }
    if (d.is<bool>()) {
      return Parameter{d.get<bool>()};
    }
    if (d.is<int>()) {
      return Parameter{d.get<int>()};
    }
    if (d.is<double>()) {
      return Parameter{d.get<double>()};
    }
    if (d.is<std::string>()) {
      return Parameter{d.get<std::string>()};
    }
    if (d.is<std::vector<Data>>()) {
      auto values = std::vector<Parameter>{};
      for (const auto& v : d.get<std::vector<Data>>()) {
        values.push_back(convertToParameter(v, n));
      }
      return Parameter{std::move(values)};
    }
    if (d.is<DataMap>()) {
      auto parameters = Parameters{};
      for (const auto& [k, v] : d.get<DataMap>()) {
        parameters.insert(may_throw, k, convertToParameter(v, k));
      }
      return Parameter{std::move(parameters)};
    }
    throw std::runtime_error(
        "unsupported value of the entry '" + n +
        "' (a map whose keys are numbers, or a structure)");
  }  // end of convertToParameter

#endif /* MFEM_MGIS_HAVE_TFEL_UTILITIES */

  LinearSolverOverride parseLinearSolverOverride(const std::string_view d) {
    const auto error = [&d](const std::string& msg) {
      return std::runtime_error("invalid linear solver override '" +
                                std::string{d} + "': " + msg);
    };
#ifdef MFEM_MGIS_HAVE_TFEL_UTILITIES
    using namespace tfel::utilities;
    auto value = Parameter{};
    try {
      auto tokenizer = CxxTokenizer{};
      tokenizer.parseString(std::string{d});
      auto p = tokenizer.begin();
      const auto data = Data::read(p, tokenizer.end());
      if (p != tokenizer.end()) {
        throw std::runtime_error("unexpected text after the description");
      }
      value = convertToParameter(data, "linear solver override");
    } catch (std::exception& e) {
      throw error(e.what());
    }
    // a name, or the entries `Name` and `Options`, as the arguments of the
    // other factories of mfem-mgis
    auto ctx = Context{};
    const auto o =
        getFactoryArgument(ctx, Parameters{{"LinearSolverOverride", value}},
                           "LinearSolverOverride");
    if (isInvalid(o)) {
      throw error(ctx.getErrorMessage());
    }
    return LinearSolverOverride{.description = std::string{d},
                                .name = o->first,
                                .parameters = o->second};
#else  /* MFEM_MGIS_HAVE_TFEL_UTILITIES */
    throw error(
        "mfem-mgis is built without the utilities of TFEL, which read the "
        "description");
#endif /* MFEM_MGIS_HAVE_TFEL_UTILITIES */
  }  // end of parseLinearSolverOverride

  void readLinearSolverOverride(const int argc,
                                const MainFunctionArguments argv) {
    auto& o = getLinearSolverOverrideStorage();
    o.reset();
    const char* d = nullptr;
    for (int i = 1; i < argc; ++i) {
      if ((std::strcmp(argv[i], linear_solver_override_option) != 0) &&
          (std::strcmp(argv[i], linear_solver_override_short_option) != 0)) {
        continue;
      }
      if (d != nullptr) {
        throw std::runtime_error(std::string{"the option "} +
                                 linear_solver_override_option +
                                 " is given twice");
      }
      // a description never starts with '-', unlike the next option
      if ((i + 1 == argc) || (argv[i + 1][0] == '-')) {
        throw std::runtime_error(std::string{"no value given to the option "} +
                                 linear_solver_override_option);
      }
      d = argv[i + 1];
    }
    if (d == nullptr) {
      d = std::getenv(linear_solver_override_variable);
    }
    if ((d == nullptr) || (std::strlen(d) == 0)) {
      return;
    }
    o = parseLinearSolverOverride(d);
  }  // end of readLinearSolverOverride

#ifdef MFEM_USE_MPI

  void checkLinearSolverOverride(const MPI_Comm comm) {
    const auto& o = getLinearSolverOverride();
    const auto description = o.has_value() ? o->description : std::string{};
    // the description of the first process
    auto size = static_cast<int>(description.size());
    MPI_Bcast(&size, 1, MPI_INT, 0, comm);
    auto first = description;
    first.resize(static_cast<std::size_t>(size));
    MPI_Bcast(first.data(), size, MPI_CHAR, 0, comm);
    auto differ = static_cast<int>(first != description);
    MPI_Allreduce(MPI_IN_PLACE, &differ, 1, MPI_INT, MPI_LOR, comm);
    if (differ != 0) {
      throw std::runtime_error(
          "the processes have different linear solver overrides (the "
          "environment variable MFEM_MGIS_LINEAR_SOLVER_OVERRIDE shall be "
          "given to all the processes)");
    }
  }  // end of checkLinearSolverOverride

#endif /* MFEM_USE_MPI */

  const std::optional<LinearSolverOverride>&
  getLinearSolverOverride() noexcept {
    return getLinearSolverOverrideStorage();
  }  // end of getLinearSolverOverride

  void reportLinearSolverOverride(const std::string_view replaced) noexcept {
    static auto reported = std::atomic<bool>{false};
    const auto& o = getLinearSolverOverride();
    if ((!o.has_value()) || reported.exchange(true)) {
      return;
    }
    try {
      mfem::out << "mfem-mgis: the linear solver override '" << o->name
                << "' replaces " << replaced << std::endl;
    } catch (...) {
    }
  }  // end of reportLinearSolverOverride

  std::pair<std::string_view, const Parameters&> selectLinearSolver(
      const std::string_view n, const Parameters& p) noexcept {
    const auto& o = getLinearSolverOverride();
    if (!o.has_value()) {
      return {n, p};
    }
    try {
      reportLinearSolverOverride("the linear solver '" + std::string{n} + "'");
    } catch (...) {
    }
    return {o->name, o->parameters};
  }  // end of selectLinearSolver

}  // end of namespace mfem_mgis::internals
