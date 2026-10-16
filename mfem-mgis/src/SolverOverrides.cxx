/*!
 * \file   src/SolverOverrides.cxx
 * \brief  This file implements the reading of the linear solver, of the
 * preconditioner and of the non linear solver which replace the ones of the
 * non linear evolution problems (see `SolverOverrides.hxx`)
 * \date   16/10/2026
 */

#include <map>
#include <set>
#include <mutex>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#ifdef MFEM_MGIS_HAVE_TFEL_UTILITIES
#include "TFEL/Utilities/CxxTokenizer.hxx"
#include "TFEL/Utilities/Data.hxx"
#endif /* MFEM_MGIS_HAVE_TFEL_UTILITIES */
#include "mfem/general/globals.hpp"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/PreconditionerFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverFactory.hxx"
#include "LinearSolverModules.hxx"
#include "SolverOverrides.hxx"

namespace mfem_mgis::internals {

  //! \return the storage of the linear solver override
  static std::optional<FactoryArgumentOverride>&
  getLinearSolverOverrideStorage() {
    static auto o = std::optional<FactoryArgumentOverride>{};
    return o;
  }  // end of getLinearSolverOverrideStorage

  //! \return the storage of the preconditioner override
  static std::optional<FactoryArgumentOverride>&
  getPreconditionerOverrideStorage() {
    static auto o = std::optional<FactoryArgumentOverride>{};
    return o;
  }  // end of getPreconditionerOverrideStorage

  //! \return the storage of the non linear solver override
  static std::optional<FactoryArgumentOverride>&
  getNonLinearSolverOverrideStorage() {
    static auto o = std::optional<FactoryArgumentOverride>{};
    return o;
  }  // end of getNonLinearSolverOverrideStorage

  /*!
   * \brief prints a message on the output stream of MFEM, unless it was
   * already printed: each replacement is reported once, whatever the number
   * of problems
   * \param[in] message: message, without the `mfem-mgis: ` prefix
   */
  static void reportOnce(const std::string& message) noexcept {
    static auto m = std::mutex{};
    static auto reported = std::set<std::string>{};
    try {
      const auto lock = std::lock_guard<std::mutex>{m};
      if (reported.insert(message).second) {
        mfem::out << "mfem-mgis: " << message << std::endl;
      }
    } catch (...) {
    }
  }  // end of reportOnce

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

  FactoryArgumentOverride parseFactoryArgumentOverride(
      const std::string_view d, const std::string_view what) {
    const auto error = [&d, &what](const std::string& msg) {
      return std::runtime_error("invalid " + std::string{what} + " '" +
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
      value = convertToParameter(data, std::string{what});
    } catch (std::exception& e) {
      throw error(e.what());
    }
    // a name, or the entries `Name` and `Options`, as the arguments of the
    // other factories of mfem-mgis
    auto ctx = Context{};
    const auto o =
        getFactoryArgument(ctx, Parameters{{"Override", value}}, "Override");
    if (isInvalid(o)) {
      throw error(ctx.getErrorMessage());
    }
    return FactoryArgumentOverride{
        .description = std::string{d}, .name = o->first, .options = o->second};
#else  /* MFEM_MGIS_HAVE_TFEL_UTILITIES */
    throw error(
        "mfem-mgis is built without the utilities of TFEL, which read the "
        "description");
#endif /* MFEM_MGIS_HAVE_TFEL_UTILITIES */
  }  // end of parseFactoryArgumentOverride

  /*!
   * \brief reads an override given by an option, or else by an environment
   * variable
   * \param[out] o: override
   * \param[in] argc: number of arguments
   * \param[in] argv: arguments
   * \param[in] option: long name of the option
   * \param[in] short_option: short name of the option
   * \param[in] variable: environment variable
   * \param[in] what: kind of override, for the error messages
   */
  static void readOverride(std::optional<FactoryArgumentOverride>& o,
                           const int argc,
                           const MainFunctionArguments argv,
                           const char* const option,
                           const char* const short_option,
                           const char* const variable,
                           const std::string_view what) {
    o.reset();
    const char* d = nullptr;
    for (int i = 1; i < argc; ++i) {
      if ((std::strcmp(argv[i], option) != 0) &&
          (std::strcmp(argv[i], short_option) != 0)) {
        continue;
      }
      if (d != nullptr) {
        throw std::runtime_error(std::string{"the option "} + option +
                                 " is given twice");
      }
      // a description never starts with '-', unlike the next option
      if ((i + 1 == argc) || (argv[i + 1][0] == '-')) {
        throw std::runtime_error(std::string{"no value given to the option "} +
                                 option);
      }
      d = argv[i + 1];
    }
    if (d == nullptr) {
      d = std::getenv(variable);
    }
    if ((d == nullptr) || (std::strlen(d) == 0)) {
      return;
    }
    o = parseFactoryArgumentOverride(d, what);
  }  // end of readOverride

  /*!
   * \return if a name is registered in a factory, in sequential or in
   * parallel
   * \param[in] sequential: names of the sequential factory
   * \param[in] parallel: names of the parallel factory
   * \param[in] n: name
   */
  static bool isRegistered(const std::vector<std::string>& sequential,
                           const std::vector<std::string>& parallel,
                           const std::string& n) {
    return (std::find(sequential.begin(), sequential.end(), n) !=
            sequential.end()) ||
           (std::find(parallel.begin(), parallel.end(), n) != parallel.end());
  }  // end of isRegistered

  /*!
   * \brief checks that the solver overrides are registered in the
   * factories, so that a misspelt name is refused at the initialization, and
   * that `None` (no preconditioner) has no options
   */
  static void checkOverridesNames() {
    const auto check = [](const std::optional<FactoryArgumentOverride>& o,
                          const std::vector<std::string>& sequential,
                          const std::vector<std::string>& parallel,
                          const std::string& what,
                          std::string (*const get_hint)(std::string_view)) {
      if ((!o.has_value()) || (isRegistered(sequential, parallel, o->name))) {
        return;
      }
      auto msg = "invalid " + what + " '" + o->description +
                 "': unknown name '" + o->name + "'";
      const auto hint = get_hint(o->name);
      if (!hint.empty()) {
        msg += " (" + hint + ")";
      }
      throw std::runtime_error(msg);
    };
#ifdef MFEM_USE_MPI
    const auto parallel_solvers =
        LinearSolverFactory<true>::getFactory().getNames();
    const auto parallel_preconditioners =
        PreconditionerFactory<true>::getFactory().getNames();
#else  /* MFEM_USE_MPI */
    const auto parallel_solvers = std::vector<std::string>{};
    const auto parallel_preconditioners = std::vector<std::string>{};
#endif /* MFEM_USE_MPI */
    check(getLinearSolverOverride(),
          LinearSolverFactory<false>::getFactory().getNames(), parallel_solvers,
          "linear solver override", getUnavailableLibraryHint);
    auto sequential_preconditioners =
        PreconditionerFactory<false>::getFactory().getNames();
    sequential_preconditioners.push_back("None");
    check(getPreconditionerOverride(), sequential_preconditioners,
          parallel_preconditioners, "preconditioner override",
          getUnavailableLibraryHint);
    check(getNonLinearSolverOverride(),
          NonLinearSolverFactory::get().getNames(), std::vector<std::string>{},
          "non linear solver override",
          getUnavailableNonLinearSolverLibraryHint);
    const auto& pc = getPreconditionerOverride();
    if (pc.has_value() && (pc->name == "None") && (!pc->options.empty())) {
      throw std::runtime_error("invalid preconditioner override '" +
                               pc->description +
                               "': no options expected for 'None'");
    }
  }  // end of checkOverridesNames

  void readSolverOverrides(const int argc, const MainFunctionArguments argv) {
    readOverride(getLinearSolverOverrideStorage(), argc, argv,
                 linear_solver_override_option,
                 linear_solver_override_short_option,
                 linear_solver_override_variable, "linear solver override");
    readOverride(getPreconditionerOverrideStorage(), argc, argv,
                 preconditioner_override_option,
                 preconditioner_override_short_option,
                 preconditioner_override_variable, "preconditioner override");
    readOverride(getNonLinearSolverOverrideStorage(), argc, argv,
                 non_linear_solver_override_option,
                 non_linear_solver_override_short_option,
                 non_linear_solver_override_variable,
                 "non linear solver override");
    checkOverridesNames();
  }  // end of readSolverOverrides

#ifdef MFEM_USE_MPI

  /*!
   * \return if a description differs from the one of the first process
   * \param[in] comm: communicator
   * \param[in] description: description of the process
   */
  static bool differsFromFirstProcess(const MPI_Comm comm,
                                      const std::string& description) {
    auto size = static_cast<int>(description.size());
    MPI_Bcast(&size, 1, MPI_INT, 0, comm);
    auto first = description;
    first.resize(static_cast<std::size_t>(size));
    MPI_Bcast(first.data(), size, MPI_CHAR, 0, comm);
    auto differ = static_cast<int>(first != description);
    MPI_Allreduce(MPI_IN_PLACE, &differ, 1, MPI_INT, MPI_LOR, comm);
    return differ != 0;
  }  // end of differsFromFirstProcess

  void checkSolverOverrides(const MPI_Comm comm) {
    const auto description =
        [](const std::optional<FactoryArgumentOverride>& o) {
          return o.has_value() ? o->description : std::string{};
        };
    const auto ls =
        differsFromFirstProcess(comm, description(getLinearSolverOverride()));
    const auto pc =
        differsFromFirstProcess(comm, description(getPreconditionerOverride()));
    const auto nls = differsFromFirstProcess(
        comm, description(getNonLinearSolverOverride()));
    if (ls || pc || nls) {
      throw std::runtime_error(
          "the processes have different solver overrides (the environment "
          "variables MFEM_MGIS_LINEAR_SOLVER_OVERRIDE, "
          "MFEM_MGIS_PRECONDITIONER_OVERRIDE and "
          "MFEM_MGIS_NON_LINEAR_SOLVER_OVERRIDE shall be given to all the "
          "processes)");
    }
  }  // end of checkSolverOverrides

#endif /* MFEM_USE_MPI */

  const std::optional<FactoryArgumentOverride>&
  getLinearSolverOverride() noexcept {
    return getLinearSolverOverrideStorage();
  }  // end of getLinearSolverOverride

  const std::optional<FactoryArgumentOverride>&
  getPreconditionerOverride() noexcept {
    return getPreconditionerOverrideStorage();
  }  // end of getPreconditionerOverride

  const std::optional<FactoryArgumentOverride>&
  getNonLinearSolverOverride() noexcept {
    return getNonLinearSolverOverrideStorage();
  }  // end of getNonLinearSolverOverride

  void reportNonLinearSolverOverride(const std::string_view replaced) noexcept {
    const auto& o = getNonLinearSolverOverride();
    if (!o.has_value()) {
      return;
    }
    try {
      reportOnce("the non linear solver override '" + o->name +
                 "' replaces the non linear solver '" + std::string{replaced} +
                 "'");
    } catch (...) {
    }
  }  // end of reportNonLinearSolverOverride

  Parameters mergeNonLinearSolverOverrideOptions(const Parameters& p) {
    auto merged = p;
    const auto& o = getNonLinearSolverOverride();
    if (o.has_value()) {
      // the other options of the override were applied when it was built
      for (const auto& [n, v] : o->options) {
        if (contains(p, n)) {
          merged.replaceOrInsert(n, v);
        }
      }
    }
    return merged;
  }  // end of mergeNonLinearSolverOverrideOptions

  std::pair<Parameters, std::vector<std::string>>
  splitNonLinearSolverParameters(const Parameters& p) {
    auto common = getIterativeSolverParametersList();
    common.push_back("DiscardLinearSolverFailure");
    const auto& o = getNonLinearSolverOverride();
    auto r = std::pair<Parameters, std::vector<std::string>>{};
    for (const auto& [n, v] : p) {
      // a parameter also given by the options of the override is accepted
      // by it
      if ((std::find(common.begin(), common.end(), n) != common.end()) ||
          (o.has_value() && contains(o->options, n))) {
        r.first.insert(may_throw, n, v);
      } else {
        r.second.push_back(n);
      }
    }
    return r;
  }  // end of splitNonLinearSolverParameters

  void reportIgnoredNonLinearSolverParameters(
      const std::vector<std::string>& ignored,
      const std::string_view reason) noexcept {
    const auto& o = getNonLinearSolverOverride();
    if ((!o.has_value()) || ignored.empty()) {
      return;
    }
    try {
      auto names = std::string{};
      for (const auto& n : ignored) {
        names += (names.empty() ? "'" : ", '") + n + "'";
      }
      reportOnce("the non linear solver override '" + o->name +
                 "' refuses the parameters given by the program to its non "
                 "linear solver: the parameters " +
                 names + " are ignored. Reason of the refusal:\n" +
                 std::string{reason});
    } catch (...) {
    }
  }  // end of reportIgnoredNonLinearSolverParameters

  void reportLinearSolverOverride(const std::string_view replaced) noexcept {
    const auto& o = getLinearSolverOverride();
    if (!o.has_value()) {
      return;
    }
    try {
      reportOnce("the linear solver override '" + o->name + "' replaces " +
                 std::string{replaced});
    } catch (...) {
    }
  }  // end of reportLinearSolverOverride

  /*!
   * \brief reports, once per linear solver, on the output stream of MFEM,
   * that the preconditioner override is given to a linear solver
   * \param[in] n: name of the linear solver
   */
  static void reportPreconditionerOverride(const std::string_view n) noexcept {
    const auto& o = getPreconditionerOverride();
    if (!o.has_value()) {
      return;
    }
    try {
      reportOnce("the preconditioner override '" + o->name +
                 "' is given to the linear solver '" + std::string{n} + "'");
    } catch (...) {
    }
  }  // end of reportPreconditionerOverride

  std::pair<std::string, Parameters> selectLinearSolver(
      const std::string_view n, const Parameters& p) {
    const auto& ls = getLinearSolverOverride();
    auto r = std::pair<std::string, Parameters>{n, p};
    if (ls.has_value()) {
      // an empty name is given when the linear solver built by the program
      // is replaced, which is reported by the caller
      if (!n.empty()) {
        reportLinearSolverOverride("the linear solver '" + std::string{n} +
                                   "'");
      }
      r = {ls->name, ls->options};
    }
    const auto& pc = getPreconditionerOverride();
    if (pc.has_value()) {
      reportPreconditionerOverride(r.first);
      // the preconditioner is replaced; `None` only removes it (the
      // solvers which accept no preconditioner refuse the option)
      auto& options = r.second;
      options = remove(options, std::vector<std::string>{"Preconditioner"});
      if (pc->name != "None") {
        options.insert(
            may_throw, "Preconditioner",
            Parameters{{"Name", pc->name}, {"Options", pc->options}});
      }
    }
    return r;
  }  // end of selectLinearSolver

}  // end of namespace mfem_mgis::internals
