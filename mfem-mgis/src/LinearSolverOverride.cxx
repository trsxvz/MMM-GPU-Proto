/*!
 * \file   src/LinearSolverOverride.cxx
 * \brief  This file implements the reading of the linear solver which
 * replaces the linear solvers of the non linear evolution problems (see
 * `LinearSolverOverride.hxx`)
 * \date   16/10/2026
 */

#include <map>
#include <algorithm>
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
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/PreconditionerFactory.hxx"
#include "LinearSolverModules.hxx"
#include "LinearSolverOverride.hxx"

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
   * \brief checks that the linear solver override and the preconditioner
   * override are registered in the factories, so that a misspelt name is
   * refused at the initialization, and that `None` (no preconditioner) has
   * no options
   */
  static void checkOverridesNames() {
    const auto check = [](const std::optional<FactoryArgumentOverride>& o,
                          const std::vector<std::string>& sequential,
                          const std::vector<std::string>& parallel,
                          const std::string& what) {
      if ((!o.has_value()) || (isRegistered(sequential, parallel, o->name))) {
        return;
      }
      auto msg = "invalid " + what + " '" + o->description +
                 "': unknown name '" + o->name + "'";
      const auto hint = getUnavailableLibraryHint(o->name);
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
          "linear solver override");
    auto sequential_preconditioners =
        PreconditionerFactory<false>::getFactory().getNames();
    sequential_preconditioners.push_back("None");
    check(getPreconditionerOverride(), sequential_preconditioners,
          parallel_preconditioners, "preconditioner override");
    const auto& pc = getPreconditionerOverride();
    if (pc.has_value() && (pc->name == "None") && (!pc->options.empty())) {
      throw std::runtime_error("invalid preconditioner override '" +
                               pc->description +
                               "': no options expected for 'None'");
    }
  }  // end of checkOverridesNames

  void readLinearSolverOverrides(const int argc,
                                 const MainFunctionArguments argv) {
    readOverride(getLinearSolverOverrideStorage(), argc, argv,
                 linear_solver_override_option,
                 linear_solver_override_short_option,
                 linear_solver_override_variable, "linear solver override");
    readOverride(getPreconditionerOverrideStorage(), argc, argv,
                 preconditioner_override_option,
                 preconditioner_override_short_option,
                 preconditioner_override_variable, "preconditioner override");
    checkOverridesNames();
  }  // end of readLinearSolverOverrides

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

  void checkLinearSolverOverrides(const MPI_Comm comm) {
    const auto description =
        [](const std::optional<FactoryArgumentOverride>& o) {
          return o.has_value() ? o->description : std::string{};
        };
    const auto ls =
        differsFromFirstProcess(comm, description(getLinearSolverOverride()));
    const auto pc =
        differsFromFirstProcess(comm, description(getPreconditionerOverride()));
    if (ls || pc) {
      throw std::runtime_error(
          "the processes have different linear solver or preconditioner "
          "overrides (the environment variables "
          "MFEM_MGIS_LINEAR_SOLVER_OVERRIDE and "
          "MFEM_MGIS_PRECONDITIONER_OVERRIDE shall be given to all the "
          "processes)");
    }
  }  // end of checkLinearSolverOverrides

#endif /* MFEM_USE_MPI */

  const std::optional<FactoryArgumentOverride>&
  getLinearSolverOverride() noexcept {
    return getLinearSolverOverrideStorage();
  }  // end of getLinearSolverOverride

  const std::optional<FactoryArgumentOverride>&
  getPreconditionerOverride() noexcept {
    return getPreconditionerOverrideStorage();
  }  // end of getPreconditionerOverride

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

  /*!
   * \brief reports, once, on the output stream of MFEM, that the
   * preconditioner override is given to a linear solver
   * \param[in] n: name of the linear solver
   */
  static void reportPreconditionerOverride(const std::string_view n) noexcept {
    static auto reported = std::atomic<bool>{false};
    const auto& o = getPreconditionerOverride();
    if ((!o.has_value()) || reported.exchange(true)) {
      return;
    }
    try {
      mfem::out << "mfem-mgis: the preconditioner override '" << o->name
                << "' is given to the linear solver '" << n << "'" << std::endl;
    } catch (...) {
    }
  }  // end of reportPreconditionerOverride

  std::pair<std::string, Parameters> selectLinearSolver(
      const std::string_view n, const Parameters& p) {
    const auto& ls = getLinearSolverOverride();
    auto r = std::pair<std::string, Parameters>{n, p};
    if (ls.has_value()) {
      reportLinearSolverOverride("the linear solver '" + std::string{n} + "'");
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
