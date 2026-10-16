/*!
 * \file   src/SolverOverrides.hxx
 * \brief  This file declares the reading of the linear solver, of the
 * preconditioner and of the non linear solver which replace the ones of the
 * non linear evolution problems, given by the options
 * `--linear-solver-override`, `--preconditioner-override` and
 * `--non-linear-solver-override` or by the environment variables
 * `MFEM_MGIS_LINEAR_SOLVER_OVERRIDE`, `MFEM_MGIS_PRECONDITIONER_OVERRIDE` and
 * `MFEM_MGIS_NON_LINEAR_SOLVER_OVERRIDE` (internal header)
 * \date   16/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_SOLVEROVERRIDES_HXX
#define LIB_MFEM_MGIS_SRC_SOLVEROVERRIDES_HXX

#include <string>
#include <vector>
#include <utility>
#include <optional>
#include <string_view>
#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis::internals {

  //! \brief long name of the option giving the linear solver override
  inline constexpr const char* linear_solver_override_option =
      "--linear-solver-override";
  //! \brief short name of the option `linear_solver_override_option`
  inline constexpr const char* linear_solver_override_short_option = "-lso";
  //! \brief environment variable giving the linear solver override
  inline constexpr const char* linear_solver_override_variable =
      "MFEM_MGIS_LINEAR_SOLVER_OVERRIDE";
  //! \brief long name of the option giving the preconditioner override
  inline constexpr const char* preconditioner_override_option =
      "--preconditioner-override";
  //! \brief short name of the option `preconditioner_override_option`
  inline constexpr const char* preconditioner_override_short_option = "-pco";
  //! \brief environment variable giving the preconditioner override
  inline constexpr const char* preconditioner_override_variable =
      "MFEM_MGIS_PRECONDITIONER_OVERRIDE";

  //! \brief long name of the option giving the non linear solver override
  inline constexpr const char* non_linear_solver_override_option =
      "--non-linear-solver-override";
  //! \brief short name of the option `non_linear_solver_override_option`
  inline constexpr const char* non_linear_solver_override_short_option =
      "-nlso";
  //! \brief environment variable giving the non linear solver override
  inline constexpr const char* non_linear_solver_override_variable =
      "MFEM_MGIS_NON_LINEAR_SOLVER_OVERRIDE";

  //! \brief object of a factory given by the user, with its description
  struct FactoryArgumentOverride {
    //! \brief description, as given
    std::string description;
    //! \brief name of the object
    std::string name;
    //! \brief options of the object
    Parameters options;
  };

  /*!
   * \return the object of a factory described in the syntax of the data of
   * TFEL (the syntax of the options of the domain specific languages of
   * MFront) as a factory argument of mfem-mgis: a name (`CGSolver` or
   * `"CGSolver"`), or a map with the entries `Name` and `Options`
   * (optional), for instance `{Name: "BelosSolver", Options: {Solver: "CG",
   * Preconditioner: {Name: "MueLu"}}}`. The integers, reals, strings,
   * booleans, lists and maps (empty braces included) are converted to
   * parameters.
   * \param[in] d: description
   * \param[in] what: kind of object, for the error messages (`linear solver
   * override`)
   * \throws std::runtime_error if the description is invalid, or if
   * mfem-mgis is built without the utilities of TFEL
   */
  [[nodiscard]] FactoryArgumentOverride parseFactoryArgumentOverride(
      std::string_view d, std::string_view what);
  /*!
   * \brief reads the linear solver override, given by the option
   * `--linear-solver-override` (or `-lso`), or else by the environment
   * variable `MFEM_MGIS_LINEAR_SOLVER_OVERRIDE`, the preconditioner
   * override, given by the option `--preconditioner-override` (or `-pco`),
   * or else by the environment variable `MFEM_MGIS_PRECONDITIONER_OVERRIDE`,
   * and the non linear solver override, given by the option
   * `--non-linear-solver-override` (or `-nlso`), or else by the environment
   * variable `MFEM_MGIS_NON_LINEAR_SOLVER_OVERRIDE` (see
   * `parseFactoryArgumentOverride`), returned afterwards by
   * `getLinearSolverOverride`, `getPreconditionerOverride` and
   * `getNonLinearSolverOverride`. An empty value cancels an override.
   * \param[in] argc: number of arguments
   * \param[in] argv: arguments
   * \throws std::runtime_error on an invalid description, or if an option
   * is given twice or without value
   */
  void readSolverOverrides(const int argc, const MainFunctionArguments argv);
#ifdef MFEM_USE_MPI
  /*!
   * \brief checks that all the processes of a communicator have the same
   * overrides (the environment variables are not always given to the
   * processes of the other nodes), collective
   * \param[in] comm: communicator
   * \throws std::runtime_error on all the processes if they differ
   */
  void checkSolverOverrides(const MPI_Comm comm);
#endif /* MFEM_USE_MPI */
  /*!
   * \return the linear solver which replaces the linear solvers of the non
   * linear evolution problems (see `setLinearSolver`), if any
   */
  [[nodiscard]] const std::optional<FactoryArgumentOverride>&
  getLinearSolverOverride() noexcept;
  /*!
   * \return the preconditioner which replaces the preconditioners of the
   * linear solvers of the non linear evolution problems (`Preconditioner`
   * option of the linear solvers set by `setLinearSolver`), if any
   */
  [[nodiscard]] const std::optional<FactoryArgumentOverride>&
  getPreconditionerOverride() noexcept;
  /*!
   * \return the non linear solver which replaces the non linear solvers of
   * the non linear evolution problems (parameter `NonLinearSolver`, Newton
   * by default), if any
   */
  [[nodiscard]] const std::optional<FactoryArgumentOverride>&
  getNonLinearSolverOverride() noexcept;
  /*!
   * \brief reports, once per replaced non linear solver, on the output
   * stream of MFEM, that the non linear solver override replaces a non
   * linear solver of the program
   * \param[in] replaced: name of the replaced non linear solver
   */
  void reportNonLinearSolverOverride(std::string_view replaced) noexcept;
  /*!
   * \return the parameters given by the program to its non linear solver
   * (`setSolverParameters`), the options of the non linear solver override
   * taking precedence over the parameters of the same name (its other
   * options were applied when it was built)
   * \param[in] p: parameters of the program
   */
  [[nodiscard]] Parameters mergeNonLinearSolverOverrideOptions(
      const Parameters& p);
  /*!
   * \return the parameters accepted by all the non linear solvers
   * (`VerbosityLevel`, `RelativeTolerance`, `AbsoluteTolerance`,
   * `MaximumNumberOfIterations`, `DiscardLinearSolverFailure`) and the ones
   * also given by the options of the non linear solver override, then the
   * names of the other parameters
   * \param[in] p: parameters given by the program to its non linear solver
   */
  [[nodiscard]] std::pair<Parameters, std::vector<std::string>>
  splitNonLinearSolverParameters(const Parameters& p);
  /*!
   * \brief reports, once per message, on the output stream of MFEM, the
   * parameters of the program ignored because the non linear solver
   * override refuses them
   * \param[in] ignored: names of the ignored parameters (nothing is
   * reported if empty)
   * \param[in] reason: message of the refusal
   */
  void reportIgnoredNonLinearSolverParameters(
      const std::vector<std::string>& ignored,
      std::string_view reason) noexcept;
  /*!
   * \brief reports, once per replaced linear solver, on the output stream
   * of MFEM, that the linear solver override replaces a linear solver of the
   * program
   * \param[in] replaced: description of the replaced linear solver
   */
  void reportLinearSolverOverride(std::string_view replaced) noexcept;
  /*!
   * \return the name and the parameters of the linear solver to build: the
   * linear solver override, if any, or else the given ones, whose option
   * `Preconditioner` is the preconditioner override, if any. The
   * replacements are reported once on the output stream of MFEM.
   * \param[in] n: name of the linear solver requested by the program
   * \param[in] p: parameters of the linear solver requested by the program
   */
  [[nodiscard]] std::pair<std::string, Parameters> selectLinearSolver(
      std::string_view n, const Parameters& p);

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_SOLVEROVERRIDES_HXX */
