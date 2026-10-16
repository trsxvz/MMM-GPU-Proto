/*!
 * \file   src/LinearSolverOverride.hxx
 * \brief  This file declares the reading of the linear solver and of the
 * preconditioner which replace the ones of the non linear evolution
 * problems, given by the options `--linear-solver-override` and
 * `--preconditioner-override` or by the environment variables
 * `MFEM_MGIS_LINEAR_SOLVER_OVERRIDE` and `MFEM_MGIS_PRECONDITIONER_OVERRIDE`
 * (internal header)
 * \date   16/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_LINEARSOLVEROVERRIDE_HXX
#define LIB_MFEM_MGIS_SRC_LINEARSOLVEROVERRIDE_HXX

#include <string>
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
   * variable `MFEM_MGIS_LINEAR_SOLVER_OVERRIDE`, and the preconditioner
   * override, given by the option `--preconditioner-override` (or `-pco`),
   * or else by the environment variable `MFEM_MGIS_PRECONDITIONER_OVERRIDE`
   * (see `parseFactoryArgumentOverride`), returned afterwards by
   * `getLinearSolverOverride` and `getPreconditionerOverride`. An empty
   * value cancels an override.
   * \param[in] argc: number of arguments
   * \param[in] argv: arguments
   * \throws std::runtime_error on an invalid description, or if an option
   * is given twice or without value
   */
  void readLinearSolverOverrides(const int argc,
                                 const MainFunctionArguments argv);
#ifdef MFEM_USE_MPI
  /*!
   * \brief checks that all the processes of a communicator have the same
   * overrides (the environment variables are not always given to the
   * processes of the other nodes), collective
   * \param[in] comm: communicator
   * \throws std::runtime_error on all the processes if they differ
   */
  void checkLinearSolverOverrides(const MPI_Comm comm);
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
   * \brief reports, once, on the output stream of MFEM, that the linear
   * solver override replaces a linear solver of the program
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

#endif /* LIB_MFEM_MGIS_SRC_LINEARSOLVEROVERRIDE_HXX */
