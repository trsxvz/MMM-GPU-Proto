/*!
 * \file   src/LinearSolverOverride.hxx
 * \brief  This file declares the reading of the linear solver which
 * replaces the linear solvers of the non linear evolution problems, given
 * by the option `--linear-solver-override` or the environment variable
 * `MFEM_MGIS_LINEAR_SOLVER_OVERRIDE` (internal header)
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

  //! \brief linear solver replacing the ones set by `setLinearSolver`
  struct LinearSolverOverride {
    //! \brief description, as given
    std::string description;
    //! \brief name of the linear solver
    std::string name;
    //! \brief parameters of the linear solver
    Parameters parameters;
  };

  /*!
   * \return the linear solver described in the syntax of the data of TFEL
   * (the syntax of the options of the domain specific languages of MFront)
   * as a factory argument of mfem-mgis: a name (`CGSolver` or `"CGSolver"`),
   * or a map with the entries `Name` and `Options` (optional), for instance
   * `{Name: "BelosSolver", Options: {Solver: "CG", Preconditioner: {Name:
   * "MueLu"}}}`. The integers, reals, strings, booleans, lists and maps
   * (empty braces included) are converted to parameters.
   * \param[in] d: description
   * \throws std::runtime_error if the description is invalid, or if
   * mfem-mgis is built without the utilities of TFEL
   */
  [[nodiscard]] LinearSolverOverride parseLinearSolverOverride(
      std::string_view d);
  /*!
   * \brief reads the linear solver override, given by the option
   * `--linear-solver-override` (or `-lso`), or else by the environment
   * variable `MFEM_MGIS_LINEAR_SOLVER_OVERRIDE` (see
   * `parseLinearSolverOverride`), returned afterwards by
   * `getLinearSolverOverride`. An empty value cancels the override.
   * \param[in] argc: number of arguments
   * \param[in] argv: arguments
   * \throws std::runtime_error on an invalid description, or if the option
   * is given twice or without value
   */
  void readLinearSolverOverride(const int argc,
                                const MainFunctionArguments argv);
#ifdef MFEM_USE_MPI
  /*!
   * \brief checks that all the processes of a communicator have the same
   * linear solver override (the environment variables are not always given
   * to the processes of the other nodes), collective
   * \param[in] comm: communicator
   * \throws std::runtime_error on all the processes if they differ
   */
  void checkLinearSolverOverride(const MPI_Comm comm);
#endif /* MFEM_USE_MPI */
  /*!
   * \return the linear solver which replaces the linear solvers of the non
   * linear evolution problems (see `setLinearSolver`), if any
   */
  [[nodiscard]] const std::optional<LinearSolverOverride>&
  getLinearSolverOverride() noexcept;
  /*!
   * \brief reports, once, on the output stream of MFEM, that the linear
   * solver override replaces a linear solver of the program
   * \param[in] replaced: description of the replaced linear solver
   */
  void reportLinearSolverOverride(std::string_view replaced) noexcept;
  /*!
   * \return the name and the parameters of the linear solver to build: the
   * linear solver override, if any (the replacement being reported, see
   * `reportLinearSolverOverride`), or the given ones
   * \param[in] n: name of the linear solver requested by the program
   * \param[in] p: parameters of the linear solver requested by the program
   */
  [[nodiscard]] std::pair<std::string_view, const Parameters&>
  selectLinearSolver(std::string_view n, const Parameters& p) noexcept;

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_LINEARSOLVEROVERRIDE_HXX */
