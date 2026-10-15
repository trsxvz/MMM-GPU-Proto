/*!
 * \file   include/MFEMMGIS/LinearSolverDescription.hxx
 * \brief  This file declares the `LinearSolverDescription` class, which
 * describes a linear solver or a preconditioner registered in the factories
 * of linear solvers and of preconditioners
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_LINEARSOLVERDESCRIPTION_HXX
#define LIB_MFEM_MGIS_LINEARSOLVERDESCRIPTION_HXX

#include <string>
#include <vector>
#include <optional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief where a linear solver or a preconditioner computes. It is compared
   * to the location of the linear solvers given by the execution settings
   * (see `ExecutionSettings::linear_solver`).
   */
  enum struct LinearSolverPlacement {
    /*!
     * \brief on the host only (direct solvers of the CPU, and the
     * preconditioners of hypre which only compute on the host)
     */
    HOST,
    //! \brief on the GPU only (AmgX, cuDSS)
    DEVICE,
    //! \brief where the device of MFEM computes (solvers of MFEM)
    MFEM,
    //! \brief where the memory of hypre is (solvers of hypre)
    HYPRE
  };  // end of LinearSolverPlacement

  //! \brief kind of linear solver
  enum struct LinearSolverKind {
    /*!
     * \brief iterative solver, converging up to its tolerances (Krylov
     * solvers)
     */
    ITERATIVE,
    //! \brief direct solver (factorization)
    DIRECT,
    /*!
     * \brief approximate solver, applying a fixed number of cycles or
     * sweeps of a preconditioner (for instance cycles of an algebraic
     * multigrid used by the fixed point solvers)
     */
    APPROXIMATE,
    //! \brief preconditioner (registered in the preconditioner factory)
    PRECONDITIONER
  };  // end of LinearSolverKind

  //! \brief preconditioners accepted by a linear solver
  enum struct AcceptedPreconditioners {
    //! \brief no preconditioner (direct solvers, self-contained solvers)
    NONE,
    //! \brief any preconditioner (seen as an `mfem::Solver`)
    ANY,
    /*!
     * \brief the preconditioners of the same library only (for instance,
     * the Krylov solvers of hypre only accept the preconditioners of hypre)
     */
    SAME_LIBRARY
  };  // end of AcceptedPreconditioners

  /*!
   * \brief description of a linear solver or of a preconditioner, given
   * when it is registered. The checks of the execution settings and of the
   * preconditioners given to a linear solver are based on this description.
   */
  struct LinearSolverDescription {
    //! \brief library providing the solver (`MFEM`, `hypre`, `MUMPS`, ...)
    std::string library;
    //! \brief kind of solver
    LinearSolverKind kind = LinearSolverKind::ITERATIVE;
    /*!
     * \brief where the solver computes. It is empty for a wrapper of another
     * solver, which computes where the wrapped solver does (for instance
     * the `AssembledOperator` preconditioner): the location of the wrapped
     * solver is then checked.
     */
    std::optional<LinearSolverPlacement> placement =
        LinearSolverPlacement::MFEM;
    /*!
     * \brief preconditioners accepted by a linear solver (meaningless for a
     * preconditioner). The linear solvers registered without description
     * accept any preconditioner (see `LinearSolverFactory::add`).
     */
    AcceptedPreconditioners preconditioners = AcceptedPreconditioners::NONE;
    //! \brief names of the options accepted by the solver
    std::vector<std::string> options;
    //! \brief short description of the solver
    std::string description;
  };  // end of LinearSolverDescription

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_LINEARSOLVERDESCRIPTION_HXX */
