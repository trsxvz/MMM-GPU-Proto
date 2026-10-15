/*!
 * \file   src/PETScSNESSolver.hxx
 * \brief  This file declares the `PETScSNESSolver` class, the non linear
 * solver based on the SNES of PETSc
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_PETSCSNESSOLVER_HXX
#define LIB_MFEM_MGIS_PETSCSNESSOLVER_HXX

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC)

#include <limits>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverBase.hxx"

namespace mfem_mgis::internals {

  // forward declaration
  struct PETScSNESSolverImplementation;

  /*!
   * \brief non linear solver based on the SNES of PETSc, through the
   * `mfem::PetscNonlinearSolver` class of MFEM, only available in parallel:
   * Newton's method with a backtracking line search, by default, or the
   * other methods of SNES (trust region, quasi-Newton, nonlinear Krylov
   * methods...), selected by the options of PETSc (option `Options`, for
   * instance `{"snes_type", "newtontr"}`, the options being given without
   * the leading dash).
   *
   * The convergence test of SNES is the one of mfem-mgis (the option
   * `snes_convergence_test` is refused), and the trust region methods,
   * which multiply the residual by the transposed tangent operator, require
   * an assembled tangent operator.
   *
   * Linear systems: by default (option `LinearSolver` set to `Problem`),
   * the tangent operator is inverted by the linear solver of the problem,
   * seen by PETSc as a shell preconditioner applied once (`ksp_type
   * preonly`, unless the options of PETSc give another Krylov solver,
   * which the linear solver of the problem then preconditions); the tangent
   * operator is seen by PETSc as a shell matrix, without copy (the line
   * search multiplies it by the Newton step). With `LinearSolver` set to
   * `PETSc`, the tangent operator, which must then be assembled, is
   * copied in the AIJ format of the host at each update, without
   * initializing the GPU (see `makeHostMatrix`), and the Krylov
   * solver and the preconditioner are the ones of PETSc given by the
   * options (`ksp_type`, `pc_type`, ...). With a GPU, SNES and the solvers
   * of PETSc compute on the host, the vectors being copied from and to the
   * vectors of mfem-mgis.
   *
   * The residual is the one of the problem. A failure of the integration of
   * the behaviours is reported to SNES as a domain error: the backtracking
   * line search then reduces its step, other methods stop. A failure of the
   * linear solver of the problem stops the resolution, unless the option
   * `DiscardLinearSolverFailure` is true.
   *
   * SNES stops when the euclidean norm of the residual is lower than
   * `max(RelativeTolerance * r_ref, AbsoluteTolerance)`, `r_ref` being the
   * reference norm of the residual (the norm of the residual of the
   * prediction, if any, or of the initial residual), as the Newton solver.
   * The resolution is stopped and accepted when the residual stagnates
   * (three iterations without a decrease of more than 10 %) below
   * `StagnationRelativeTolerance` times the reference norm. The other
   * criteria of SNES (length of the step, divergence) are not used.
   *
   * The behaviours are integrated again at the solution returned by SNES if
   * their last integration was made elsewhere (a trial point of the line
   * search), so that the state of the materials matches the solution.
   *
   * With `LinearSolver` set to `PETSc`, the tangent operator copied in the
   * AIJ format receives the rigid body modes as near null space and
   * the number of components as block size, for the algebraic multigrid of
   * PETSc, when the number of components of the unknowns is the space
   * dimension and the unknowns are ordered by vector dimension (option
   * `RigidBodyModes`, true by default, see `PETScRigidBodyModes`).
   *
   * Options: the ones of the iterative solvers (`RelativeTolerance`,
   * `AbsoluteTolerance`, `MaximumNumberOfIterations`, strictly positive,
   * `VerbosityLevel`), `DiscardLinearSolverFailure`, `RigidBodyModes`,
   * `LinearSolver`
   * (`Problem`, by default, or `PETSc`), `Options` (the options of PETSc,
   * for SNES, its line search and its Krylov solver; the options of the
   * preconditioner, `pc_...`, are refused with the linear solver of the
   * problem), `StagnationRelativeTolerance` (10^4 times the machine epsilon
   * by default, as for the Newton solver).
   */
  struct PETScSNESSolver final : public NonLinearSolverBase {
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    PETScSNESSolver(NonLinearEvolutionProblemImplementation<true>& p);
    /*!
     * \brief constructor, which throws: SNES is only available in parallel
     * \param[in] p: non linear evolution problem
     */
    PETScSNESSolver(NonLinearEvolutionProblemImplementation<false>& p);
    /*!
     * \brief solve the non linear problem
     * \param[in] b: right hand side, unused
     * \param[in, out] x: initial guess, then solution
     */
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override;
    [[nodiscard]] bool setSolverParameters(Context&,
                                           const Parameters&) noexcept override;
    //! \brief destructor
    ~PETScSNESSolver() override;

   private:
    friend struct PETScSNESSolverImplementation;
    friend struct ProblemLinearSolverPreconditioner;
    friend struct TangentOperatorForPETSc;
    //! \brief options of PETSc
    Parameters options;
    /*!
     * \brief if true, the linear systems are solved by the Krylov solver and
     * the preconditioner of PETSc, and by the linear solver of the problem
     * otherwise
     */
    bool petsc_linear_solver = false;
    /*!
     * \brief relative residual below which a stagnating resolution is
     * stopped and accepted (see the option `StagnationRelativeTolerance` of
     * `NewtonSolver`)
     */
    real stagnation_relative_tolerance =
        10000 * std::numeric_limits<real>::epsilon();
    /*!
     * \brief if true, the rigid body modes are given to the tangent
     * operator converted for the linear solvers of PETSc (see
     * `PETScRigidBodyModes`)
     */
    bool rigid_body_modes = true;
    //! \brief finite element space of the problem
    FiniteElementSpace<true>* const fespace = nullptr;
  };  // end of struct PETScSNESSolver

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC) */

#endif /* LIB_MFEM_MGIS_PETSCSNESSOLVER_HXX */
