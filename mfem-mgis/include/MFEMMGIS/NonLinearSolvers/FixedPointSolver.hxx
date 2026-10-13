/*!
 * \file   MFEMMGIS/NonLinearSolvers/FixedPointSolver.hxx
 * \brief  This file declares the `FixedPointSolver` class
 * \date   10/10/2026
 */

#ifndef LIB_MFEM_MGIS_NONLINEARSOLVERS_FIXEDPOINTSOLVER_HXX
#define LIB_MFEM_MGIS_NONLINEARSOLVERS_FIXEDPOINTSOLVER_HXX

#include <memory>
#include <string_view>
#include <functional>
#include "mfem/linalg/operator.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/LinearSolverHandler.hxx"
#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverBase.hxx"
#include "MFEMMGIS/NonLinearSolvers/AccelerationAlgorithm.hxx"
#include "MFEMMGIS/NonLinearSolvers/AdaptiveLinearSolverTolerance.hxx"

namespace mfem_mgis {

  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

  /*!
   * \brief fixed point iterations `x_{i+1} = x_i - S r(x_i)`, `r` being the
   * residual, possibly accelerated (see `AccelerationAlgorithm`), declared
   * as `FixedPoint` in the factory of the non linear solvers.
   *
   * Methods (option `Method`):
   *
   * - `ElasticOperator` (default): `S` is an exact or approximate inverse of
   *   the elastic operator `K0` (rigidity of the initial state, as in
   *   Cast3M). `K0` is assembled from the elastic operators of the
   *   behaviours (option `ElasticOperatorComputation`: `Prediction`, by
   *   default, uses the prediction operator of the behaviours;
   *   `Integration` integrates the behaviours and uses their elastic
   *   operator, for the behaviours which do not provide the prediction
   *   operator), at each resolution or at the first one only (option
   *   `ElasticOperatorUpdate`: `EachResolution`, by default, or `Once`).
   *   The behaviours are integrated without tangent operator at each
   *   iteration.
   * - `ModifiedNewton`: `S` is the inverse of the consistent tangent
   *   operator at the first iteration of each resolution, given to the
   *   linear solver of the problem.
   * - `NewtonK0`: Newton's method, the consistent tangent operator being
   *   inverted by an iterative solver of MFEM preconditioned by the linear
   *   solver of the elastic operator (option `ElasticOperatorLinearSolver`,
   *   required), whose setup is kept: for instance BoomerAMG built on `K0`
   *   instead of a new setup at each iteration. `K0` is computed as with the
   *   `ElasticOperator` method (options `ElasticOperatorComputation` and
   *   `ElasticOperatorUpdate`): with `EachResolution`, the setup of its
   *   linear solver is made at each resolution, with `Once` at the first
   *   one only. The iterative solver is given by the option
   *   `TangentLinearSolver` (see below). The tangent operator may be applied
   *   without assembly (parameter `TangentOperatorComputation` of the
   *   problem set to `MatrixFree`).
   *
   * Linear solver of `K0` (option `ElasticOperatorLinearSolver`, a linear
   * solver of the factory of the linear solvers, `LinearSolverFactory`, for
   * instance `HypreBoomerAMG` for one cycle of BoomerAMG per iteration in a
   * parallel computation):
   * owned by the solver, so that its setup is kept between resolutions with
   * `ElasticOperatorUpdate` set to `Once`. The value `None` removes it. If
   * not given, the linear solver of the problem is used, and set up with
   * `K0` at each resolution (the prediction of the problem gives it its own
   * operator).
   *
   * Acceleration (option `Acceleration`, see `AccelerationAlgorithmFactory`):
   * none by default. The value `None` removes it. The linear solver of the
   * elastic operator and the acceleration are given by a name, by a
   * dictionary associating a name with its options, or by a dictionary with
   * the entries `Name` and `Options`. The linear solver of the elastic
   * operator is refused with the `ModifiedNewton` method, and the
   * acceleration with the `NewtonK0` method.
   *
   * Iterative solver of the tangent operator of the `NewtonK0` method
   * (option `TangentLinearSolver`, ignored by the other methods): an
   * iterative solver of MFEM declared in the factory of the linear solvers
   * (`CGSolver`, by default, `GMRESSolver`, `FGMRESSolver`,
   * `BiCGSTABSolver`, `MINRESSolver` or `SLISolver`), given as the linear
   * solver of the elastic operator, with its options except
   * `Preconditioner`. The `RelativeTolerance`, the `AbsoluteTolerance` and
   * the `MaximumNumberOfIterations` are 1e-4, 0 and 1000 if not given. The
   * conjugate gradient requires a symmetric positive definite tangent
   * operator: `GMRESSolver` or `BiCGSTABSolver` are suited to the other
   * ones (for instance non associated behaviours). The warnings of the
   * iterative solvers (for instance an operator which is not positive
   * definite) are printed if their `VerbosityLevel` is given (0: warnings
   * and errors, 1: each iteration).
   *
   * Relative tolerance of the iterative solver of the tangent operator of
   * the `NewtonK0` method adapted to the decrease of the residual (inexact
   * Newton method, option `AdaptiveLinearSolverTolerance`, see
   * `AdaptiveLinearSolverTolerance`): refused with the other methods. The
   * tolerance given by `TangentLinearSolver` is restored at the end of each
   * resolution.
   *
   * Linear solver of the elastic prediction of the problem (option
   * `PredictionLinearSolver`): `Problem` (default), the linear solver of
   * the problem, or an iterative solver of MFEM preconditioned by the
   * linear solver of the elastic operator, given as for
   * `TangentLinearSolver` (`ElasticOperator` being a synonym of
   * `CGSolver`), with a `RelativeTolerance`, an `AbsoluteTolerance` and a
   * `MaximumNumberOfIterations` of 1e-8, 0 and 1000 if not given. The latter
   * is used by the elastic predictions made at the beginning of the time
   * steps which follow the first resolution, with the `ElasticOperator`
   * and `NewtonK0` methods and a linear solver of the elastic operator;
   * otherwise, the linear solver of the problem is used. The linear solver
   * of the elastic operator shall be a fixed linear operator (a direct
   * solver, a fixed number of cycles of BoomerAMG...) when it preconditions
   * an iterative solver other than `FGMRESSolver`: an iterative solver with
   * a tolerance makes the preconditioner vary.
   *
   * The other options are the ones of the iterative solvers
   * (`RelativeTolerance`, `AbsoluteTolerance`, `MaximumNumberOfIterations`,
   * `VerbosityLevel`) and `DiscardLinearSolverFailure`. The iterations stop
   * if the norm of the residual is not finite. When the
   * integration of the behaviours fails at a new estimate, it is replaced
   * by the middle of this estimate and of the previous one, ten times at
   * most.
   */
  struct MFEM_MGIS_EXPORT FixedPointSolver final : NonLinearSolverBase {
    //! \brief methods
    enum struct Method { ELASTIC_OPERATOR, MODIFIED_NEWTON, NEWTON_K0 };
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    FixedPointSolver(NonLinearEvolutionProblemImplementation<true>& p);
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    FixedPointSolver(NonLinearEvolutionProblemImplementation<false>& p);
    [[nodiscard]] bool setSolverParameters(Context&,
                                           const Parameters&) noexcept override;
    /*!
     * \brief solve the non linear problem
     * \param[in] b: right hand side, unused
     * \param[in, out] x: initial guess, then solution
     */
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override;
    [[nodiscard]] LinearSolver* getElasticPredictionLinearSolver() noexcept
        override;
    //! \brief destructor
    ~FixedPointSolver() noexcept override;

   private:
    /*!
     * \brief integrates the behaviours at an estimate of the unknowns, then
     * performs the actions associated with a new estimate
     * \return true on success
     * \param[in] u: estimate of the unknowns
     * \param[in] it: integration type
     */
    [[nodiscard]] bool integrate(const mfem::Vector&,
                                 const IntegrationType) const;
    /*!
     * \brief assembles the elastic operator at an estimate of the unknowns
     * and sets up its linear solver, if required
     * \return true on success
     * \param[in] u: estimate of the unknowns
     */
    [[nodiscard]] bool setElasticOperator(const mfem::Vector&) const;
    //! \return the linear solver giving the corrections
    [[nodiscard]] mfem::Solver& getCorrectionSolver() const;
    /*!
     * \return an iterative solver of MFEM given by the factory of the linear
     * solvers, to be preconditioned by the linear solver of the elastic
     * operator (see `getPreconditionedSolver`), or an invalid value on error
     * \param[in, out] ctx: execution context
     * \param[in] o: name of the option, for the error messages
     * \param[in] n: name of the iterative solver
     * \param[in] opts: options of the iterative solver
     * \param[in] tol: default relative tolerance
     */
    [[nodiscard]] LinearSolverHandler generateIterativeSolver(
        Context&,
        std::string_view,
        std::string_view,
        const Parameters&,
        const real) const noexcept;
    /*!
     * \return an iterative solver, after having given it the linear solver
     * of the elastic operator as preconditioner, whose setup is kept, if not
     * already done
     * \param[in] s: iterative solver
     * \param[in, out] p: preconditioner
     */
    [[nodiscard]] mfem::Solver& getPreconditionedSolver(
        const LinearSolverHandler&, std::unique_ptr<mfem::Solver>&) const;
    //! \brief integration of the behaviours
    std::function<bool(const mfem::Vector&, const IntegrationType)>
        integrate_behaviours;
    //! \brief generator of the linear solver of the elastic operator
    std::function<LinearSolverHandler(
        Context&, std::string_view, const Parameters&)>
        generate_linear_solver;
    //! \brief method
    Method method = Method::ELASTIC_OPERATOR;
    /*!
     * \brief settings of the relative tolerance of the solver of the
     * tangent operator adapted to the decrease of the residual (option
     * `AdaptiveLinearSolverTolerance`), disabled by default
     */
    AdaptiveLinearSolverTolerance adaptive_linear_solver_tolerance;
    //! \brief if true, the elastic operator is assembled once
    bool elastic_operator_assembled_once = false;
    /*!
     * \brief if true, the elastic operator is computed by an integration of
     * the behaviours, otherwise by their prediction operator
     */
    bool elastic_operator_by_integration = false;
    //! \brief acceleration algorithm, if any
    std::unique_ptr<AccelerationAlgorithm> acceleration;
    //! \brief linear solver of the elastic operator, if any
    LinearSolverHandler elastic_operator_linear_solver;
    //! \brief elastic operator
    mutable std::unique_ptr<mfem::Operator> elastic_operator;
    /*!
     * \brief iterative solver of the elastic predictions, preconditioned by
     * the linear solver of the elastic operator, if any
     */
    LinearSolverHandler prediction_linear_solver;
    //! \brief preconditioner of the elastic predictions
    std::unique_ptr<mfem::Solver> prediction_preconditioner;
    //! \brief iterative solver of the tangent operator (`NewtonK0`)
    LinearSolverHandler tangent_linear_solver;
    //! \brief preconditioner of the iterative solver of the tangent operator
    mutable std::unique_ptr<mfem::Solver> tangent_preconditioner;
  };  // end of struct FixedPointSolver

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEARSOLVERS_FIXEDPOINTSOLVER_HXX */
