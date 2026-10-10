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
   *
   * Linear solver of `K0` (option `ElasticOperatorLinearSolver`, `Name` and
   * `Options` as the linear solvers of the problems, for instance
   * `HypreBoomerAMG` for one cycle of BoomerAMG per iteration): owned by
   * the solver, so that its setup is kept between resolutions with
   * `ElasticOperatorUpdate` set to `Once`. If not given, the linear solver
   * of the problem is used, and set up with `K0` at each resolution (the
   * prediction of the problem gives it its own operator).
   *
   * Acceleration (option `Acceleration`, a name or a `Name` and `Options`,
   * see `AccelerationAlgorithmFactory`): none by default.
   *
   * Linear solver of the elastic prediction of the problem (option
   * `PredictionLinearSolver`): `Problem` (default), the linear solver of
   * the problem, or `ElasticOperator`, the conjugate gradient of MFEM
   * preconditioned by the linear solver of the elastic operator, whose
   * setup is kept, with the options `RelativeTolerance` (1e-8 by default)
   * and `MaximumNumberOfIterations` (1000 by default). The latter is used
   * by the predictions which follow the first resolution, with the
   * `ElasticOperator` method and a linear solver of the elastic operator.
   *
   * The other options are the ones of the iterative solvers
   * (`RelativeTolerance`, `AbsoluteTolerance`, `MaximumNumberOfIterations`,
   * `VerbosityLevel`) and `DiscardLinearSolverFailure`. When the
   * integration of the behaviours fails at a new estimate, it is replaced
   * by the middle of this estimate and of the previous one, ten times at
   * most.
   */
  struct MFEM_MGIS_EXPORT FixedPointSolver final : NonLinearSolverBase {
    //! \brief methods
    enum struct Method { ELASTIC_OPERATOR, MODIFIED_NEWTON };
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    FixedPointSolver(NonLinearEvolutionProblemImplementation<true>&);
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    FixedPointSolver(NonLinearEvolutionProblemImplementation<false>&);
    [[nodiscard]] bool setSolverParameters(Context&,
                                           const Parameters&) noexcept override;
    void Mult(const mfem::Vector&, mfem::Vector&) const override;
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
    //! \brief integration of the behaviours
    std::function<bool(const mfem::Vector&, const IntegrationType)>
        integrate_behaviours;
    //! \brief generator of the linear solver of the elastic operator
    std::function<LinearSolverHandler(
        Context&, std::string_view, const Parameters&)>
        generate_linear_solver;
    //! \brief method
    Method method = Method::ELASTIC_OPERATOR;
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
     * \brief if true, the elastic predictions are solved by the conjugate
     * gradient preconditioned by the linear solver of the elastic operator
     */
    bool prediction_by_elastic_operator = false;
    //! \brief relative tolerance of the elastic predictions
    real prediction_relative_tolerance = 1e-8;
    //! \brief maximum number of iterations of the elastic predictions
    int prediction_maximum_number_of_iterations = 1000;
    //! \brief linear solver of the elastic predictions
    std::unique_ptr<mfem::IterativeSolver> prediction_linear_solver;
    /*!
     * \brief preconditioner of the elastic predictions, applying the linear
     * solver of the elastic operator
     */
    std::unique_ptr<mfem::Solver> prediction_preconditioner;
  };  // end of struct FixedPointSolver

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEARSOLVERS_FIXEDPOINTSOLVER_HXX */
