/*!
 * \file   MFEMMGIS/AbstractBehaviourIntegrator.hxx
 * \brief  This file declares the `AbstractBehaviourIntegrator` class
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#ifndef LIB_MFEM_MGIS_ABSTRACTBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_ABSTRACTBEHAVIOURINTEGRATOR_HXX

#include <span>
#include <array>
#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#include "MFEMMGIS/Batch/InnerForces.hxx"
#include "MFEMMGIS/Batch/StiffnessMatrices.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Material;
  struct PartialQuadratureSpace;
  struct ImmutablePartialQuadratureFunctionView;
  struct AbstractQPEvaluator;
  enum struct IntegrationType;

  /*!
   * \brief abstract class for all behaviour integrators
   *
   * This class provides methods to:
   *
   * - integrate the behaviour over the time step
   * - compute the nodal forces due to the material reaction (see the
   *   `updateResidual` method).
   * - compute the stiffness matrix (see the `updateJacobian` method).
   */
  struct MFEM_MGIS_EXPORT AbstractBehaviourIntegrator {
    //! \return the current time increment
    virtual real getTimeIncrement() const noexcept = 0;
    /*!
     * \brief set the time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real dt) = 0;
    //! \return the partial quadrature space
    virtual const PartialQuadratureSpace &getPartialQuadratureSpace()
        const noexcept = 0;
    /*!
     * \brief return the integration rule for the given element and element
     * transformation
     * \param[in] e: element
     * \param[in] tr: element transformation
     * \return the integration rule
     */
    virtual const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &e,
        const mfem::ElementTransformation &tr) const = 0;
    /*!
     * \brief return the weight of the integration point, taking the
     * modelling hypothesis into account
     * \param[in] tr: element transformation
     * \param[in] ip: integration point
     * \return the weight of the integration point
     */
    virtual real getIntegrationPointWeight(
        mfem::ElementTransformation &tr,
        const mfem::IntegrationPoint &ip) const = 0;
    /*!
     * \brief method called at the beginning of each resolution
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] virtual bool setup(Context &ctx,
                                     const real t,
                                     const real dt) noexcept = 0;
    /*!
     * \brief integrate the behaviour over the time step
     *
     * If successful, the values of the thermodynamic forces, consistent
     * tangent operator and internal state variables are updated.
     *
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \param[in] it: integration type
     * \return true on success
     */
    [[nodiscard]] virtual bool integrate(const mfem::FiniteElement &e,
                                         mfem::ElementTransformation &tr,
                                         const mfem::Vector &u,
                                         const IntegrationType it) = 0;
    /*!
     * \return true if the behaviour integrator supports the integration of
     * the behaviour by batches (see `computeGradients` and
     * `integrateBatch`). The default implementation returns false.
     */
    [[nodiscard]] virtual bool supportsBatchIntegration() const noexcept;
    /*!
     * \brief compute the gradients at the integration points of an element,
     * first step of the integration by batches: the operations of
     * `integrate` which precede the integration of the behaviour at each
     * point (including the update of the external state variables computed
     * from the unknowns), the gradients being left in the global frame.
     *
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \return true on success
     * \note the default implementation returns false
     */
    [[nodiscard]] virtual bool computeGradients(const mfem::FiniteElement &e,
                                                mfem::ElementTransformation &tr,
                                                const mfem::Vector &u);
    /*!
     * \brief integrate the behaviour at all the integration points, as a
     * batch, after the computation of the gradients on all the elements
     * (see `computeGradients`), or without it for a prediction. For an
     * orthotropic behaviour, the gradients are rotated in the material frame
     * before the integration, and the tangent operator blocks in the global
     * frame after it, at all the points at once.
     *
     * For a behaviour generated by the `generic-parallel` interface, all the
     * points are integrated, even if the integration fails at some of them
     * (`integrate` stops at the first failure). For a behaviour generated by
     * the `generic` interface, MGIS also stops at the first failure.
     *
     * \param[in] it: integration type
     * \return true on success
     * \note the default implementation returns false
     */
    [[nodiscard]] virtual bool integrateBatch(const IntegrationType it);
    /*!
     * \return true if `integrateBatch` integrates the behaviour on the
     * copies of the arrays of the material on the device, which are then
     * the valid ones (see `MaterialArrays`)
     * \note the default implementation returns false
     */
    [[nodiscard]] virtual bool integratesBatchesOnDevice() const noexcept;
    /*!
     * \return true if the gradients, the inner forces and the stiffness
     * matrices of the behaviour integrator can be computed by kernels (see
     * `computeGradientsByKernels`, `addInnerForcesByKernels` and
     * `addStiffnessMatricesByKernels`). The default implementation returns
     * false.
     */
    [[nodiscard]] virtual bool supportsKernels() const noexcept;
    /*!
     * \return true if the computation of the gradients by kernels requires
     * the values of the unknowns at the integration points (see
     * `batch::UnknownsDerivatives::point_values`), for instance the
     * temperature of the stationary heat transfer behaviours, an external
     * state variable of the behaviour, or the temperature of the transient
     * term of the heat equation, the gradient of the behaviour. The default
     * implementation returns false.
     */
    [[nodiscard]] virtual bool requiresUnknownsValuesForKernels()
        const noexcept;
    /*!
     * \brief compute the gradients at all the integration points by
     * kernels, from the derivatives of the unknowns at the points of the
     * integration rule of the elements (see
     * `PartialQuadratureSpace::getUniformIntegrationRule`): first step of
     * the integration by batches, replacing the calls to `computeGradients`
     * element by element. As with `computeGradients`, the gradients are left
     * in the global frame.
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] d: derivatives of the unknowns, which must contain the ones
     * of the integration rule of the elements
     * \param[in] g: geometric data, which must contain the ones of the
     * integration rule of the elements
     * \note the default implementation reports an error
     */
    [[nodiscard]] virtual bool computeGradientsByKernels(
        Context&,
        const batch::UnknownsDerivativesSet&,
        const batch::ElementsGeometrySet&) noexcept;
    /*!
     * \brief add the inner forces of all the elements, computed by kernels
     * from the thermodynamic forces at the end of the time step, to values
     * element by element in the layout of the element restriction of MFEM
     * (see `batch::addInnerForces`)
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in, out] f: values element by element
     * \param[in] g: geometric data, which must contain the ones of the
     * integration rule of the elements
     * \note the default implementation reports an error
     */
    [[nodiscard]] virtual bool addInnerForcesByKernels(
        Context&, mfem::Vector&, const batch::ElementsGeometrySet&) noexcept;
    /*!
     * \brief add the stiffness matrices of all the elements, computed by
     * kernels from the tangent operator blocks, to matrices stored element by
     * element (see `batch::addStiffnessMatrices`)
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in, out] m: matrices element by element
     * \param[in] g: geometric data, which must contain the ones of the
     * integration rule of the elements
     * \note the default implementation reports an error
     */
    [[nodiscard]] virtual bool addStiffnessMatricesByKernels(
        Context&, mfem::Vector&, const batch::ElementsGeometrySet&) noexcept;
    /*!
     * \brief add the stiffness matrices of all the elements, computed by
     * kernels from the tangent operator blocks, to an assembled sparse
     * matrix, without storing them (see `batch::addStiffnessMatrices`)
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in, out] m: assembled matrix
     * \param[in] g: geometric data
     * \note the default implementation reports an error
     */
    [[nodiscard]] virtual bool addStiffnessMatricesByKernels(
        Context&,
        const batch::AssembledStiffnessMatrix&,
        const batch::ElementsGeometrySet&) noexcept;
    /*!
     * \brief add the action of the stiffness matrices of all the elements,
     * computed by kernels from the tangent operator blocks without
     * assembling them, on increments of the unknowns given element by
     * element: `Σ Bᵀ K B x w det J`, in the layout of the element restriction
     * of MFEM (see `addInnerForcesByKernels`)
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in, out] y: values element by element
     * \param[in] x: increments of the unknowns element by element (layout
     * of the element restriction of MFEM)
     * \param[in] g: geometric data
     * \note the default implementation reports an error
     */
    [[nodiscard]] virtual bool addTangentOperatorActionByKernels(
        Context&,
        mfem::Vector&,
        const mfem::Vector&,
        const batch::ElementsGeometrySet&) noexcept;
    /*!
     * \brief add the diagonals of the stiffness matrices of all the elements,
     * computed by kernels, in the layout of the element restriction of MFEM
     * (see `addInnerForcesByKernels`)
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in, out] d: values element by element
     * \param[in] g: geometric data
     * \note the default implementation reports an error
     */
    [[nodiscard]] virtual bool addTangentOperatorDiagonalByKernels(
        Context&, mfem::Vector&, const batch::ElementsGeometrySet&) noexcept;
    /*!
     * \brief compute the contribution of the given element to the inner forces
     * \param[out] Fe: inner forces
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     */
    virtual void computeInnerForces(mfem::Vector &Fe,
                                    const mfem::FiniteElement &e,
                                    mfem::ElementTransformation &tr) = 0;
    /*!
     * \brief compute the contribution of the given element to the residual
     * \param[out] Fe: element contribution to the residual
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     */
    virtual void updateResidual(mfem::Vector &Fe,
                                const mfem::FiniteElement &e,
                                mfem::ElementTransformation &tr,
                                const mfem::Vector &u) = 0;
    /*!
     * \brief compute the contribution of the given element to the jacobian
     * \param[out] Ke: element stiffness matrix
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     */
    virtual void updateJacobian(mfem::DenseMatrix &Ke,
                                const mfem::FiniteElement &e,
                                mfem::ElementTransformation &tr,
                                const mfem::Vector &u) = 0;
    /*!
     * \brief clean-up at the end of a resolution
     *
     * \param[in, out] ctx: execution context
     * \return true on success
     *
     * \note this method is the counterpart of the `setup` method and is meant
     * to release any memory allocated by this method.
     */
    [[nodiscard]] virtual bool cleanup(Context &ctx) noexcept = 0;
    /*!
     * \brief revert the internal state variables.
     *
     * The values of the internal state variables at the beginning of the time
     * step are copied on the values of the internal state variables at
     * end of the time step.
     */
    virtual void revert() = 0;
    /*!
     * \brief update the internal state variables.
     *
     * The values of the internal state variables at the end of the time step
     * are copied on the values of the internal state variables at beginning of
     * the time step.
     */
    virtual void update() = 0;
    /*!
     * \return if the call to getMaterial is valid
     *
     * This has been introduced to be able to build behaviour integrators not
     * built on MGIS and MFront.
     */
    virtual bool hasMaterial() const noexcept = 0;
    /*!
     * \brief return the underlying material
     * \param[in, out] ctx: execution context
     * \return the underlying material, if any
     */
    [[nodiscard]] virtual OptionalReference<Material> getMaterial(
        Context &ctx) noexcept = 0;
    /*!
     * \brief return the underlying material
     * \param[in, out] ctx: execution context
     * \return the underlying material, if any
     */
    [[nodiscard]] virtual OptionalReference<const Material> getMaterial(
        Context &ctx) const noexcept = 0;
    /*!
     * \brief set the macroscopic gradients
     * \param[in] g: macroscopic gradients
     */
    virtual void setMacroscopicGradients(std::span<const real> g) = 0;
    /*!
     * \return if the current solution is required for assembling the residual
     * \note this is required when creating a linear operator evaluating the
     * residual
     */
    [[nodiscard]] virtual bool requiresCurrentSolutionForResidualAssembly()
        const noexcept = 0;
    /*!
     * \return if the current solution is required for assembling the jacobian
     * \note this is required when creating a linear operator evaluating the
     * jacobian
     */
    [[nodiscard]] virtual bool requiresCurrentSolutionForJacobianAssembly()
        const noexcept = 0;
    /*!
     * \brief set the value of a material property
     *
     * \param[in, out] ctx: execution context
     * \param[in] name: name of the material property
     * \param[in] e: evaluator of the material property
     * \param[in] ts: time step stage
     * \return true on success
     */
    [[nodiscard]] virtual bool setMaterialProperty(
        Context &ctx,
        std::string_view name,
        std::shared_ptr<const AbstractQPEvaluator> e,
        const TimeStepStage ts) noexcept = 0;
    /*!
     * \brief set the value of an external state variable
     *
     * \param[in, out] ctx: execution context
     * \param[in] name: name of the external state variable
     * \param[in] e: evaluator of the external state variable
     * \param[in] ts: time step stage
     * \return true on success
     */
    [[nodiscard]] virtual bool setExternalStateVariable(
        Context &ctx,
        std::string_view name,
        std::shared_ptr<const AbstractQPEvaluator> e,
        const TimeStepStage ts) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractBehaviourIntegrator();
  };  // end of struct AbstractBehaviourIntegrator

  /*!
   * \brief return the measure of the mesh on which the behaviour integrator is
   * built: volume in 3D and axisymmetrical hypotheses, surface in other
   * bidimensional hypotheses.
   * \param[in] bi: behaviour integrator
   * \return the measure of the mesh
   */
  MFEM_MGIS_EXPORT real computeMeasure(const AbstractBehaviourIntegrator &bi);

  /*!
   * \brief return the integral of a partial quadrature function
   * \param[in] bi: behaviour integrator
   * \param[in] f: function
   * \return the integral of the function
   */
  template <typename ValueType>
  ValueType computeIntegral(const AbstractBehaviourIntegrator &bi,
                            const ImmutablePartialQuadratureFunctionView &f);
  /*!
   * \return the integral of a partial quadrature function
   * \param[in] bi: behaviour integrator
   * \param[in] f: function
   */
  template <>
  MFEM_MGIS_EXPORT real
  computeIntegral(const AbstractBehaviourIntegrator &bi,
                  const ImmutablePartialQuadratureFunctionView &f);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACTBEHAVIOURINTEGRATOR_HXX */
