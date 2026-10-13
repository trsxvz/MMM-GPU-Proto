/*!
 * \file   include/MFEMMGIS/FBarBehaviourIntegratorCRTPBase.hxx
 * \brief  This file declares the `FBarBehaviourIntegratorCRTPBase` class
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#ifndef LIB_MFEM_MGIS_ISOTROPICFBARBEHAVIOURINTEGRATORCRTPBASE_HXX
#define LIB_MFEM_MGIS_ISOTROPICFBARBEHAVIOURINTEGRATORCRTPBASE_HXX

#include <vector>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Batch/FBar.hxx"
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  /*!
   * \brief a base class for FBar behaviour integrators, factorizing code by
   * static polymorphism using the CRTP idiom.
   *
   * This class provides a way to optimise dynamic memory allocations.
   *
   * The `Child` class must provide:
   *
   * - a method called `getIntegrationRule`
   * - a method called `getIntegrationPointWeight`
   * - a method called `updateGradients`
   * - a method called `updateInnerForces`
   * - a method called `updateStiffnessMatrix`
   * - a method called `getRotationMatrix`
   * - a method called `rotateGradients`
   * - a method called `rotateThermodynamicForces`
   * - a method called `rotateTangentOperatorBlocks`
   * - if the external state variables are updated from the values of the
   *   unknowns (see `BehaviourIntegratorTraits`), a method called
   *   `updateExternalStateVariablesFromUnknownsValues`, which sets them at
   *   one integration point (the computations by kernels, which do not call
   *   it, are then not supported)
   *
   * The thermodynamic forces of an orthotropic behaviour without functions
   * rotating the values of one point (`generic-parallel` interface) are
   * rotated in the global frame at all the points at once (see
   * `BehaviourIntegratorBase::getRotatedThermodynamicForces`).
   *
   * The gradients, the inner forces and the stiffness matrices can also be
   * computed by kernels (see `batch::FBarData`): the deformation gradients
   * at the center of the elements are computed with the integration rule of
   * order 0 of their geometry (see `getAdditionalIntegrationRulesForKernels`).
   */
  template <typename Child, Hypothesis H>
  struct FBarBehaviourIntegratorCRTPBase : BehaviourIntegratorBase {
    /*!
     * \return true: the behaviour can be integrated by batches (see
     * `AbstractBehaviourIntegrator::supportsBatchIntegration`)
     */
    [[nodiscard]] bool supportsBatchIntegration() const noexcept override;
    /*!
     * \brief compute the gradients of the FBar formulation at the
     * integration points of an element (see
     * `AbstractBehaviourIntegrator::computeGradients`): the data of the FBar
     * formulation computed by kernels no longer match the gradients (see
     * `batch::FBarData::valid`)
     */
    [[nodiscard]] bool computeGradients(const mfem::FiniteElement&,
                                        mfem::ElementTransformation&,
                                        const mfem::Vector&) override;
    /*!
     * \brief integrate the behaviour at all the integration points, as a
     * batch (see `AbstractBehaviourIntegrator::integrateBatch`): the tangent
     * operator of the FBar formulation computed by kernels must then be
     * computed again
     */
    [[nodiscard]] bool integrateBatch(const IntegrationType it) override;
    /*!
     * \brief revert the state of the material and the data of the FBar
     * formulation computed by kernels to the beginning of the time step (see
     * `batch::revertFBarData`)
     */
    void revert() override;
    /*!
     * \brief update the state of the material and the data of the FBar
     * formulation computed by kernels at the end of the time step (see
     * `batch::updateFBarData`)
     */
    void update() override;
    /*!
     * \return true for the deformation gradients in 3D and in plane strain,
     * if the quadrature space is supported by the kernels (see
     * `isQuadratureSpaceSupportedByKernels`)
     */
    [[nodiscard]] bool supportsKernels() const noexcept override;
    //! \return the integration rule of the center of the elements
    [[nodiscard]] std::vector<const mfem::IntegrationRule*>
    getAdditionalIntegrationRulesForKernels() const override;
    /*!
     * \brief compute the deformation gradients at all the integration
     * points by kernels, then the gradients of the FBar formulation (see
     * `batch::applyFBar`)
     */
    [[nodiscard]] bool computeGradientsByKernels(
        Context& ctx,
        const batch::UnknownsDerivativesSet& d,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the inner forces of all the elements, computed by kernels
     * from the thermodynamic forces (see
     * `addInnerForcesFromThermodynamicForcesByKernels`)
     */
    [[nodiscard]] bool addInnerForcesByKernels(
        Context& ctx,
        mfem::Vector& f,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the stiffness matrices of the FBar formulation of all the
     * elements, computed by kernels, to matrices stored element by element
     * (see `batch::addFBarStiffnessMatrices`)
     */
    [[nodiscard]] bool addStiffnessMatricesByKernels(
        Context& ctx,
        mfem::Vector& m,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the stiffness matrices of the FBar formulation of all the
     * elements, computed by kernels, to an assembled sparse matrix (see
     * `batch::addFBarStiffnessMatrices`)
     */
    [[nodiscard]] bool addStiffnessMatricesByKernels(
        Context& ctx,
        const batch::AssembledStiffnessMatrix& m,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the action of the stiffness matrices of the FBar
     * formulation of all the elements, computed by kernels without
     * assembling them (see `batch::addFBarTangentOperatorAction`)
     */
    [[nodiscard]] bool addTangentOperatorActionByKernels(
        Context& ctx,
        mfem::Vector& y,
        const mfem::Vector& x,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the diagonals of the stiffness matrices of the FBar
     * formulation of all the elements, computed by kernels (see
     * `batch::addFBarStiffnessMatricesDiagonals`)
     */
    [[nodiscard]] bool addTangentOperatorDiagonalByKernels(
        Context& ctx,
        mfem::Vector& d,
        const batch::ElementsGeometrySet& g) noexcept override;

   protected:
    static_assert(((H == Hypothesis::PLANESTRAIN) ||
                   (H == Hypothesis::TRIDIMENSIONAL)),
                  "invalid modelling hypothesis");
    //! \brief data of the FBar formulation computed by kernels
    batch::FBarData fbar_data;
    /*!
     * \brief gets the data of the computations by kernels and computes the
     * tangent operator of the FBar formulation if it is not up to date (see
     * `batch::computeFBarTangentOperator`)
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[out] data: data of the computations by kernels
     * \param[in] name: name of the calling method, for the messages
     * \param[in] g: geometric data
     */
    [[nodiscard]] bool prepareFBarTangentOperatorForKernels(
        Context&,
        KernelsData&,
        const char* const,
        const batch::ElementsGeometrySet&) noexcept;
    /*!
     * \brief calls a function with the tangent operator blocks of the
     * behaviour in the global frame, readable on the device if a device is
     * enabled, the exceptions being reported in the execution context
     * \return the result of the function, false on exception
     * \param[in, out] ctx: execution context
     * \param[in] name: name of the calling method, for the messages
     * \param[in] f: function taking the tangent operator blocks
     */
    template <typename Function>
    [[nodiscard]] bool callWithTangentOperatorBlocks(Context&,
                                                     const char* const,
                                                     const Function&) noexcept;
    /*!
     * \brief constructor
     * \param[in] s: quadrature space
     * \param[in] b_ptr: behaviour
     */
    FBarBehaviourIntegratorCRTPBase(
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::unique_ptr<const Behaviour> b_ptr);
    /*!
     * \brief integrate the mechanical behaviour over the time step
     *
     * If successful, the values of the stress, consistent tangent
     * operator and internal state variables are updated.
     *
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \param[in] it: integration type
     * \return true on success
     */
    bool implementIntegrate(const mfem::FiniteElement& e,
                            mfem::ElementTransformation& tr,
                            const mfem::Vector& u,
                            const IntegrationType it);
    /*!
     * \brief compute the gradients at the integration points of an element:
     * the operations of `implementIntegrate` which precede the integration of
     * the behaviour, except the rotation of the gradients of an orthotropic
     * behaviour in the material frame, made at all the points by
     * `integrateBatch`
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \return true on success
     */
    bool implementComputeGradients(const mfem::FiniteElement& e,
                                   mfem::ElementTransformation& tr,
                                   const mfem::Vector& u);
    /*!
     * \return the determinant of the deformation gradient at the center of
     * the element
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \param[in, out] shape_functions: values of the shape functions
     * \param[in, out] shape_functions_derivatives: derivatives of the shape
     * functions
     */
    real computeDeformationGradientDeterminantAtCenter(
        const mfem::FiniteElement& e,
        mfem::ElementTransformation& tr,
        const mfem::Vector& u,
        mfem::Vector& shape_functions,
        mfem::DenseMatrix& shape_functions_derivatives);
    /*!
     * \brief compute the modified deformation gradient at an integration
     * point in the global frame (macroscopic gradients, contributions of the
     * nodes, external state variables computed from the unknowns, scaling
     * by the determinant at the center)
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \param[in] ip: integration point
     * \param[in] o: offset of the integration point
     * \param[in] J0: determinant of the deformation gradient at the center of
     * the element
     * \param[in, out] shape_functions: values of the shape functions
     * \param[in, out] shape_functions_derivatives: derivatives of the shape
     * functions
     */
    void computeGradientsAtIntegrationPoint(
        const mfem::FiniteElement& e,
        mfem::ElementTransformation& tr,
        const mfem::Vector& u,
        const mfem::IntegrationPoint& ip,
        const size_type o,
        const real J0,
        mfem::Vector& shape_functions,
        mfem::DenseMatrix& shape_functions_derivatives);
    /*!
     * \brief compute the contribution of the element to the residual
     * \param[out] Fe: element contribution to the residual
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimation of the displacement field
     *
     * \note Thanks to the CRTP idiom, this implementation can call
     * the `updateInnerForces` method defined
     * in the derived class without a virtual call. This call may
     * even be inlined.
     * \note The implementation of the `updateResidual` in the
     * `Child` class trivially calls this method. This indirection is made to
     * control where the code associated to the `implementUpdateResidual` is
     * generated.
     */
    void implementUpdateResidual(mfem::Vector& Fe,
                                 const mfem::FiniteElement& e,
                                 mfem::ElementTransformation& tr,
                                 const mfem::Vector& u);
    /*!
     * \brief compute the contribution of the element to the jacobian
     *
     * \param[out] Ke: element stiffness matrix
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     *
     * \note The implementation of the `updateJacobian` in the
     * `Child` class trivially calls this method. This indirection is made to
     * control where the code associated to the
     * `implementUpdateJacobian` is generated.
     */
    void implementUpdateJacobian(mfem::DenseMatrix& Ke,
                                 const mfem::FiniteElement& e,
                                 mfem::ElementTransformation& tr,
                                 const mfem::Vector& u);
    /*!
     * \brief compute the contribution of the element to the inner forces
     * \param[out] Fe: inner forces
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     */
    void implementComputeInnerForces(mfem::Vector& Fe,
                                     const mfem::FiniteElement& e,
                                     mfem::ElementTransformation& tr);
    //! \brief destructor
    ~FBarBehaviourIntegratorCRTPBase() override;

#ifndef MFEM_THREAD_SAFE
   protected:
    /*!
     * \brief vector used to store the value of the shape functions at the
     * center of the element
     */
    mfem::Vector shape0;
    /*!
     * \brief matrix used to store the derivatives of the shape functions at the
     * center of the element
     */
    mfem::DenseMatrix dshape0;
    //! \brief vector used to store the value of the shape functions
    mfem::Vector shape;
    //! \brief matrix used to store the derivatives of the shape functions
    mfem::DenseMatrix dshape;
#endif

  };  // end of FBarBehaviourIntegratorCRTPBase

}  // end of namespace mfem_mgis

#include "MFEMMGIS/FBarBehaviourIntegratorCRTPBase.ixx"

#endif /* LIB_MFEM_MGIS_ISOTROPICFBARBEHAVIOURINTEGRATORCRTPBASE_HXX */
