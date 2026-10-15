/*!
 * \file   include/MFEMMGIS/StandardBehaviourIntegratorCRTPBase.hxx
 * \brief  This file declares the `StandardBehaviourIntegratorCRTPBase` class
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#ifndef LIB_MFEM_MGIS_ISOTROPICSTANDARDBEHAVIOURINTEGRATORCRTPBASE_HXX
#define LIB_MFEM_MGIS_ISOTROPICSTANDARDBEHAVIOURINTEGRATORCRTPBASE_HXX

#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  /*!
   * \brief a base class for standard behaviour integrators, factorizing code
   * by static polymorphism using the CRTP idiom.
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
   *   one integration point, and a method called
   *   `setExternalStateVariablesFromUnknownsValues`, which sets them at all
   *   the integration points (computations by kernels)
   *
   * The thermodynamic forces of an orthotropic behaviour without functions
   * rotating the values of one point (`generic-parallel` interface) are
   * rotated in the global frame at all the points at once (see
   * `BehaviourIntegratorBase::getRotatedThermodynamicForces`).
   */
  template <typename Child>
  struct StandardBehaviourIntegratorCRTPBase : BehaviourIntegratorBase {
    /*!
     * \return true: the behaviour can be integrated by batches (see
     * `AbstractBehaviourIntegrator::supportsBatchIntegration`)
     */
    [[nodiscard]] bool supportsBatchIntegration() const noexcept override;
    /*!
     * \brief compute the gradients at the integration points of an element
     * (see `AbstractBehaviourIntegrator::computeGradients` and
     * `implementComputeGradients`)
     */
    [[nodiscard]] bool computeGradients(const mfem::FiniteElement&,
                                        mfem::ElementTransformation&,
                                        const mfem::Vector&) override;
    /*!
     * \return true if the kind of the gradients of the behaviour matches the
     * computations of the behaviour integrator (see
     * `BehaviourIntegratorTraits`), if neither the residual nor the jacobian
     * depend on the current solution, and if the quadrature space is
     * supported by the kernels (see `isQuadratureSpaceSupportedByKernels`)
     */
    [[nodiscard]] bool supportsKernels() const noexcept override;
    /*!
     * \return true if the external state variables are updated from the
     * values of the unknowns or if the gradients are computed from the
     * values of the shape functions (see `BehaviourIntegratorTraits`)
     */
    [[nodiscard]] bool requiresUnknownsValuesForKernels()
        const noexcept override;
    /*!
     * \brief compute the gradients at all the integration points by kernels
     * (see `batch::computeGradients`) and, if the external state variables
     * are updated from the values of the unknowns, set them (see
     * `AbstractBehaviourIntegrator::computeGradientsByKernels`)
     */
    [[nodiscard]] bool computeGradientsByKernels(
        Context& ctx,
        const batch::UnknownsDerivativesSet& d,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the inner forces of all the elements, computed by kernels
     * from the thermodynamic forces (see
     * `AbstractBehaviourIntegrator::addInnerForcesByKernels` and
     * `addInnerForcesFromThermodynamicForcesByKernels`)
     */
    [[nodiscard]] bool addInnerForcesByKernels(
        Context& ctx,
        mfem::Vector& f,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the stiffness matrices of all the elements, computed by
     * kernels, to matrices stored element by element (see
     * `AbstractBehaviourIntegrator::addStiffnessMatricesByKernels` and
     * `batch::addStiffnessMatrices`)
     */
    [[nodiscard]] bool addStiffnessMatricesByKernels(
        Context& ctx,
        mfem::Vector& m,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the stiffness matrices of all the elements, computed by
     * kernels, to an assembled sparse matrix (see
     * `AbstractBehaviourIntegrator::addStiffnessMatricesByKernels` and
     * `batch::addStiffnessMatrices`)
     */
    [[nodiscard]] bool addStiffnessMatricesByKernels(
        Context& ctx,
        const batch::AssembledStiffnessMatrix& m,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the action of the stiffness matrices of all the elements,
     * computed by kernels without assembling them (see
     * `AbstractBehaviourIntegrator::addTangentOperatorActionByKernels` and
     * `batch::addTangentOperatorAction`)
     */
    [[nodiscard]] bool addTangentOperatorActionByKernels(
        Context& ctx,
        mfem::Vector& y,
        const mfem::Vector& x,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the diagonals of the stiffness matrices of all the
     * elements, computed by kernels (see
     * `AbstractBehaviourIntegrator::addTangentOperatorDiagonalByKernels` and
     * `batch::addStiffnessMatricesDiagonals`)
     */
    [[nodiscard]] bool addTangentOperatorDiagonalByKernels(
        Context& ctx,
        mfem::Vector& d,
        const batch::ElementsGeometrySet& g) noexcept override;
    /*!
     * \brief add the nodal blocks of the stiffness matrices of all the
     * elements, computed by kernels (see
     * `AbstractBehaviourIntegrator::addTangentOperatorNodalBlocksByKernels`
     * and `batch::addStiffnessMatricesNodalBlocks`)
     */
    [[nodiscard]] bool addTangentOperatorNodalBlocksByKernels(
        Context& ctx,
        mfem::Vector& b,
        const batch::ElementsGeometrySet& g) noexcept override;

   protected:
    /*!
     * \brief values of the unknowns at the integration points of the
     * material, computed by kernels (workspace, see
     * `requiresUnknownsValuesForKernels`)
     */
    mfem::Vector unknowns_values;
    /*!
     * \brief implementation of `addStiffnessMatricesByKernels`
     * \tparam Destination: matrices element by element (`mfem::Vector`) or
     * assembled matrix (`batch::AssembledStiffnessMatrix`)
     */
    template <typename Destination>
    [[nodiscard]] bool implementAddStiffnessMatricesByKernels(
        Context&, Destination&, const batch::ElementsGeometrySet&) noexcept;
    /*!
     * \brief constructor
     * \param[in] s: quadrature space
     * \param[in] b_ptr: behaviour
     */
    StandardBehaviourIntegratorCRTPBase(
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::unique_ptr<const Behaviour> b_ptr);
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
     * \brief compute the gradients at an integration point in the global
     * frame: macroscopic gradients, contributions of the nodes, external
     * state variables computed from the unknowns
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \param[in] ip: integration point
     * \param[in] o: offset of the integration point
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
        mfem::Vector& shape_functions,
        mfem::DenseMatrix& shape_functions_derivatives);
    /*!
     * \brief compute the contribution of the element to the residual
     * \param[out] Fe: element contribution to the residual
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns, unused
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
     * \param[out] Ke: element stiffness matrix
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     *
     * \note The implementation of the `updateJacobian` in the
     * `Child` class trivially calls this method. This indirection is made to
     * control where the code associated to the
     * `implementUpdateJacobian` is generated.
     */
    void implementUpdateJacobian(mfem::DenseMatrix& Ke,
                                 const mfem::FiniteElement& e,
                                 mfem::ElementTransformation& tr);
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
    ~StandardBehaviourIntegratorCRTPBase() override;

#ifndef MFEM_THREAD_SAFE
   protected:
    //! \brief vector used to store the value of the shape functions
    mfem::Vector shape;
    //! \brief matrix used to store the derivatives of the shape functions
    mfem::DenseMatrix dshape;
#endif

  };  // end of StandardBehaviourIntegratorCRTPBase

}  // end of namespace mfem_mgis

#include "MFEMMGIS/StandardBehaviourIntegratorCRTPBase.ixx"

#endif /* LIB_MFEM_MGIS_ISOTROPICSTANDARDBEHAVIOURINTEGRATORCRTPBASE_HXX */
