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
   *
   * The thermodynamic forces of an orthotropic behaviour without functions
   * rotating the values of one point (`generic-parallel` interface) are
   * rotated in the global frame at all the points at once (see
   * `BehaviourIntegratorBase::getRotatedThermodynamicForces`).
   */
  template <typename Child>
  struct StandardBehaviourIntegratorCRTPBase : BehaviourIntegratorBase {
    [[nodiscard]] bool supportsBatchIntegration() const noexcept override;
    [[nodiscard]] bool computeGradients(const mfem::FiniteElement&,
                                        mfem::ElementTransformation&,
                                        const mfem::Vector&) override;
    [[nodiscard]] bool supportsKernels() const noexcept override;
    [[nodiscard]] bool computeGradientsByKernels(
        Context&,
        const batch::UnknownsDerivativesSet&,
        const batch::ElementsGeometrySet&) noexcept override;
    [[nodiscard]] bool addInnerForcesByKernels(
        Context&,
        mfem::Vector&,
        const batch::ElementsGeometrySet&) noexcept override;
    [[nodiscard]] bool addStiffnessMatricesByKernels(
        Context&,
        mfem::Vector&,
        const batch::ElementsGeometrySet&) noexcept override;
    [[nodiscard]] bool addStiffnessMatricesByKernels(
        Context&,
        const batch::AssembledStiffnessMatrix&,
        const batch::ElementsGeometrySet&) noexcept override;

    [[nodiscard]] bool addTangentOperatorActionByKernels(
        Context&,
        mfem::Vector&,
        const mfem::Vector&,
        const batch::ElementsGeometrySet&) noexcept override;
    [[nodiscard]] bool addTangentOperatorDiagonalByKernels(
        Context&,
        mfem::Vector&,
        const batch::ElementsGeometrySet&) noexcept override;

   protected:
    //! \brief data of the computations by kernels on the material
    struct KernelsData {
      //! \brief kind of the gradients
      batch::GradientsKind kind = batch::GradientsKind::SMALL_STRAIN_3D;
      //! \brief integration rule of the material, null if it has no element
      const mfem::IntegrationRule* rule = nullptr;
      //! \brief geometric data of this rule
      const batch::ElementsGeometry* geometry = nullptr;
    };
    /*!
     * \brief checks that the material can be processed by kernels (same
     * integration rule for all its elements, supported behaviour, geometric
     * data of its rule computed) and gets the data of the computations
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[out] data: data, whose rule is null if the material has no
     * element (nothing to compute)
     * \param[in] name: name of the calling method, for the messages
     * \param[in] g: geometric data of the integration rules
     */
    [[nodiscard]] bool getKernelsData(
        Context&,
        KernelsData&,
        const char* const,
        const batch::ElementsGeometrySet&) const noexcept;
    /*!
     * \return a view of the tangent operator blocks, in the global frame
     * (rotated after the integration of the behaviour for an orthotropic
     * behaviour), readable on the device if `d` is true (see
     * `MaterialArrays`)
     * \param[in] d: if true, the view is located on the device (if a device
     * is enabled)
     */
    [[nodiscard]] mgis::FieldView<const real> readTangentOperatorBlocks(
        const bool) const;
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
