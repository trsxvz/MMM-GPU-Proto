/*!
 * \file   include/MFEMMGIS/FBarBehaviourIntegratorCRTPBase.hxx
 * \brief  This file declares the `FBarBehaviourIntegratorCRTPBase` class
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#ifndef LIB_MFEM_MGIS_ISOTROPICFBARBEHAVIOURINTEGRATORCRTPBASE_HXX
#define LIB_MFEM_MGIS_ISOTROPICFBARBEHAVIOURINTEGRATORCRTPBASE_HXX

#include <mfem/linalg/densemat.hpp>
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
   */
  template <typename Child, Hypothesis H>
  struct FBarBehaviourIntegratorCRTPBase : BehaviourIntegratorBase {
    [[nodiscard]] bool supportsBatchIntegration() const noexcept override;
    [[nodiscard]] bool computeGradients(const mfem::FiniteElement&,
                                        mfem::ElementTransformation&,
                                        const mfem::Vector&) override;
    [[nodiscard]] bool integrateBatch(const IntegrationType) override;

   protected:
    static_assert(((H == Hypothesis::PLANESTRAIN) ||
                   (H == Hypothesis::TRIDIMENSIONAL)),
                  "invalid modelling hypothesis");
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
     * the behaviour
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
     * point (macroscopic gradients, contributions of the nodes, external
     * state variables computed from the unknowns, scaling by the
     * determinant at the center, rotation in the material frame)
     * \return the rotation matrix at the integration point
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
    auto computeGradientsAtIntegrationPoint(
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
