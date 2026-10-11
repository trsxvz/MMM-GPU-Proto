/*!
 * \file MFEMMGIS/Faltus2026RegularizedBehaviourIntegrators.hxx
 * \brief  This file declares the
 * `Faltus2026RegularizedIsotropicBehaviourIntegrator` class and the functions
 * generating it
 * \author Thomas Helfer
 * \date   17/03/2026
 */

#ifndef LIB_MFEM_MGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_HXX
#define LIB_MFEM_MGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"
#include "MFEMMGIS/Batch/Faltus.hxx"
#include "IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \brief select the base class of the regularized behaviour integrator
   * \tparam H: modelling hypothesis
   */
  template <Hypothesis H>
  struct Faltus2026RegularizedIsotropicBehaviourIntegratorBaseDispatch;

  //! \brief plane strain specialisation
  template <>
  struct Faltus2026RegularizedIsotropicBehaviourIntegratorBaseDispatch<
      Hypothesis::PLANESTRAIN> {
    //! \brief base class of the regularized behaviour integrator
    using type =
        IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator;
  };

  //! \brief plane stress specialisation
  template <>
  struct Faltus2026RegularizedIsotropicBehaviourIntegratorBaseDispatch<
      Hypothesis::PLANESTRESS> {
    //! \brief base class of the regularized behaviour integrator
    using type =
        IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator;
  };

  //! \brief tridimensional specialisation
  template <>
  struct Faltus2026RegularizedIsotropicBehaviourIntegratorBaseDispatch<
      Hypothesis::TRIDIMENSIONAL> {
    //! \brief base class of the regularized behaviour integrator
    using type =
        IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator;
  };

  //! \brief base class of the regularized behaviour integrator
  template <Hypothesis H>
  using Faltus2026RegularizedIsotropicBehaviourIntegratorBase =
      typename Faltus2026RegularizedIsotropicBehaviourIntegratorBaseDispatch<
          H>::type;

  /*!
   * \brief a behaviour integrator which enhances a standard finite strain
   * behaviour with the regularization proposed by Faltus et al.
   *
   * This regularization adds a term penalizing the difference between the
   * deformation gradient \f$\underline{F}\f$ and its value
   * \f$\bar{\underline{F}}\f$ at the centroid of the element:
   *
   * \f[
   * W\left(\underline{F}, \bar{\underline{F}}\right) =
   * \frac{\alpha}{2}\,\left(\underline{F}-\bar{\underline{F}}\right)\,\colon\,
   * \left(\underline{F}-\bar{\underline{F}}\right)
   * \f]
   */
  template <Hypothesis H>
  struct Faltus2026RegularizedIsotropicBehaviourIntegrator final
      : Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H> {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     * \param[in] params: parameters defining the penalization coefficient
     */
    Faltus2026RegularizedIsotropicBehaviourIntegrator(
        const FiniteElementDiscretization& fed,
        const size_type m,
        std::unique_ptr<const Behaviour> b_ptr,
        const Parameters& params);
    /*!
     * \brief compute the contribution of the given element to the residual,
     * including the regularization term
     * \param[out] Fe: element contribution to the residual
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     */
    void updateResidual(mfem::Vector& Fe,
                        const mfem::FiniteElement& e,
                        mfem::ElementTransformation& tr,
                        const mfem::Vector& u) override;
    /*!
     * \brief compute the contribution of the given element to the jacobian,
     * including the regularization term
     * \param[out] Je: element stiffness matrix
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     */
    void updateJacobian(mfem::DenseMatrix& Je,
                        const mfem::FiniteElement& e,
                        mfem::ElementTransformation& tr,
                        const mfem::Vector& u) override;
    //! \return if the current solution is required for assembling the residual
    [[nodiscard]] bool requiresCurrentSolutionForResidualAssembly()
        const noexcept override;
    //! \return if the current solution is required for assembling the jacobian
    [[nodiscard]] bool requiresCurrentSolutionForJacobianAssembly()
        const noexcept override;
    /*!
     * \brief compute the gradients of an element: the deformation gradients
     * at the center of the elements computed by kernels no longer match the
     * gradients (see `batch::FaltusData::valid`)
     */
    [[nodiscard]] bool computeGradients(const mfem::FiniteElement&,
                                        mfem::ElementTransformation&,
                                        const mfem::Vector&) override;
    /*!
     * \brief integrate the behaviour on an element: the deformation gradients
     * at the center of the elements computed by kernels no longer match the
     * gradients (see `batch::FaltusData::valid`)
     */
    [[nodiscard]] bool integrate(const mfem::FiniteElement&,
                                 mfem::ElementTransformation&,
                                 const mfem::Vector&,
                                 const IntegrationType) override;
    void revert() override;
    void update() override;
    /*!
     * \return true in 3D and in plane strain (deformation gradient of the
     * kind of the hypothesis), if the kernels are available (see
     * `batch::areGradientsKernelsAvailable`) and if all the elements of the
     * material use the same integration rule: the regularization is
     * computed by the kernels of `MFEMMGIS/Batch/Faltus.hxx`
     */
    [[nodiscard]] bool supportsKernels() const noexcept override;
    //! \return the integration rule of the center of the elements
    [[nodiscard]] std::vector<const mfem::IntegrationRule*>
    getAdditionalIntegrationRulesForKernels() const override;
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

   private:
    /*!
     * \brief gets the data of the computations by kernels (see
     * `getKernelsData`) and the geometric data of the center of the elements
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[out] data: data of the computations
     * \param[out] g0: geometric data of the center of the elements, null if
     * the material has no element
     * \param[in] name: name of the calling method
     * \param[in] g: geometric data of the integration rules
     */
    [[nodiscard]] bool getFaltusKernelsData(
        Context&,
        BehaviourIntegratorBase::KernelsData&,
        const batch::ElementsGeometry*&,
        const char* const,
        const batch::ElementsGeometrySet&) const noexcept;
    //! \brief penalization coefficient
    const real alpha;
    //! \brief data of the regularization computed by kernels
    batch::FaltusData faltus_data;
#ifndef MFEM_THREAD_SAFE
    /*!
     * \brief matrix used to store the derivatives of the shape functions at
     * the center of the element
     */
    mfem::DenseMatrix dshape0;
#endif
  };

  /*!
   * \brief generate a behaviour integrator with the regularization proposed
   * by Faltus et al. under the plane strain hypothesis
   * \param[in, out] ctx: execution context
   * \param[in] fed: finite element discretization
   * \param[in] m: material attribute
   * \param[in] b: behaviour
   * \param[in] params: parameters defining the penalization coefficient
   * \return the new behaviour integrator, a null pointer on failure
   */
  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStrainFaltus2026RegularizedMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept;

  /*!
   * \brief generate a behaviour integrator with the regularization proposed
   * by Faltus et al. under the plane stress hypothesis
   * \param[in, out] ctx: execution context
   * \param[in] fed: finite element discretization
   * \param[in] m: material attribute
   * \param[in] b: behaviour
   * \param[in] params: parameters defining the penalization coefficient
   * \return the new behaviour integrator, a null pointer on failure
   */
  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStressFaltus2026RegularizedMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept;

  /*!
   * \brief generate a behaviour integrator with the regularization proposed
   * by Faltus et al. under the tridimensional hypothesis
   * \param[in, out] ctx: execution context
   * \param[in] fed: finite element discretization
   * \param[in] m: material attribute
   * \param[in] b: behaviour
   * \param[in] params: parameters defining the penalization coefficient
   * \return the new behaviour integrator, a null pointer on failure
   */
  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generateTridimensionalFaltus2026RegularizedMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Faltus2026RegularizedBehaviourIntegrators.ixx"

#endif /* LIB_MFEM_MGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_HXX */
