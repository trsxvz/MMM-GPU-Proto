/*!
 * \file include/MFEMMGIS/BidimensionalMicromorphicDamageBehaviourIntegrator.hxx
 * \brief  This file declares the
 * `BidimensionalMicromorphicDamageBehaviourIntegrator` class
 * \date   07/12/2021
 */

#ifndef LIB_MFEM_MGIS_BIDIMENSIONALMICROMORPHICDAMAGEBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_BIDIMENSIONALMICROMORPHICDAMAGEBEHAVIOURINTEGRATOR_HXX

#include <span>
#include <array>
#include "mfem/linalg/densemat.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"
#include "MFEMMGIS/StandardBehaviourIntegratorCRTPBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;
  // forward declaration
  struct BidimensionalMicromorphicDamageBehaviourIntegrator;

  /*!
   * \brief specialisation of the `BehaviourIntegratorTraits` class for the
   * `BidimensionalMicromorphicDamageBehaviourIntegrator` behaviour integrator
   */
  template <>
  struct BehaviourIntegratorTraits<
      BidimensionalMicromorphicDamageBehaviourIntegrator> {
    //! \brief number of components of the unknowns
    static constexpr size_type unknownsSize = 1;
    //! \brief if the computation of the gradients requires the shape functions
    static constexpr bool gradientsComputationRequiresShapeFunctions = true;
    /*!
     * \brief if the computation of the gradients requires the derivatives of
     * the shape functions
     */
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = true;
    //! \brief if the external state variables are updated from the unknowns
    static constexpr bool updateExternalStateVariablesFromUnknownsValues =
        false;
  };  // end of struct BehaviourIntegratorTraits<>

  /*!
   * \brief behaviour integrator dedicated to micromorphic damage in two
   dimensions (the modelling hypothesis has no effect on this specific
   behaviour as out of plane damage gradients are assumed to be zero). The
   gradients are the micromorphic damage and its gradient, the thermodynamic
   forces their dual forces, and the tangent operator blocks the derivative
   of the first dual force with respect to the micromorphic damage and the
   derivative of the second dual force with respect to the gradient.
   */
  struct MFEM_MGIS_EXPORT
      BidimensionalMicromorphicDamageBehaviourIntegrator final
      : StandardBehaviourIntegratorCRTPBase<
            BidimensionalMicromorphicDamageBehaviourIntegrator> {
    //! \brief a dummy structure (isotropic behaviour)
    struct RotationMatrix {};
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] m: material attribute
     * \param[in] b_ptr: behaviour
     */
    BidimensionalMicromorphicDamageBehaviourIntegrator(
        const FiniteElementDiscretization &fed,
        const size_type m,
        std::unique_ptr<const Behaviour> b_ptr);
    /*!
     * \return the rotation matrix associated with the given integration
     * point
     * \param[in] i: offset of the integration point
     */
    inline RotationMatrix getRotationMatrix(const size_type i) const;
    /*!
     * \brief rotate the gradients in the material frame
     * \param[in, out] g: gradients
     * \param[in] r: rotation matrix
     */
    inline void rotateGradients(std::span<real> g, const RotationMatrix &r);
    /*!
     * \return the thermodynamic forces in the global frame
     * \param[in] s: thermodynamic forces
     * \param[in] r: rotation matrix
     */
    inline std::span<const real> rotateThermodynamicForces(
        std::span<const real> s, const RotationMatrix &r);
    /*!
     * \brief rotate the tangent operator blocks in the global frame
     * \param[in, out] Kip: tangent operator blocks
     * \param[in] r: rotation matrix
     */
    inline void rotateTangentOperatorBlocks(std::span<real> Kip,
                                            const RotationMatrix &r);
    const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &e,
        const mfem::ElementTransformation &tr) const override;
    real getIntegrationPointWeight(
        mfem::ElementTransformation &tr,
        const mfem::IntegrationPoint &ip) const noexcept override;
    bool integrate(const mfem::FiniteElement &e,
                   mfem::ElementTransformation &tr,
                   const mfem::Vector &u,
                   const IntegrationType it) override;
    void updateResidual(mfem::Vector &Fe,
                        const mfem::FiniteElement &e,
                        mfem::ElementTransformation &tr,
                        const mfem::Vector &u) override;
    void updateJacobian(mfem::DenseMatrix &Ke,
                        const mfem::FiniteElement &e,
                        mfem::ElementTransformation &tr,
                        const mfem::Vector &u) override;
    void computeInnerForces(mfem::Vector &Fe,
                            const mfem::FiniteElement &e,
                            mfem::ElementTransformation &tr) override;
    //! \brief destructor
    ~BidimensionalMicromorphicDamageBehaviourIntegrator() override;

   protected:
    //! \brief allow the CRTP base class to access the protected members
    friend struct StandardBehaviourIntegratorCRTPBase<
        BidimensionalMicromorphicDamageBehaviourIntegrator>;
    /*!
     * \brief select the integration rule for the given element and element
     * transformation
     * \param[in] e: element
     * \param[in] t: element transformation
     * \return the integration rule
     */
    static const mfem::IntegrationRule &selectIntegrationRule(
        const mfem::FiniteElement &e, const mfem::ElementTransformation &t);
    /*!
     * \brief build the quadrature space for the given material
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \return the partial quadrature space
     */
    static std::shared_ptr<const PartialQuadratureSpace> buildQuadratureSpace(
        const FiniteElementDiscretization &fed, const size_type m);
    /*!
     * \brief update the gradients at the integration point with the
     * contribution of the given node: micromorphic damage and its gradient
     * \param[in, out] g: gradients
     * \param[in] u: nodal values of the micromorphic damage
     * \param[in] N: shape functions
     * \param[in] dN: derivatives of the shape functions
     * \param[in] ni: node index
     */
    void updateGradients(std::span<real> &g,
                         const mfem::Vector &u,
                         const mfem::Vector &N,
                         const mfem::DenseMatrix &dN,
                         const size_type ni) noexcept;
    /*!
     * \brief update the inner forces of the given node with the contribution
     * of the thermodynamic forces of an integration point
     * \param[in, out] Fe: inner forces
     * \param[in] s: thermodynamic forces
     * \param[in] N: shape functions
     * \param[in] dN: derivatives of the shape functions
     * \param[in] w: weight of the integration point
     * \param[in] ni: node index
     */
    void updateInnerForces(mfem::Vector &Fe,
                           const std::span<const real> &s,
                           const mfem::Vector &N,
                           const mfem::DenseMatrix &dN,
                           const real w,
                           const size_type ni) const noexcept;
    /*!
     * \brief update the stiffness matrix of the given node with the
     * contribution of the tangent operator blocks of an integration point
     * \param[in, out] Ke: stiffness matrix
     * \param[in] Kip: tangent operator blocks of the integration point
     * \param[in] N: shape functions
     * \param[in] dN: derivatives of the shape functions
     * \param[in] w: weight of the integration point
     * \param[in] ni: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                               const std::span<const real> &Kip,
                               const mfem::Vector &N,
                               const mfem::DenseMatrix &dN,
                               const real w,
                               const size_type ni) const noexcept;
  };  // end of struct BidimensionalMicromorphicDamageBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_BIDIMENSIONALMICROMORPHICDAMAGEBEHAVIOURINTEGRATOR_HXX \
        */
