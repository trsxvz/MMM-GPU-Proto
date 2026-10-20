/*!
 * \file   src/OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator.cxx
 * \brief  This file implements the
 * `OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator` class
 * \author Thomas Helfer
 * \date   07/12/2021
 */

#include <algorithm>
#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      selectIntegrationRule(const mfem::FiniteElement &e,
                            const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }  // end of selectIntegrationRule

  PartialQuadratureSpace
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      buildQuadratureSpace(attributes::Throwing,
                           const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [&fed](const mfem::FiniteElement &e,
                           const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return applyIntegrationRuleOrder(fed, e, selectIntegrationRule(e, tr));
    };  // end of selector
    return {fed, LocationIdentifier{MaterialIdentifier{m}}, selector};
  }  // end of buildQuadratureSpace

  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator>(
            buildQuadratureSpace(throwing, fed, m), std::move(b_ptr)) {
    this->checkBehaviourSymmetry(throwing, Behaviour::ORTHOTROPIC);
    this->checkBehaviourVariablesSizes(
        throwing, {.gradients_sizes = {1, 2},
                   .thermodynamic_forces_sizes = {1, 2},
                   .tangent_operator_blocks_sizes = {{1, 1}, {2, 2}}});
  }  // end of OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator

  real OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      getIntegrationPointWeight(
          mfem::ElementTransformation &tr,
          const mfem::IntegrationPoint &ip) const noexcept {
    return ip.weight * tr.Weight();
  }  // end of getIntegrationPointWeight

  const mfem::IntegrationRule &
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &e,
                         const mfem::ElementTransformation &t) const {
    return this->quadrature_space.getIntegrationRule(e, t);
  }  // end of getIntegrationRule

  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::RotationMatrix
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      getRotationMatrix(const size_type i) const {
    return this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
  }  // end of getRotationMatrix

  void OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      rotateGradients(std::span<real> g, const RotationMatrix &r) {
    this->b.rotate_gradients_ptr(g.data(), g.data(), r.data());
  }  // end of rotateGradients

  std::array<real, 3>
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      rotateThermodynamicForces(std::span<const real> s,
                                const RotationMatrix &r) {
    std::array<real, 3> rs;
    std::copy(s.begin(), s.end(), rs.begin());
    this->b.rotate_thermodynamic_forces_ptr(rs.data(), rs.data(), r.data());
    return rs;
  }  // end of rotateThermodynamicForces

  void OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      rotateTangentOperatorBlocks(std::span<real> Kip,
                                  const RotationMatrix &r) {
    this->b.rotate_tangent_operator_blocks_ptr(Kip.data(), Kip.data(),
                                               r.data());
  }  // end of rotateTangentOperatorBlocks

  void OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      updateGradients(std::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::Vector &N,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    // micromorphic damage and the components of its gradient
    g[0] += u[ni] * N[ni];
    g[1] += u[ni] * dN(ni, 0);
    g[2] += u[ni] * dN(ni, 1);
  }  // end of updateGradients

  void OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const std::span<const real> &s,
                        const mfem::Vector &N,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    // s contains the dual force conjugated with the micromorphic damage
    // (scalar) and the dual force conjugated with its gradient (vector)
    Fe[ni] += w * (s[0] * N[ni] +      //
                   s[1] * dN(ni, 0) +  //
                   s[2] * dN(ni, 1));
  }  // end of updateInnerForces

  void OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const std::span<const real> &Kip,
                            const mfem::Vector &N,
                            const mfem::DenseMatrix &dN,
                            const real w,
                            const size_type ni) const noexcept {
    // Kip contains the derivative of the first dual force with respect to
    // the micromorphic damage (scalar) and the derivative of the second dual
    // force with respect to the gradient of the micromorphic damage (matrix)
    const auto nnodes = N.Size();
    const auto Ni = N[ni];
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto Nj = N[nj];
      const auto dNj_0 = dN(nj, 0);
      const auto dNj_1 = dN(nj, 1);
      Ke(ni, nj) += w * (Kip[0] * Ni * Nj +        //
                         Kip[2] * dNi_0 * dNj_1 +  //
                         dNj_0 * Kip[3] * dNi_1 +  //
                         dNj_1 * dNi_1 * Kip[4] +  //
                         dNj_0 * dNi_0 * Kip[1]);
    }
  }  // end of updateStiffnessMatrix

  bool OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::integrate(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u,
      const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::updateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::updateJacobian(
      mfem::DenseMatrix &Ke,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &) {
    this->implementUpdateJacobian(Ke, e, tr);
  }  // end of updateJacobian

  void OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      ~OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
