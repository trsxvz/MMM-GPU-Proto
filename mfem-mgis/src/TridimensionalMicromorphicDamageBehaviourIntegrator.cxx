/*!
 * \file   src/TridimensionalMicromorphicDamageBehaviourIntegrator.cxx
 * \brief  This file implements the
 * `TridimensionalMicromorphicDamageBehaviourIntegrator` class
 * \author Thomas Helfer
 * \date   31/03/2023
 */

#include <algorithm>
#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/TridimensionalMicromorphicDamageBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  TridimensionalMicromorphicDamageBehaviourIntegrator::selectIntegrationRule(
      const mfem::FiniteElement &e, const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }  // end of selectIntegrationRule

  PartialQuadratureSpace
  TridimensionalMicromorphicDamageBehaviourIntegrator::buildQuadratureSpace(
      attributes::Throwing,
      const FiniteElementDiscretization &fed,
      const size_type m) {
    auto selector = [&fed](const mfem::FiniteElement &e,
                           const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return applyIntegrationRuleOrder(fed, e, selectIntegrationRule(e, tr));
    };  // end of selector
    return {fed, LocationIdentifier{MaterialIdentifier{m}}, selector};
  }  // end of buildQuadratureSpace

  TridimensionalMicromorphicDamageBehaviourIntegrator::
      TridimensionalMicromorphicDamageBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            TridimensionalMicromorphicDamageBehaviourIntegrator>(
            buildQuadratureSpace(throwing, fed, m), std::move(b_ptr)) {
    this->checkBehaviourSymmetry(throwing, Behaviour::ISOTROPIC);
    this->checkHypothesis(throwing, Hypothesis::TRIDIMENSIONAL);
    this->checkBehaviourVariablesSizes(
        throwing, {.gradients_sizes = {1, 3},
                   .thermodynamic_forces_sizes = {1, 3},
                   .tangent_operator_blocks_sizes = {{1, 1}, {3, 3}}});
  }  // end of TridimensionalMicromorphicDamageBehaviourIntegrator

  real TridimensionalMicromorphicDamageBehaviourIntegrator::
      getIntegrationPointWeight(
          mfem::ElementTransformation &tr,
          const mfem::IntegrationPoint &ip) const noexcept {
    return ip.weight * tr.Weight();
  }  // end of getIntegrationPointWeight

  const mfem::IntegrationRule &
  TridimensionalMicromorphicDamageBehaviourIntegrator::getIntegrationRule(
      const mfem::FiniteElement &e,
      const mfem::ElementTransformation &t) const {
    return this->quadrature_space.getIntegrationRule(e, t);
  }  // end of getIntegrationRule

  TridimensionalMicromorphicDamageBehaviourIntegrator::RotationMatrix
  TridimensionalMicromorphicDamageBehaviourIntegrator::getRotationMatrix(
      const size_type) const {
    return RotationMatrix{};
  }  // end of getRotationMatrix

  void TridimensionalMicromorphicDamageBehaviourIntegrator::rotateGradients(
      std::span<real>, const RotationMatrix &) {}  // end of rotateGradients

  std::span<const real> TridimensionalMicromorphicDamageBehaviourIntegrator::
      rotateThermodynamicForces(std::span<const real> s,
                                const RotationMatrix &) {
    return s;
  }  // end of rotateThermodynamicForces

  void TridimensionalMicromorphicDamageBehaviourIntegrator::
      rotateTangentOperatorBlocks(std::span<real>, const RotationMatrix &) {
  }  // end of rotateTangentOperatorBlocks

  void TridimensionalMicromorphicDamageBehaviourIntegrator::updateGradients(
      std::span<real> &g,
      const mfem::Vector &u,
      const mfem::Vector &N,
      const mfem::DenseMatrix &dN,
      const size_type ni) noexcept {
    // micromorphic damage and the components of its gradient
    g[0] += u[ni] * N[ni];
    g[1] += u[ni] * dN(ni, 0);
    g[2] += u[ni] * dN(ni, 1);
    g[3] += u[ni] * dN(ni, 2);
  }  // end of updateGradients

  void TridimensionalMicromorphicDamageBehaviourIntegrator::updateInnerForces(
      mfem::Vector &Fe,
      const std::span<const real> &s,
      const mfem::Vector &N,
      const mfem::DenseMatrix &dN,
      const real w,
      const size_type ni) const noexcept {
    // s contains the dual force conjugated with the micromorphic damage
    // (scalar) and the dual force conjugated with its gradient (vector)
    Fe[ni] += w * (s[0] * N[ni] +      //
                   s[1] * dN(ni, 0) +  //
                   s[2] * dN(ni, 1) +  //
                   s[3] * dN(ni, 2));
  }  // end of updateInnerForces

  void
  TridimensionalMicromorphicDamageBehaviourIntegrator::updateStiffnessMatrix(
      mfem::DenseMatrix &Ke,
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
    const auto dNi_2 = dN(ni, 2);
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto Nj = N[nj];
      const auto dNj_0 = dN(nj, 0);
      const auto dNj_1 = dN(nj, 1);
      const auto dNj_2 = dN(nj, 2);
      Ke(ni, nj) += w * (Kip[0] * Ni * Nj +  //
                         dNi_1 * dNj_0 * Kip[4] + dNi_1 * dNj_1 * Kip[5] +
                         Kip[3] * dNj_2 * dNi_0 + dNi_2 * dNj_1 * Kip[8] +
                         dNj_0 * Kip[1] * dNi_0 + Kip[2] * dNi_0 * dNj_1 +
                         dNi_1 * Kip[6] * dNj_2 + Kip[9] * dNj_2 * dNi_2 +
                         dNi_2 * dNj_0 * Kip[7]);
    }
  }  // end of updateStiffnessMatrix

  bool TridimensionalMicromorphicDamageBehaviourIntegrator::integrate(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u,
      const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void TridimensionalMicromorphicDamageBehaviourIntegrator::updateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void TridimensionalMicromorphicDamageBehaviourIntegrator::updateJacobian(
      mfem::DenseMatrix &Ke,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &) {
    this->implementUpdateJacobian(Ke, e, tr);
  }  // end of updateJacobian

  void TridimensionalMicromorphicDamageBehaviourIntegrator::computeInnerForces(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  TridimensionalMicromorphicDamageBehaviourIntegrator::
      ~TridimensionalMicromorphicDamageBehaviourIntegrator() = default;

}  // end of namespace mfem_mgis
