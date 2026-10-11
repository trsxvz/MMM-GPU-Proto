/*!
 * \file   include/MFEMMGIS/FBarBehaviourIntegratorCRTPBase.ixx
 * \brief  This file implements the inline methods of the
 * `FBarBehaviourIntegratorCRTPBase` class
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#ifndef LIB_MFEM_MGIS_FBARBEHAVIOURINTEGRATORCRTPBASE_IXX
#define LIB_MFEM_MGIS_FBARBEHAVIOURINTEGRATORCRTPBASE_IXX

#include <string>
#include <exception>
#include "mfem/general/device.hpp"
#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MGIS/Raise.hxx"
#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/t2tot2.hxx"
#include "TFEL/Math/Array/View.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"

namespace mfem_mgis {

  template <typename Child, Hypothesis H>
  FBarBehaviourIntegratorCRTPBase<Child, H>::FBarBehaviourIntegratorCRTPBase(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::unique_ptr<const Behaviour> b_ptr)
      : BehaviourIntegratorBase(s, std::move(b_ptr)) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr auto usize = static_cast<size_type>(Traits::unknownsSize);
    const auto &fed = s->getFiniteElementDiscretization();
    const auto nc = getNumberOfComponents(fed);
    if (nc != usize) {
      raise("the number of components of the finite element space (" +
            std::to_string(nc) +
            ") does not match "
            "the number of components expected by the behaviour integrator (" +
            std::to_string(usize) + ")");
    }
  }  // end of FBarBehaviourIntegratorCRTPBase

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::implementIntegrate(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u,
      const IntegrationType it) {
    //
    using Traits = BehaviourIntegratorTraits<Child>;
    static_assert(
        ((Traits::gradientsComputationRequiresShapeFunctions) ||
         (Traits::gradientsComputationRequiresShapeFunctionsDerivatives)),
        "neither the shape functions or their derivatives are required to "
        "compute the gradients, this is not supported");
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    static_assert(evaluateShapeFunctionsDerivatives,
                  "the derivatives of the shape functions are required to "
                  "compute the gradients");
    //
    auto &child = static_cast<Child &>(*this);
    // the tangent operator blocks change
    this->fbar_data.has_tangent_operator = false;
    // element offset
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    const auto &ir = child.getIntegrationRule(e, tr);
    if ((it == IntegrationType::PREDICTION_TANGENT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_SECANT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_ELASTIC_OPERATOR)) {
      for (size_type i = 0; i != ir.GetNPoints(); ++i) {
        // offset of the integration point
        const auto o = eoffset + i;
        if (!this->performsLocalBehaviourIntegration(o, it)) {
          return false;
        }
        // rotate the tangent operator blocks
        const auto r = child.getRotationMatrix(o);
        auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
        child.rotateTangentOperatorBlocks(Kip, r);
      }
      return true;
    }
    // the gradients are computed element by element: the data of the FBar
    // formulation computed by kernels no longer match the unknowns
    this->fbar_data.valid = false;
    this->fbar_data.has_tangent_operator = false;
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape0;
    mfem::DenseMatrix dshape0;
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
    if constexpr (evaluateShapeFunctions) {
      shape0.SetSize(e.GetDof());
      shape.SetSize(e.GetDof());
    }
    dshape0.SetSize(e.GetDof(), e.GetDim());
    dshape.SetSize(e.GetDof(), e.GetDim());
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape0.SetSize(e.GetDof());
      this->shape.SetSize(e.GetDof());
    }
    this->dshape0.SetSize(e.GetDof(), e.GetDim());
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto J0 = this->computeDeformationGradientDeterminantAtCenter(
        e, tr, u, shape0, dshape0);
    const auto gsize = this->s1.gradients_stride;
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      // offset of the integration point
      const auto o = eoffset + i;
      this->computeGradientsAtIntegrationPoint(e, tr, u, ir.IntPoint(i), o, J0,
                                               shape, dshape);
      const auto r = child.getRotationMatrix(o);
      child.rotateGradients(this->s1.gradients.subspan(o * gsize, gsize), r);
      if (!this->performsLocalBehaviourIntegration(o, it)) {
        return false;
      }
      // Here we rotate the tangent operator blocks but not the thermodynamic
      // forces.
      if (it != IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR) {
        auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
        child.rotateTangentOperatorBlocks(Kip, r);
      }
    }
    return true;
  }  // end of implementIntegrate

  template <typename Child, Hypothesis H>
  real FBarBehaviourIntegratorCRTPBase<Child, H>::
      computeDeformationGradientDeterminantAtCenter(
          const mfem::FiniteElement &e,
          mfem::ElementTransformation &tr,
          const mfem::Vector &u,
          mfem::Vector &shape_functions,
          mfem::DenseMatrix &shape_functions_derivatives) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr auto dimension = (H == Hypothesis::PLANESTRAIN) ? 2 : 3;
    auto &child = static_cast<Child &>(*this);
    const auto nnodes = e.GetDof();
    auto F0 = tfel::math::tensor<dimension, real>{};
    std::copy(this->macroscopic_gradients.begin(),
              this->macroscopic_gradients.end(), F0.begin());
    // computation of the mean value of the deformation gradient
    auto F0v = std::span<real>(F0.data(), F0.size());
    const auto &ir0 = mfem::IntRules.Get(e.GetGeomType(), 0);
    const auto &ip = ir0.IntPoint(0);
    tr.SetIntPoint(&ip);
    if constexpr (evaluateShapeFunctions) {
      e.CalcPhysShape(tr, shape_functions);
    }
    e.CalcPhysDShape(tr, shape_functions_derivatives);
    for (size_type ni = 0; ni != nnodes; ++ni) {
      if constexpr (Traits::gradientsComputationRequiresShapeFunctions) {
        child.updateGradients(F0v, u, shape_functions,
                              shape_functions_derivatives, ni);
      } else {
        child.updateGradients(F0v, u, shape_functions_derivatives, ni);
      }
    }
    return tfel::math::det(F0);
  }  // end of computeDeformationGradientDeterminantAtCenter

  template <typename Child, Hypothesis H>
  void FBarBehaviourIntegratorCRTPBase<Child, H>::
      computeGradientsAtIntegrationPoint(
          const mfem::FiniteElement &e,
          mfem::ElementTransformation &tr,
          const mfem::Vector &u,
          const mfem::IntegrationPoint &ip,
          const size_type o,
          const real J0,
          mfem::Vector &shape_functions,
          mfem::DenseMatrix &shape_functions_derivatives) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr auto dimension = (H == Hypothesis::PLANESTRAIN) ? 2 : 3;
    // inverse of the exponent of the FBar regularisation
    constexpr auto ie = (H == Hypothesis::PLANESTRAIN) ? 2 : 3;
    constexpr auto gsize = tfel::math::TensorDimeToSize<dimension>::value;
    auto &child = static_cast<Child &>(*this);
    const auto nnodes = e.GetDof();
    tr.SetIntPoint(&ip);
    if constexpr (evaluateShapeFunctions) {
      // get the gradients of the shape functions
      e.CalcPhysShape(tr, shape_functions);
    }
    if constexpr (Traits::updateExternalStateVariablesFromUnknownsValues) {
      child.updateExternalStateVariablesFromUnknownsValues(u, shape_functions,
                                                           o);
    }
    // get the gradients of the shape functions
    e.CalcPhysDShape(tr, shape_functions_derivatives);
    auto F = tfel::math::tensor<dimension, real>{};
    std::copy(this->macroscopic_gradients.begin(),
              this->macroscopic_gradients.end(), F.begin());
    auto Fv = std::span<real>(F.data(), F.size());
    for (size_type ni = 0; ni != nnodes; ++ni) {
      if constexpr (Traits::gradientsComputationRequiresShapeFunctions) {
        child.updateGradients(Fv, u, shape_functions,
                              shape_functions_derivatives, ni);
      } else {
        child.updateGradients(Fv, u, shape_functions_derivatives, ni);
      }
    }
    const auto J = tfel::math::det(F);
    const auto alpha = tfel::math::power<1, ie>(J0 / J);
    auto g = this->s1.gradients.subspan(o * gsize, gsize);
    if constexpr (H == Hypothesis::PLANESTRAIN) {
      for (size_type c = 0; c != gsize; ++c) {
        g[c] = alpha * F[c];
      }
      g[2] = 1;
    } else {
      for (size_type c = 0; c != gsize; ++c) {
        g[c] = alpha * F[c];
      }
    }
  }  // end of computeGradientsAtIntegrationPoint

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::implementComputeGradients(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    auto &child = static_cast<Child &>(*this);
    // the data of the FBar formulation computed by kernels no longer match
    // the unknowns
    this->fbar_data.valid = false;
    this->fbar_data.has_tangent_operator = false;
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    const auto &ir = child.getIntegrationRule(e, tr);
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape0;
    mfem::DenseMatrix dshape0;
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
    if constexpr (evaluateShapeFunctions) {
      shape0.SetSize(e.GetDof());
      shape.SetSize(e.GetDof());
    }
    dshape0.SetSize(e.GetDof(), e.GetDim());
    dshape.SetSize(e.GetDof(), e.GetDim());
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape0.SetSize(e.GetDof());
      this->shape.SetSize(e.GetDof());
    }
    this->dshape0.SetSize(e.GetDof(), e.GetDim());
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto J0 = this->computeDeformationGradientDeterminantAtCenter(
        e, tr, u, shape0, dshape0);
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      this->computeGradientsAtIntegrationPoint(e, tr, u, ir.IntPoint(i),
                                               eoffset + i, J0, shape, dshape);
    }
    return true;
  }  // end of implementComputeGradients

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::supportsBatchIntegration()
      const noexcept {
    return true;
  }  // end of supportsBatchIntegration

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::computeGradients(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    return this->implementComputeGradients(e, tr, u);
  }  // end of computeGradients

  template <typename Child, Hypothesis H>
  void FBarBehaviourIntegratorCRTPBase<Child, H>::implementUpdateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of implementUpdateResidual

  template <typename Child, Hypothesis H>
  void FBarBehaviourIntegratorCRTPBase<Child, H>::implementComputeInnerForces(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    static_assert(evaluateShapeFunctionsDerivatives,
                  "the derivatives of the shape functions are required to "
                  "compute the inner forces");
    auto &child = static_cast<Child &>(*this);
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
    if constexpr (evaluateShapeFunctions) {
      shape.SetSize(e.GetDof());
    }
    dshape.SetSize(e.GetDof(), e.GetDim());
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape.SetSize(e.GetDof());
    }
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto nnodes = e.GetDof();
    const auto thsize = this->s1.thermodynamic_forces_stride;
    // thermodynamic forces in the global frame: rotated point by point, or,
    // for an orthotropic behaviour without functions rotating the values of
    // one point (generic-parallel interface), rotated at all the points at
    // once by the integration by batches
    const auto rotated_thforces = this->getRotatedThermodynamicForces();
    // element offset
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Fe.SetSize(e.GetDof() * Traits::unknownsSize);
    Fe = 0.;
    const auto &ir = child.getIntegrationRule(e, tr);
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      if constexpr (evaluateShapeFunctions) {
        // get the shape functions
        e.CalcPhysShape(tr, shape);
      }
      if constexpr (evaluateShapeFunctionsDerivatives) {
        // get the gradients of the shape functions
        e.CalcPhysDShape(tr, dshape);
      }
      // get the weights associated to point ip
      const auto w = child.getIntegrationPointWeight(tr, ip);
      // offset of the integration point
      const auto o = eoffset + i;
      auto update = [&](const std::span<const real>& rs) {
        for (size_type ni = 0; ni != nnodes; ++ni) {
          if constexpr (evaluateShapeFunctionsDerivatives) {
            child.updateInnerForces(Fe, rs, dshape, w, ni);
          } else {
            child.updateInnerForces(Fe, rs, shape, w, ni);
          }
        }
      };
      if (!rotated_thforces.empty()) {
        update(rotated_thforces.subspan(o * thsize, thsize));
      } else {
        const auto r = child.getRotationMatrix(o);
        const auto s =
            this->s1.thermodynamic_forces.subspan(o * thsize, thsize);
        update(child.rotateThermodynamicForces(s, r));
      }
    }
  }  // end of implementComputeInnerForces

  template <typename Child, Hypothesis H>
  void FBarBehaviourIntegratorCRTPBase<Child, H>::implementUpdateJacobian(
      mfem::DenseMatrix &Ke,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    using Traits = BehaviourIntegratorTraits<Child>;
    auto &child = static_cast<Child &>(*this);
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    constexpr auto dimension = []() constexpr->unsigned short {
      return (H == Hypothesis::PLANESTRAIN) ? 2 : 3;
    }
    ();
    // inverse of the exponent of the FBar regularisation
    constexpr auto ie = []() constexpr->unsigned short {
      return (H == Hypothesis::PLANESTRAIN) ? 2 : 3;
    }
    ();
    constexpr auto gsize = tfel::math::TensorDimeToSize<dimension>::value;
    static_assert(
        ((evaluateShapeFunctions) || (evaluateShapeFunctionsDerivatives)),
        "neither the shape functions or their derivatives are required to "
        "compute the jacobian, this is not supported");
    static_assert(evaluateShapeFunctionsDerivatives,
                  "the derivatives of the shape functions are required to "
                  "compute the jacobian, this is not supported");
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape0, shape;
    mfem::DenseMatrix dshape0, dshape;
    if constexpr (evaluateShapeFunctions) {
      shape0.SetSize(e.GetDof());
      shape.SetSize(e.GetDof());
    }
    dshape0.SetSize(e.GetDof(), e.GetDim());
    dshape.SetSize(e.GetDof(), e.GetDim());
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape0.SetSize(e.GetDof());
      this->shape.SetSize(e.GetDof());
    }
    this->dshape0.SetSize(e.GetDof(), e.GetDim());
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto nnodes = e.GetDof();
    auto F0 = tfel::math::tensor<dimension, real>{};
    std::copy(this->macroscopic_gradients.begin(),
              this->macroscopic_gradients.end(), F0.begin());
    {  // computation of the mean value of the deformation gradient
      auto F0v = std::span<real>(F0.data(), F0.size());
      const auto &ir = mfem::IntRules.Get(e.GetGeomType(), 0);
      const auto &ip = ir.IntPoint(0);
      tr.SetIntPoint(&ip);
      if constexpr (evaluateShapeFunctions) {
        e.CalcPhysShape(tr, shape0);
      }
      e.CalcPhysDShape(tr, dshape0);
      for (size_type ni = 0; ni != nnodes; ++ni) {
        if constexpr (Traits::gradientsComputationRequiresShapeFunctions) {
          child.updateGradients(F0v, u, shape0, dshape0, ni);
        } else {
          child.updateGradients(F0v, u, dshape0, ni);
        }
      }
    }
    const auto J0 = tfel::math::det(F0);
    // element offset
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Ke.SetSize(e.GetDof() * Traits::unknownsSize,
               e.GetDof() * Traits::unknownsSize);
    Ke = 0.;
    const auto &ir = child.getIntegrationRule(e, tr);
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      // get the gradients of the shape functions
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      if constexpr (evaluateShapeFunctions) {
        // get the shape functions
        e.CalcPhysShape(tr, shape);
      }
      // get the derivatives of the shape functions
      e.CalcPhysDShape(tr, dshape);
      // get the weights associated to point ip
      const auto w = child.getIntegrationPointWeight(tr, ip);
      // offset of the integration point
      const auto o = eoffset + i;
      auto F = tfel::math::tensor<dimension, real>{};
      std::copy(this->macroscopic_gradients.begin(),
                this->macroscopic_gradients.end(), F.begin());
      auto Fv = std::span<real>(F.data(), F.size());
      for (size_type ni = 0; ni != nnodes; ++ni) {
        if constexpr (Traits::gradientsComputationRequiresShapeFunctions) {
          child.updateGradients(Fv, u, shape, dshape, ni);
        } else {
          child.updateGradients(Fv, u, dshape, ni);
        }
      }
      const auto J = tfel::math::det(F);
      const auto alpha = tfel::math::power<1, ie>(J0 / J);
      auto dJ_dF = tfel::math::computeDeterminantDerivative(F);
      auto dJ_dF0 = tfel::math::computeDeterminantDerivative(F0);
      const auto dalpha_dF0 = (alpha / (ie * J0)) * dJ_dF0;
      const auto dalpha_dF = (-alpha / (ie * J)) * dJ_dF;
      auto dFbar_dF0 = tfel::math::eval(F ^ dalpha_dF0);
      auto dFbar_dF = tfel::math::eval(
          alpha * tfel::math::t2tot2<dimension, real>::Id() + (F ^ dalpha_dF));
      if constexpr (ie == 2u) {
        // those are brutal corrections that seems working
        // we shall verify this someday
        for (size_type idx = 0; idx != gsize; ++idx) {
          dFbar_dF0(2, idx) = 0;
          dFbar_dF0(idx, 2) = 0;
          dFbar_dF(2, idx) = 0;
          dFbar_dF(idx, 2) = 0;
        }
        dFbar_dF0(2, 2) = 1;
        dFbar_dF(2, 2) = 1;
      }
      // tangent operator
      const auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
      const auto Kb =
          tfel::math::map<tfel::math::t2tot2<dimension, real>>(Kip.data());
      const auto K_F0 = tfel::math::eval(Kb * dFbar_dF0);
      const auto K_F = tfel::math::eval(Kb * dFbar_dF);
      const auto K_F0v = std::span<const real>(K_F0.data(), gsize * gsize);
      const auto K_Fv = std::span<const real>(K_F.data(), gsize * gsize);
      // assembly of the stiffness matrix
      for (size_type ni = 0; ni != nnodes; ++ni) {
        if constexpr (evaluateShapeFunctions) {
          child.updateStiffnessMatrix(Ke, K_F0v, shape0, dshape, dshape0, w,
                                      ni);
          child.updateStiffnessMatrix(Ke, K_Fv, shape, dshape, dshape, w, ni);
        } else {
          child.updateStiffnessMatrix(Ke, K_F0v, dshape, dshape0, w, ni);
          child.updateStiffnessMatrix(Ke, K_Fv, dshape, dshape, w, ni);
        }
      }
    }
  }  // end of implementUpdateJacobian

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::integrateBatch(
      const IntegrationType it) {
    // the tangent operator blocks change
    this->fbar_data.has_tangent_operator = false;
    return BehaviourIntegratorBase::integrateBatch(it);
  }  // end of integrateBatch

  template <typename Child, Hypothesis H>
  void FBarBehaviourIntegratorCRTPBase<Child, H>::revert() {
    // data of the FBar formulation at the beginning of the time step
    batch::revertFBarData(this->fbar_data);
    BehaviourIntegratorBase::revert();
  }  // end of revert

  template <typename Child, Hypothesis H>
  void FBarBehaviourIntegratorCRTPBase<Child, H>::update() {
    // the data of the converged unknowns are the data at the beginning of
    // the next time step
    batch::updateFBarData(this->fbar_data);
    BehaviourIntegratorBase::update();
  }  // end of update

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::supportsKernels()
      const noexcept {
    // deformation gradient in 3D or in plane strain, computed by the
    // kernels as by the FBar behaviour integrator
    if (!batch::areGradientsKernelsAvailable()) {
      return false;
    }
    const auto k = batch::getGradientsKind(this->b);
    constexpr auto expected =
        (H == Hypothesis::PLANESTRAIN)
            ? batch::GradientsKind::DEFORMATION_GRADIENT_2D
            : batch::GradientsKind::DEFORMATION_GRADIENT_3D;
    if ((!k.has_value()) || (*k != expected)) {
      return false;
    }
    // all the elements must use the same integration rule
    return (this->quadrature_space->getNumberOfElements() == 0) ||
           (this->quadrature_space->getUniformIntegrationRule() != nullptr);
  }  // end of supportsKernels

  template <typename Child, Hypothesis H>
  std::vector<const mfem::IntegrationRule*> FBarBehaviourIntegratorCRTPBase<
      Child,
      H>::getAdditionalIntegrationRulesForKernels() const {
    const auto* const r0 =
        batch::getElementsCenterIntegrationRule(*(this->quadrature_space));
    if (r0 == nullptr) {
      return {};
    }
    return {r0};
  }  // end of getAdditionalIntegrationRulesForKernels

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::computeGradientsByKernels(
      Context& ctx,
      const batch::UnknownsDerivativesSet& d,
      const batch::ElementsGeometrySet& g) noexcept {
    auto data = KernelsData{};
    if (!this->getKernelsData(ctx, data, "computeGradientsByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    try {
      const auto p = d.find(data.rule);
      const auto* const r0 =
          batch::getElementsCenterIntegrationRule(*(this->quadrature_space));
      const auto p0 = d.find(r0);
      const auto g0 = g.find(r0);
      if ((p == d.end()) || (p0 == d.end()) || (g0 == g.end())) {
        return ctx.registerErrorMessage(
            "computeGradientsByKernels: the derivatives of the unknowns or the "
            "geometric data are not computed for the integration rules of the "
            "material (points and center of the elements)");
      }
      // deformation gradients, then gradients of the FBar formulation
      if (!batch::computeGradients(ctx, *this, p->second, *(data.geometry),
                                   data.kind)) {
        return false;
      }
      return batch::applyFBar(ctx, this->fbar_data, *this, p0->second,
                              g0->second, data.kind);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"computeGradientsByKernels: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "computeGradientsByKernels: unknown exception");
    }
  }  // end of computeGradientsByKernels

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::addInnerForcesByKernels(
      Context& ctx,
      mfem::Vector& f,
      const batch::ElementsGeometrySet& g) noexcept {
    // inner forces of the standard finite strain formulation
    return this->addInnerForcesFromThermodynamicForcesByKernels(ctx, f, g);
  }  // end of addInnerForcesByKernels

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::
      prepareFBarTangentOperatorForKernels(
          Context& ctx,
          KernelsData& data,
          const char* const name,
          const batch::ElementsGeometrySet& g) noexcept {
    if (!this->getKernelsData(ctx, data, name, g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    // safeguard: this case is refused by the setup of the non linear
    // evolution problem
    if (!this->fbar_data.valid) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the tangent operator of the FBar formulation can only be "
          "computed by kernels from gradients computed by kernels "
          "(GradientsComputation set to Kernels)");
    }
    if (this->fbar_data.has_tangent_operator) {
      return true;
    }
    try {
      if (this->fbar_data.coefficients.Size() != this->n) {
        // before the first computation of the gradients (prediction of the
        // first time step): null unknowns
        const auto* const r0 =
            batch::getElementsCenterIntegrationRule(*(this->quadrature_space));
        const auto g0 = g.find(r0);
        if (g0 == g.end()) {
          return ctx.registerErrorMessage(
              std::string{name} +
              ": the geometric data are not computed for the center of the "
              "elements");
        }
        if (!batch::initializeFBar(ctx, this->fbar_data, *this, g0->second,
                                   data.kind)) {
          return false;
        }
      }
      const auto use_device = mfem::Device::IsEnabled();
      return batch::computeFBarTangentOperator(
          ctx, this->fbar_data, *(this->quadrature_space),
          this->readTangentOperatorBlocks(use_device), *(data.geometry),
          data.kind);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{name} + ": " + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(std::string{name} +
                                      ": unknown exception");
    }
  }  // end of prepareFBarTangentOperatorForKernels

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::addStiffnessMatricesByKernels(
      Context& ctx,
      mfem::Vector& m,
      const batch::ElementsGeometrySet& g) noexcept {
    auto data = KernelsData{};
    if (!this->prepareFBarTangentOperatorForKernels(
            ctx, data, "addStiffnessMatricesByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    return batch::addFBarStiffnessMatrices(ctx, m, this->fbar_data,
                                           *(this->quadrature_space),
                                           *(data.geometry), data.kind);
  }  // end of addStiffnessMatricesByKernels

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::addStiffnessMatricesByKernels(
      Context& ctx,
      const batch::AssembledStiffnessMatrix& m,
      const batch::ElementsGeometrySet& g) noexcept {
    auto data = KernelsData{};
    if (!this->prepareFBarTangentOperatorForKernels(
            ctx, data, "addStiffnessMatricesByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    return batch::addFBarStiffnessMatrices(ctx, m, this->fbar_data,
                                           *(this->quadrature_space),
                                           *(data.geometry), data.kind);
  }  // end of addStiffnessMatricesByKernels

  template <typename Child, Hypothesis H>
  bool
  FBarBehaviourIntegratorCRTPBase<Child, H>::addTangentOperatorActionByKernels(
      Context& ctx,
      mfem::Vector& y,
      const mfem::Vector& x,
      const batch::ElementsGeometrySet& g) noexcept {
    auto data = KernelsData{};
    if (!this->prepareFBarTangentOperatorForKernels(
            ctx, data, "addTangentOperatorActionByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    return batch::addFBarTangentOperatorAction(ctx, y, this->fbar_data,
                                               *(this->quadrature_space), x,
                                               *(data.geometry), data.kind);
  }  // end of addTangentOperatorActionByKernels

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::
      addTangentOperatorDiagonalByKernels(
          Context& ctx,
          mfem::Vector& d,
          const batch::ElementsGeometrySet& g) noexcept {
    auto data = KernelsData{};
    if (!this->prepareFBarTangentOperatorForKernels(
            ctx, data, "addTangentOperatorDiagonalByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    return batch::addFBarStiffnessMatricesDiagonals(
        ctx, d, this->fbar_data, *(this->quadrature_space), *(data.geometry),
        data.kind);
  }  // end of addTangentOperatorDiagonalByKernels

  template <typename Child, Hypothesis H>
  FBarBehaviourIntegratorCRTPBase<Child,
                                  H>::~FBarBehaviourIntegratorCRTPBase() =
      default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_FBARBEHAVIOURINTEGRATORCRTPBASE_IXX */
