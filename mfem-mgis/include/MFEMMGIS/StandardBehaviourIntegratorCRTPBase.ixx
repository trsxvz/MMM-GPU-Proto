/*!
 * \file   include/MFEMMGIS/StandardBehaviourIntegratorCRTPBase.ixx
 * \brief  This file implements the inline methods of the
 * `StandardBehaviourIntegratorCRTPBase` class
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#ifndef LIB_MFEM_MGIS_STANDARDBEHAVIOURINTEGRATORCRTPBASE_IXX
#define LIB_MFEM_MGIS_STANDARDBEHAVIOURINTEGRATORCRTPBASE_IXX

#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"

namespace mfem_mgis {

  template <typename Child>
  StandardBehaviourIntegratorCRTPBase<Child>::
      StandardBehaviourIntegratorCRTPBase(
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
  }  // end of StandardBehaviourIntegratorCRTPBase

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::implementIntegrate(
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
    //
    auto &child = static_cast<Child &>(*this);
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
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
    if constexpr (evaluateShapeFunctions) {
      shape.SetSize(e.GetDof());
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      dshape.SetSize(e.GetDof(), e.GetDim());
    }
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape.SetSize(e.GetDof());
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      this->dshape.SetSize(e.GetDof(), e.GetDim());
    }
#endif
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      // offset of the integration point
      const auto o = eoffset + i;
      const auto r = this->computeGradientsAtIntegrationPoint(
          e, tr, u, ir.IntPoint(i), o, shape, dshape);
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

  template <typename Child>
  auto StandardBehaviourIntegratorCRTPBase<Child>::
      computeGradientsAtIntegrationPoint(
          const mfem::FiniteElement &e,
          mfem::ElementTransformation &tr,
          const mfem::Vector &u,
          const mfem::IntegrationPoint &ip,
          const size_type o,
          mfem::Vector &shape_functions,
          mfem::DenseMatrix &shape_functions_derivatives) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    auto &child = static_cast<Child&>(*this);
    const auto nnodes = e.GetDof();
    const auto gsize = this->s1.gradients_stride;
    tr.SetIntPoint(&ip);
    if constexpr (evaluateShapeFunctions) {
      // get the gradients of the shape functions
      e.CalcPhysShape(tr, shape_functions);
    }
    if constexpr (Traits::updateExternalStateVariablesFromUnknownsValues) {
      child.updateExternalStateVariablesFromUnknownsValues(u, shape_functions,
                                                           o);
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      // get the gradients of the shape functions
      e.CalcPhysDShape(tr, shape_functions_derivatives);
    }
    auto g = this->s1.gradients.subspan(o * gsize, gsize);
    std::copy(this->macroscopic_gradients.begin(),
              this->macroscopic_gradients.end(), g.begin());
    for (size_type ni = 0; ni != nnodes; ++ni) {
      if constexpr (
          (Traits::gradientsComputationRequiresShapeFunctions) &&
          (Traits::gradientsComputationRequiresShapeFunctionsDerivatives)) {
        child.updateGradients(g, u, shape_functions,
                              shape_functions_derivatives, ni);
      } else if constexpr (
          Traits::gradientsComputationRequiresShapeFunctionsDerivatives) {
        child.updateGradients(g, u, shape_functions_derivatives, ni);
      } else {
        child.updateGradients(g, u, shape_functions, ni);
      }
    }
    const auto r = child.getRotationMatrix(o);
    child.rotateGradients(g, r);
    return r;
  }  // end of computeGradientsAtIntegrationPoint

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::implementComputeGradients(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    auto &child = static_cast<Child&>(*this);
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    const auto &ir = child.getIntegrationRule(e, tr);
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
    if constexpr (evaluateShapeFunctions) {
      shape.SetSize(e.GetDof());
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      dshape.SetSize(e.GetDof(), e.GetDim());
    }
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape.SetSize(e.GetDof());
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      this->dshape.SetSize(e.GetDof(), e.GetDim());
    }
#endif
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      static_cast<void>(this->computeGradientsAtIntegrationPoint(
          e, tr, u, ir.IntPoint(i), eoffset + i, shape, dshape));
    }
    return true;
  }  // end of implementComputeGradients

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::supportsBatchIntegration()
      const noexcept {
    // the behaviours generated by the generic-parallel interface only
    // provide the rotations of batches, not the rotations of a point used
    // by `computeGradientsAtIntegrationPoint` and `integrateBatch`
    return !((this->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) &&
             (this->b.batch_b != nullptr));
  }  // end of supportsBatchIntegration

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::computeGradients(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    return this->implementComputeGradients(e, tr, u);
  }  // end of computeGradients

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::integrateBatch(
      const IntegrationType it) {
    if (!this->integrateAllIntegrationPoints(it)) {
      return false;
    }
    if (it == IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR) {
      return true;
    }
    // rotation of the tangent operator blocks, as in `implementIntegrate`
    auto &child = static_cast<Child&>(*this);
    const auto np = static_cast<size_type>(this->n);
    for (size_type o = 0; o != np; ++o) {
      const auto r = child.getRotationMatrix(o);
      auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
      child.rotateTangentOperatorBlocks(Kip, r);
    }
    return true;
  }  // end of integrateBatch

  template <typename Child>
  void StandardBehaviourIntegratorCRTPBase<Child>::implementUpdateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of implementUpdateResidual

  template <typename Child>
  void StandardBehaviourIntegratorCRTPBase<Child>::implementComputeInnerForces(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    static_assert(
        ((evaluateShapeFunctions) || (evaluateShapeFunctionsDerivatives)),
        "neither the shape functions or their derivatives are required to "
        "compute the gradients, this is not supported");
    auto &child = static_cast<Child &>(*this);
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
    if constexpr (evaluateShapeFunctions) {
      shape.SetSize(e.GetDof());
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      dshape.SetSize(e.GetDof(), e.GetDim());
    }
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape.SetSize(e.GetDof());
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      this->dshape.SetSize(e.GetDof(), e.GetDim());
    }
#endif
    const auto nnodes = e.GetDof();
    const auto thsize = this->s1.thermodynamic_forces_stride;
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
      const auto r = child.getRotationMatrix(o);
      const auto s = this->s1.thermodynamic_forces.subspan(o * thsize, thsize);
      const auto &rs = child.rotateThermodynamicForces(s, r);
      for (size_type ni = 0; ni != nnodes; ++ni) {
        if constexpr (evaluateShapeFunctionsDerivatives) {
          child.updateInnerForces(Fe, rs, dshape, w, ni);
        } else {
          child.updateInnerForces(Fe, rs, shape, w, ni);
        }
      }
    }
  }  // end of implementComputeInnerForces

  template <typename Child>
  void StandardBehaviourIntegratorCRTPBase<Child>::implementUpdateJacobian(
      mfem::DenseMatrix &Ke,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr) {
    using Traits = BehaviourIntegratorTraits<Child>;
    auto &child = static_cast<Child &>(*this);
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    static_assert(
        ((evaluateShapeFunctions) || (evaluateShapeFunctionsDerivatives)),
        "neither the shape functions or their derivatives are required to "
        "compute the gradients, this is not supported");
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
    if constexpr (evaluateShapeFunctions) {
      shape.SetSize(e.GetDof());
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      dshape.SetSize(e.GetDof(), e.GetDim());
    }
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape.SetSize(e.GetDof());
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      this->dshape.SetSize(e.GetDof(), e.GetDim());
    }
#endif
    // element offset
    const auto nnodes = e.GetDof();
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
      if constexpr (evaluateShapeFunctionsDerivatives) {
        // get the derivatives of the shape functions
        e.CalcPhysDShape(tr, dshape);
      }
      // get the weights associated to point ip
      const auto w = child.getIntegrationPointWeight(tr, ip);
      // offset of the integration point
      const auto o = eoffset + i;
      const auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
      // assembly of the stiffness matrix
      for (size_type ni = 0; ni != nnodes; ++ni) {
        if constexpr ((evaluateShapeFunctions) &&
                      (evaluateShapeFunctionsDerivatives)) {
          child.updateStiffnessMatrix(Ke, Kip, shape, dshape, w, ni);
        } else if constexpr (evaluateShapeFunctionsDerivatives) {
          child.updateStiffnessMatrix(Ke, Kip, dshape, w, ni);
        } else {
          child.updateStiffnessMatrix(Ke, Kip, shape, w, ni);
        }
      }
    }
  }  // end of implementUpdateJacobian

  template <typename Child>
  StandardBehaviourIntegratorCRTPBase<
      Child>::~StandardBehaviourIntegratorCRTPBase() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_STANDARDBEHAVIOURINTEGRATORCRTPBASE_IXX */
