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
#include "mfem/general/device.hpp"
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
    const auto gsize = this->s1.gradients_stride;
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      // offset of the integration point
      const auto o = eoffset + i;
      this->computeGradientsAtIntegrationPoint(e, tr, u, ir.IntPoint(i), o,
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

  template <typename Child>
  void StandardBehaviourIntegratorCRTPBase<Child>::
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
      this->computeGradientsAtIntegrationPoint(e, tr, u, ir.IntPoint(i),
                                               eoffset + i, shape, dshape);
    }
    return true;
  }  // end of implementComputeGradients

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::supportsBatchIntegration()
      const noexcept {
    // the gradients and the tangent operator blocks of an orthotropic
    // behaviour are rotated at all the points by `integrateBatch`, with the
    // rotation functions of the batches for a behaviour generated by the
    // generic-parallel interface
    return true;
  }  // end of supportsBatchIntegration

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::computeGradients(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    return this->implementComputeGradients(e, tr, u);
  }  // end of computeGradients

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::supportsKernels()
      const noexcept {
    using Traits = BehaviourIntegratorTraits<Child>;
    // the kernels compute the gradients of the standard mechanical
    // behaviours from the derivatives of the unknowns only
    if constexpr ((Traits::updateExternalStateVariablesFromUnknownsValues) ||
                  (Traits::gradientsComputationRequiresShapeFunctions)) {
      return false;
    } else {
      if ((!batch::areGradientsKernelsAvailable()) ||
          (!batch::getGradientsKind(this->b).has_value())) {
        return false;
      }
      // the kernels compute the inner forces and the stiffness matrices from
      // the thermodynamic forces and the tangent operator blocks only: the
      // behaviour integrators whose residual or jacobian also depend on the
      // unknowns (regularisations, such as Faltus2026) are refused
      if ((this->requiresCurrentSolutionForResidualAssembly()) ||
          (this->requiresCurrentSolutionForJacobianAssembly())) {
        return false;
      }
      // all the elements must use the same integration rule
      return (this->quadrature_space->getNumberOfElements() == 0) ||
             (this->quadrature_space->getUniformIntegrationRule() != nullptr);
    }
  }  // end of supportsKernels

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::computeGradientsByKernels(
      Context& ctx,
      const batch::UnknownsDerivativesSet& d,
      const batch::ElementsGeometrySet& g) noexcept {
    const auto* const r = this->quadrature_space->getUniformIntegrationRule();
    if (r == nullptr) {
      if (this->quadrature_space->getNumberOfElements() == 0) {
        return true;
      }
      return ctx.registerErrorMessage(
          "computeGradientsByKernels: the elements of the material do not "
          "use the same integration rule");
    }
    const auto k = batch::getGradientsKind(this->b);
    if (!k.has_value()) {
      return ctx.registerErrorMessage(
          "computeGradientsByKernels: unsupported behaviour");
    }
    const auto p = d.find(r);
    const auto pg = g.find(r);
    if ((p == d.end()) || (pg == g.end())) {
      return ctx.registerErrorMessage(
          "computeGradientsByKernels: the derivatives of the unknowns or the "
          "geometric data are not computed for the integration rule of the "
          "material");
    }
    return batch::computeGradients(ctx, *this, p->second, pg->second, *k);
  }  // end of computeGradientsByKernels

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::addInnerForcesByKernels(
      Context& ctx,
      mfem::Vector& f,
      const batch::ElementsGeometrySet& g) noexcept {
    const auto* const r = this->quadrature_space->getUniformIntegrationRule();
    if (r == nullptr) {
      if (this->quadrature_space->getNumberOfElements() == 0) {
        return true;
      }
      return ctx.registerErrorMessage(
          "addInnerForcesByKernels: the elements of the material do not use "
          "the same integration rule");
    }
    const auto k = batch::getGradientsKind(this->b);
    if (!k.has_value()) {
      return ctx.registerErrorMessage(
          "addInnerForcesByKernels: unsupported behaviour");
    }
    const auto p = g.find(r);
    if (p == g.end()) {
      return ctx.registerErrorMessage(
          "addInnerForcesByKernels: the geometric data are not computed for "
          "the integration rule of the material");
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      const auto np = static_cast<size_type>(this->n);
      const auto nc = this->s1.thermodynamic_forces_stride;
      if (this->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {
        // thermodynamic forces in the global frame at all the points
        auto& v = this->wks.rotated_thermodynamic_forces;
        if (this->getRotatedThermodynamicForces().empty()) {
          this->computeThermodynamicForcesInGlobalFrame(v);
        }
        const auto s = std::span<const real>(
            v.Read(use_device), static_cast<std::size_t>(v.Size()));
        return batch::addInnerForces(
            ctx, f, *(this->quadrature_space),
            mgis::makePointByPointFieldView(s, np, nc, nc), p->second, *k);
      }
      const auto& v = this->getArrays().getArray(
          mfem_mgis::MaterialArrays::StateArray::THERMODYNAMIC_FORCES, ets);
      const auto s = std::span<const real>(v.Read(use_device),
                                           static_cast<std::size_t>(v.Size()));
      return batch::addInnerForces(
          ctx, f, *(this->quadrature_space),
          mgis::makeFieldView(s, this->getArrays().array_layout, np, nc),
          p->second, *k);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addInnerForcesByKernels: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addInnerForcesByKernels: unknown exception");
    }
  }  // end of addInnerForcesByKernels

  template <typename Child>
  template <typename Destination>
  bool StandardBehaviourIntegratorCRTPBase<Child>::
      implementAddStiffnessMatricesByKernels(
          Context& ctx,
          Destination& m,
          const batch::ElementsGeometrySet& g) noexcept {
    const auto* const r = this->quadrature_space->getUniformIntegrationRule();
    if (r == nullptr) {
      if (this->quadrature_space->getNumberOfElements() == 0) {
        return true;
      }
      return ctx.registerErrorMessage(
          "addStiffnessMatricesByKernels: the elements of the material do not "
          "use the same integration rule");
    }
    const auto k = batch::getGradientsKind(this->b);
    if (!k.has_value()) {
      return ctx.registerErrorMessage(
          "addStiffnessMatricesByKernels: unsupported behaviour");
    }
    const auto p = g.find(r);
    if (p == g.end()) {
      return ctx.registerErrorMessage(
          "addStiffnessMatricesByKernels: the geometric data are not computed "
          "for the integration rule of the material");
    }
    try {
      // the tangent operator blocks are in the global frame (rotated after
      // the integration of the behaviour for an orthotropic behaviour)
      const auto use_device = mfem::Device::IsEnabled();
      const auto& v = this->getArrays().getTangentOperatorBlocks();
      const auto s = std::span<const real>(v.Read(use_device),
                                           static_cast<std::size_t>(v.Size()));
      return batch::addStiffnessMatrices(
          ctx, m, *(this->quadrature_space),
          mgis::makeFieldView(s, this->getArrays().array_layout,
                              static_cast<size_type>(this->n), this->K_stride),
          p->second, *k);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addStiffnessMatricesByKernels: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addStiffnessMatricesByKernels: unknown exception");
    }
  }  // end of implementAddStiffnessMatricesByKernels

  template <typename Child>
  bool
  StandardBehaviourIntegratorCRTPBase<Child>::addStiffnessMatricesByKernels(
      Context& ctx,
      mfem::Vector& m,
      const batch::ElementsGeometrySet& g) noexcept {
    return this->implementAddStiffnessMatricesByKernels(ctx, m, g);
  }  // end of addStiffnessMatricesByKernels

  template <typename Child>
  bool
  StandardBehaviourIntegratorCRTPBase<Child>::addStiffnessMatricesByKernels(
      Context& ctx,
      const batch::AssembledStiffnessMatrix& m,
      const batch::ElementsGeometrySet& g) noexcept {
    return this->implementAddStiffnessMatricesByKernels(ctx, m, g);
  }  // end of addStiffnessMatricesByKernels

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
