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
    // the kernels compute the gradients of the standard mechanical and
    // stationary heat transfer behaviours from the derivatives of the
    // unknowns, the ones of the transient term of the heat equation from the
    // values of the shape functions (temperature at the point), and the ones
    // of the micromorphic damage behaviours from both; the
    // values of the unknowns at the integration points (temperature) are
    // also computed for the behaviours which update their external state
    // variables from them. The kind of the gradients of the behaviour must
    // match the computations of the behaviour integrator.
    const auto k = batch::getGradientsKind(this->b);
    if ((!k.has_value()) ||
        (batch::gradientsRequireUnknownsValues(*k) !=
         Traits::gradientsComputationRequiresShapeFunctions) ||
        (batch::gradientsRequireUnknownsDerivatives(*k) !=
         Traits::gradientsComputationRequiresShapeFunctionsDerivatives) ||
        (batch::areUnknownsValuesExternalStateVariables(*k) !=
         Traits::updateExternalStateVariablesFromUnknownsValues)) {
      return false;
    }
    // the kernels compute the inner forces and the stiffness matrices from
    // the thermodynamic forces and the tangent operator blocks only: the
    // behaviour integrators whose residual or jacobian also depend on the
    // unknowns (regularisations) are refused, unless they override this
    // method (FBar formulation and regularisation of Faltus et al.)
    if ((this->requiresCurrentSolutionForResidualAssembly()) ||
        (this->requiresCurrentSolutionForJacobianAssembly())) {
      return false;
    }
    return this->isQuadratureSpaceSupportedByKernels();
  }  // end of supportsKernels

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<
      Child>::requiresUnknownsValuesForKernels() const noexcept {
    using Traits = BehaviourIntegratorTraits<Child>;
    return Traits::updateExternalStateVariablesFromUnknownsValues ||
           Traits::gradientsComputationRequiresShapeFunctions;
  }  // end of requiresUnknownsValuesForKernels

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::computeGradientsByKernels(
      Context &ctx,
      const batch::UnknownsDerivativesSet &d,
      const batch::ElementsGeometrySet &g) noexcept {
    auto data = KernelsData{};
    if (!this->getKernelsData(ctx, data, "computeGradientsByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    const auto p = d.find(data.rule);
    if (p == d.end()) {
      return ctx.registerErrorMessage(
          "computeGradientsByKernels: the derivatives of the unknowns are not "
          "computed for the integration rule of the material");
    }
    if (!batch::computeGradients(ctx, *this, p->second, *(data.geometry),
                                 data.kind)) {
      return false;
    }
    using Traits = BehaviourIntegratorTraits<Child>;
    if constexpr (Traits::updateExternalStateVariablesFromUnknownsValues) {
      // values of the unknowns at the integration points, copied in the
      // external state variable of the behaviour, stored on the host (the
      // behaviours integrated on the device read the copies of the external
      // state variables, see `BehaviourIntegratorBase`)
      if (!batch::getUnknownsValues(ctx, this->unknowns_values,
                                    *(this->quadrature_space), p->second)) {
        return false;
      }
      try {
        const auto size =
            static_cast<std::size_t>(this->unknowns_values.Size());
        static_cast<Child&>(*this).setExternalStateVariablesFromUnknownsValues(
            std::span<const real>(this->unknowns_values.HostRead(), size));
      } catch (...) {
        return registerExceptionInErrorBacktrace(ctx);
      }
    }
    return true;
  }  // end of computeGradientsByKernels

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::addInnerForcesByKernels(
      Context &ctx,
      mfem::Vector &f,
      const batch::ElementsGeometrySet &g) noexcept {
    return this->addInnerForcesFromThermodynamicForcesByKernels(ctx, f, g);
  }  // end of addInnerForcesByKernels

  template <typename Child>
  template <typename Destination>
  bool StandardBehaviourIntegratorCRTPBase<Child>::
      implementAddStiffnessMatricesByKernels(
          Context &ctx,
          Destination &m,
          const batch::ElementsGeometrySet &g) noexcept {
    auto data = KernelsData{};
    if (!this->getKernelsData(ctx, data, "addStiffnessMatricesByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      return batch::addStiffnessMatrices(
          ctx, m, *(this->quadrature_space),
          this->readTangentOperatorBlocks(use_device), *(data.geometry),
          data.kind);
    } catch (std::exception &e) {
      return ctx.registerErrorMessage(
          std::string{"addStiffnessMatricesByKernels: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addStiffnessMatricesByKernels: unknown exception");
    }
  }  // end of implementAddStiffnessMatricesByKernels

  template <typename Child>
  bool
  StandardBehaviourIntegratorCRTPBase<Child>::addTangentOperatorActionByKernels(
      Context &ctx,
      mfem::Vector &y,
      const mfem::Vector &x,
      const batch::ElementsGeometrySet &g) noexcept {
    auto data = KernelsData{};
    if (!this->getKernelsData(ctx, data, "addTangentOperatorActionByKernels",
                              g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      return batch::addTangentOperatorAction(
          ctx, y, *(this->quadrature_space),
          this->readTangentOperatorBlocks(use_device), x, *(data.geometry),
          data.kind);
    } catch (std::exception &e) {
      return ctx.registerErrorMessage(
          std::string{"addTangentOperatorActionByKernels: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addTangentOperatorActionByKernels: unknown exception");
    }
  }  // end of addTangentOperatorActionByKernels

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::
      addTangentOperatorDiagonalByKernels(
          Context &ctx,
          mfem::Vector &d,
          const batch::ElementsGeometrySet &g) noexcept {
    auto data = KernelsData{};
    if (!this->getKernelsData(ctx, data, "addTangentOperatorDiagonalByKernels",
                              g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      return batch::addStiffnessMatricesDiagonals(
          ctx, d, *(this->quadrature_space),
          this->readTangentOperatorBlocks(use_device), *(data.geometry),
          data.kind);
    } catch (std::exception &e) {
      return ctx.registerErrorMessage(
          std::string{"addTangentOperatorDiagonalByKernels: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addTangentOperatorDiagonalByKernels: unknown exception");
    }
  }  // end of addTangentOperatorDiagonalByKernels

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
          if constexpr (
              (Traits::gradientsComputationRequiresShapeFunctions) &&
              (Traits::gradientsComputationRequiresShapeFunctionsDerivatives)) {
            // gradients depending on the values of the unknowns and on their
            // derivatives (micromorphic damage)
            child.updateInnerForces(Fe, rs, shape, dshape, w, ni);
          } else if constexpr (evaluateShapeFunctionsDerivatives) {
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
