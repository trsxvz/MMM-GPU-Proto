/*!
 * \file   MFEMMGIS/Faltus2026RegularizedBehaviourIntegrators.ixx
 * \brief  This file implements the inline functions declared in
 * `MFEMMGIS/Faltus2026RegularizedBehaviourIntegrators.hxx`
 * \author Thomas Helfer
 * \date   19/03/2026
 */

#ifndef LIB_MFEMMGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_IXX
#define LIB_MFEMMGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_IXX

#include <array>
#include <string>
#include <vector>
#include <cstdint>
#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Batch/FBar.hxx"

namespace mfem_mgis {

  template <Hypothesis H>
  Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::
      Faltus2026RegularizedIsotropicBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr,
          const Parameters &params)
      : Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>(
            fed, m, std::move(b_ptr)),
        alpha(get<real>(may_throw, params, "PenalizationCoefficient")) {
    checkParameters(
        may_throw, params,
        std::map<std::string, std::string>{
            {"PenalizationCoefficient",
             "penalization coefficient of the difference between the "
             "deformation gradient at the centroid of the element and the "
             "deformation gradient at integration points"}});
  }  // end of Faltus2026RegularizedIsotropicBehaviourIntegrator

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<
      H>::requiresCurrentSolutionForResidualAssembly() const noexcept {
    return true;
  }  // end of requiresCurrentSolutionForResidualAssembly

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<
      H>::requiresCurrentSolutionForJacobianAssembly() const noexcept {
    return false;
  }  // end of requiresCurrentSolutionForJacobianAssembly

  template <Hypothesis H>
  void Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::updateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>::updateResidual(
        Fe, e, tr, u);
    constexpr auto gsize = []() constexpr {
      if constexpr ((H == Hypothesis::PLANESTRAIN) ||
                    (H == Hypothesis::PLANESTRESS)) {
        return size_type{5};
      } else {
        return size_type{9};
      }
    }
    ();
#ifdef MFEM_THREAD_SAFE
    mfem::DenseMatrix dshape, dshape0;
    dshape.SetSize(e.GetDof(), e.GetDim());
    dshape0.SetSize(e.GetDof(), e.GetDim());
#else
    this->dshape0.SetSize(e.GetDof(), e.GetDim());
    // this->dshape is updated by
    // Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>
#endif
    const auto nnodes = e.GetDof();
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    auto F0 = std::array<real, gsize>{};
    std::copy(this->macroscopic_gradients.begin(),
              this->macroscopic_gradients.end(), F0.begin());
    {  // computation of the mean value of the deformation gradient
      auto F0v = std::span<real>(F0);
      const auto &ir = mfem::IntRules.Get(e.GetGeomType(), 0);
      const auto &ip = ir.IntPoint(0);
      tr.SetIntPoint(&ip);
      e.CalcPhysDShape(tr, dshape0);
      for (size_type ni = 0; ni != nnodes; ++ni) {
        this->updateGradients(F0v, u, dshape0, ni);
      }
    }
    // contribution to the residual
    auto dF = std::array<real, gsize>{};
    const auto &ir = this->getIntegrationRule(e, tr);
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      // get the gradients of the shape functions
#ifdef MFEM_THREAD_SAFE
      e.CalcPhysDShape(tr, dshape);
#else
      e.CalcPhysDShape(tr, this->dshape);
#endif
      // offset of the integration point
      const auto o = eoffset + i;
      // current estimation of the deformation gradient at the end of the time
      // step
      auto F = this->s1.gradients.subspan(o * gsize, gsize);
      for (size_type c = 0; c != gsize; ++c) {
        dF[c] = this->alpha * (F[c] - F0[c]);
      }
      const auto w = this->getIntegrationPointWeight(tr, ip);
      for (size_type ni = 0; ni != nnodes; ++ni) {
#ifdef MFEM_THREAD_SAFE
        this->updateInnerForces(Fe, dF, dshape, w, ni);
#else
        this->updateInnerForces(Fe, dF, this->dshape, w, ni);
#endif
        this->updateInnerForces(Fe, dF, dshape0, -w, ni);
      }
    }
  }  // end of updateResidual

  template <Hypothesis H>
  void Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::updateJacobian(
      mfem::DenseMatrix &Je,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>::updateJacobian(
        Je, e, tr, u);
    constexpr auto gsize = []() constexpr {
      if constexpr ((H == Hypothesis::PLANESTRAIN) ||
                    (H == Hypothesis::PLANESTRESS)) {
        return size_type{5};
      } else {
        return size_type{9};
      }
    }
    ();
#ifdef MFEM_THREAD_SAFE
    mfem::DenseMatrix dshape, dshape0;
    dshape.SetSize(e.GetDof(), e.GetDim());
    dshape0.SetSize(e.GetDof(), e.GetDim());
#else
    this->dshape0.SetSize(e.GetDof(), e.GetDim());
    // this->dshape is updated by
    // Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>
#endif
    {  // computation of the mean value of the deformation gradient
      const auto &ir = mfem::IntRules.Get(e.GetGeomType(), 0);
      const auto &ip = ir.IntPoint(0);
      tr.SetIntPoint(&ip);
      e.CalcPhysDShape(tr, dshape0);
    }
    const auto &ir = this->getIntegrationRule(e, tr);
    auto Km = std::array<real, gsize * gsize>{};
    for (size_type i = 0; i != gsize; ++i) {
      for (size_type j = 0; j != gsize; ++j) {
        Km[i * gsize + j] = 0;
      }
      Km[i * gsize + i] = this->alpha;
    }
    const auto nnodes = e.GetDof();
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      // get the gradients of the shape functions
#ifdef MFEM_THREAD_SAFE
      e.CalcPhysDShape(tr, dshape);
#else
      e.CalcPhysDShape(tr, this->dshape);
#endif
      const auto w = this->getIntegrationPointWeight(tr, ip);
      for (size_type ni = 0; ni != nnodes; ++ni) {
#ifdef MFEM_THREAD_SAFE
        this->updateStiffnessMatrix(Je, Km, dshape, dshape, w, ni);
        this->updateStiffnessMatrix(Je, Km, dshape0, dshape0, w, ni);
        this->updateStiffnessMatrix(Je, Km, dshape, dshape0, -w, ni);
        this->updateStiffnessMatrix(Je, Km, dshape0, dshape, -w, ni);
#else
        this->updateStiffnessMatrix(Je, Km, this->dshape, this->dshape, w, ni);
        this->updateStiffnessMatrix(Je, Km, this->dshape0, this->dshape0, w,
                                    ni);
        this->updateStiffnessMatrix(Je, Km, this->dshape, this->dshape0, -w,
                                    ni);
        this->updateStiffnessMatrix(Je, Km, this->dshape0, this->dshape, -w,
                                    ni);
#endif
      }
    }
  }  // end of updateJacobian

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::computeGradients(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& u) {
    this->faltus_data.valid = false;
    return Faltus2026RegularizedIsotropicBehaviourIntegratorBase<
        H>::computeGradients(e, tr, u);
  }  // end of computeGradients

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::integrate(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& u,
      const IntegrationType it) {
    this->faltus_data.valid = false;
    return Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>::integrate(
        e, tr, u, it);
  }  // end of integrate

  template <Hypothesis H>
  void Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::revert() {
    // deformation gradients at the center of the elements at the beginning
    // of the time step
    batch::revertFaltusData(this->faltus_data);
    Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>::revert();
  }  // end of revert

  template <Hypothesis H>
  void Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::update() {
    batch::updateFaltusData(this->faltus_data);
    Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>::update();
  }  // end of update

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::supportsKernels()
      const noexcept {
    // in plane stress, the axial component of the deformation gradient is
    // not given by the unknowns
    if constexpr (H == Hypothesis::PLANESTRESS) {
      return false;
    } else {
      // deformation gradient in 3D or in plane strain
      const auto k = batch::getGradientsKind(this->b);
      constexpr auto expected =
          (H == Hypothesis::PLANESTRAIN)
              ? batch::GradientsKind::DEFORMATION_GRADIENT_2D
              : batch::GradientsKind::DEFORMATION_GRADIENT_3D;
      if ((!k.has_value()) || (*k != expected)) {
        return false;
      }
      return this->isQuadratureSpaceSupportedByKernels();
    }
  }  // end of supportsKernels

  template <Hypothesis H>
  std::vector<const mfem::IntegrationRule*>
  Faltus2026RegularizedIsotropicBehaviourIntegrator<
      H>::getAdditionalIntegrationRulesForKernels() const {
    const auto* const r0 =
        batch::getElementsCenterIntegrationRule(*(this->quadrature_space));
    if (r0 == nullptr) {
      return {};
    }
    return {r0};
  }  // end of getAdditionalIntegrationRulesForKernels

  template <Hypothesis H>
  bool
  Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::getFaltusKernelsData(
      Context& ctx,
      BehaviourIntegratorBase::KernelsData& data,
      const batch::ElementsGeometry*& g0,
      const char* const name,
      const batch::ElementsGeometrySet& g) const noexcept {
    g0 = nullptr;
    if (!this->getKernelsData(ctx, data, name, g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    const auto p0 = g.find(
        batch::getElementsCenterIntegrationRule(*(this->quadrature_space)));
    if (p0 == g.end()) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the geometric data are not computed for the center of the "
          "elements");
    }
    g0 = &(p0->second);
    return true;
  }  // end of getFaltusKernelsData

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::
      computeGradientsByKernels(Context& ctx,
                                const batch::UnknownsDerivativesSet& d,
                                const batch::ElementsGeometrySet& g) noexcept {
    // deformation gradients at the integration points
    if (!Faltus2026RegularizedIsotropicBehaviourIntegratorBase<
            H>::computeGradientsByKernels(ctx, d, g)) {
      return false;
    }
    // deformation gradients at the center of the elements
    auto data = BehaviourIntegratorBase::KernelsData{};
    const batch::ElementsGeometry* g0 = nullptr;
    if (!this->getFaltusKernelsData(ctx, data, g0, "computeGradientsByKernels",
                                    g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    const auto p0 = d.find(g0->rule);
    if (p0 == d.end()) {
      return ctx.registerErrorMessage(
          "computeGradientsByKernels: the derivatives of the unknowns are not "
          "computed for the center of the elements");
    }
    return batch::computeFaltusCenters(ctx, this->faltus_data, *this,
                                       &(p0->second), *g0, data.kind);
  }  // end of computeGradientsByKernels

  template <Hypothesis H>
  bool
  Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::addInnerForcesByKernels(
      Context& ctx,
      mfem::Vector& f,
      const batch::ElementsGeometrySet& g) noexcept {
    // inner forces of the standard finite strain formulation
    if (!Faltus2026RegularizedIsotropicBehaviourIntegratorBase<
            H>::addInnerForcesByKernels(ctx, f, g)) {
      return false;
    }
    auto data = BehaviourIntegratorBase::KernelsData{};
    const batch::ElementsGeometry* g0 = nullptr;
    if (!this->getFaltusKernelsData(ctx, data, g0, "addInnerForcesByKernels",
                                    g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    const auto gsize = static_cast<std::int64_t>(this->s1.gradients_stride);
    if ((this->faltus_data.valid) &&
        (static_cast<std::int64_t>(this->faltus_data.centers.Size()) !=
         gsize * this->quadrature_space->getNumberOfElements())) {
      // before the first computation of the gradients (prediction of the
      // first time step): null unknowns
      if (!batch::computeFaltusCenters(ctx, this->faltus_data, *this, nullptr,
                                       *g0, data.kind)) {
        return false;
      }
    }
    return batch::addFaltusInnerForces(ctx, f, this->faltus_data, *this,
                                       this->alpha, *(data.geometry), *g0,
                                       data.kind);
  }  // end of addInnerForcesByKernels

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<
      H>::addStiffnessMatricesByKernels(Context& ctx,
                                        mfem::Vector& m,
                                        const batch::ElementsGeometrySet&
                                            g) noexcept {
    if (!Faltus2026RegularizedIsotropicBehaviourIntegratorBase<
            H>::addStiffnessMatricesByKernels(ctx, m, g)) {
      return false;
    }
    auto data = BehaviourIntegratorBase::KernelsData{};
    const batch::ElementsGeometry* g0 = nullptr;
    if (!this->getFaltusKernelsData(ctx, data, g0,
                                    "addStiffnessMatricesByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    return batch::addFaltusStiffnessMatrices(ctx, m, *(this->quadrature_space),
                                             this->alpha, *(data.geometry), *g0,
                                             data.kind);
  }  // end of addStiffnessMatricesByKernels

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::
      addStiffnessMatricesByKernels(
          Context& ctx,
          const batch::AssembledStiffnessMatrix& m,
          const batch::ElementsGeometrySet& g) noexcept {
    if (!Faltus2026RegularizedIsotropicBehaviourIntegratorBase<
            H>::addStiffnessMatricesByKernels(ctx, m, g)) {
      return false;
    }
    auto data = BehaviourIntegratorBase::KernelsData{};
    const batch::ElementsGeometry* g0 = nullptr;
    if (!this->getFaltusKernelsData(ctx, data, g0,
                                    "addStiffnessMatricesByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    return batch::addFaltusStiffnessMatrices(ctx, m, *(this->quadrature_space),
                                             this->alpha, *(data.geometry), *g0,
                                             data.kind);
  }  // end of addStiffnessMatricesByKernels

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<
      H>::addTangentOperatorActionByKernels(Context& ctx,
                                            mfem::Vector& y,
                                            const mfem::Vector& x,
                                            const batch::ElementsGeometrySet&
                                                g) noexcept {
    if (!Faltus2026RegularizedIsotropicBehaviourIntegratorBase<
            H>::addTangentOperatorActionByKernels(ctx, y, x, g)) {
      return false;
    }
    auto data = BehaviourIntegratorBase::KernelsData{};
    const batch::ElementsGeometry* g0 = nullptr;
    if (!this->getFaltusKernelsData(ctx, data, g0,
                                    "addTangentOperatorActionByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    return batch::addFaltusTangentOperatorAction(
        ctx, y, *(this->quadrature_space), x, this->alpha, *(data.geometry),
        *g0, data.kind);
  }  // end of addTangentOperatorActionByKernels

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<
      H>::addTangentOperatorDiagonalByKernels(Context& ctx,
                                              mfem::Vector& d,
                                              const batch::ElementsGeometrySet&
                                                  g) noexcept {
    if (!Faltus2026RegularizedIsotropicBehaviourIntegratorBase<
            H>::addTangentOperatorDiagonalByKernels(ctx, d, g)) {
      return false;
    }
    auto data = BehaviourIntegratorBase::KernelsData{};
    const batch::ElementsGeometry* g0 = nullptr;
    if (!this->getFaltusKernelsData(ctx, data, g0,
                                    "addTangentOperatorDiagonalByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    return batch::addFaltusStiffnessMatricesDiagonals(
        ctx, d, *(this->quadrature_space), this->alpha, *(data.geometry), *g0,
        data.kind);
  }  // end of addTangentOperatorDiagonalByKernels

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<
      H>::addTangentOperatorNodalBlocksByKernels(Context& ctx,
                                                 mfem::Vector& b,
                                                 const batch::
                                                     ElementsGeometrySet&
                                                         g) noexcept {
    if (!Faltus2026RegularizedIsotropicBehaviourIntegratorBase<
            H>::addTangentOperatorNodalBlocksByKernels(ctx, b, g)) {
      return false;
    }
    auto data = BehaviourIntegratorBase::KernelsData{};
    const batch::ElementsGeometry* g0 = nullptr;
    if (!this->getFaltusKernelsData(
            ctx, data, g0, "addTangentOperatorNodalBlocksByKernels", g)) {
      return false;
    }
    if (data.rule == nullptr) {
      return true;
    }
    return batch::addFaltusStiffnessMatricesNodalBlocks(
        ctx, b, *(this->quadrature_space), this->alpha, *(data.geometry), *g0,
        data.kind);
  }  // end of addTangentOperatorNodalBlocksByKernels

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_IXX */
