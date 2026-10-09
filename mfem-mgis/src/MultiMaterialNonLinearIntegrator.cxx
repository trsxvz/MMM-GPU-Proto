/*!
 * \file   src/MultiMaterialNonLinearIntegrator.cxx
 * \brief  This file implements the `MultiMaterialNonLinearIntegrator` class
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <limits>
#include <string>
#include <cstdint>
#include <utility>
#include <algorithm>
#include <exception>
#include "MGIS/Raise.hxx"
#include "mfem/fem/linearform.hpp"
#include "mfem/fem/restriction.hpp"
#include "mfem/fem/quadinterpolator.hpp"
#include "mfem/fem/bilinearform.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/mesh/pmesh.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/BehaviourIntegratorFactory.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"

namespace mfem_mgis {

  /*!
   * \brief a simple test
   * \param[in] bis: behaviour integrators per material
   * \param[in] n: name of the calling method
   * \param[in] m: material id
   */
  static void checkIfBehaviourIntegratorsAreDefined(
      const std::vector<
          std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>& bis,
      const char* const n,
      const size_type m) {
    if ((bis.empty()) || (m >= static_cast<size_type>(bis.size()))) {
      abort("MultiMaterialNonLinearIntegrator::" + std::string(n) +
            ": no behaviour integrator associated with material '" +
            std::to_string(m) + "'");
    }
    if (bis.at(m).empty()) {
      abort("MultiMaterialNonLinearIntegrator::" + std::string(n) +
            ": no behaviour integrator associated with material '" +
            std::to_string(m) + "'");
    }
    // ok this is paranoïac, but does not hurt
    for (const auto& bi : bis.at(m)) {
      if (bi.get() == nullptr) {
        abort("MultiMaterialNonLinearIntegrator::" + std::string(n) +
              ": invalid behaviour integrator associated with material '" +
              std::to_string(m) + "'");
      }
    }
  }  // end if checkIfBehaviourIntegratorsAreDefined

  struct InternalForcesIntegrator final : public LinearFormIntegrator {
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] u: unknown
     * \param[in] fes: finite element space
     * \param[in] bis: behaviour integrators
     */
    InternalForcesIntegrator(
        const FiniteElementSpace<true>& fes,
        const mfem::Vector& u,
        std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
            bis)
        : pfespace(&fes), unknowns(u), behaviour_integrators(bis) {}
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] u: unknown
     * \param[in] fes: finite element space
     * \param[in] bis: behaviour integrators
     */
    InternalForcesIntegrator(
        const FiniteElementSpace<false>& fes,
        const mfem::Vector& u,
        std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
            bis)
        : fespace(&fes), unknowns(u), behaviour_integrators(bis) {}
    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& F) override {
      const auto m = tr.Attribute;
      checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                            "AssembleElementVector", m);
      auto update = [this, &e, &tr](mfem::Vector& Fe,
                                    AbstractBehaviourIntegrator& bi) {
        if (bi.requiresCurrentSolutionForResidualAssembly()) {
          auto vdofs = mfem::Array<int>{};
#ifdef MFEM_USE_MPI
          if (this->pfespace != nullptr) {
            this->pfespace->GetElementVDofs(tr.ElementNo, vdofs);
          } else {
            this->fespace->GetElementVDofs(tr.ElementNo, vdofs);
          }
#else  /* MFEM_USE_MPI */
          this->fespace->GetElementVDofs(tr.ElementNo, vdofs);
#endif /* MFEM_USE_MPI */
          this->unknowns.GetSubVector(vdofs, this->Ue);
        }
        bi.updateResidual(Fe, e, tr, this->Ue);
      };
      //
      const auto& integrators = this->behaviour_integrators.at(m);
      if (integrators.size() == 1) {
        update(F, *(integrators.front()));
        return;
      }
      const auto ndofs = [this, &e] {
#ifdef MFEM_USE_MPI
        if (this->pfespace != nullptr) {
          return e.GetDof() * (this->pfespace->GetVDim());
        }
#endif /* MFEM_USE_MPI */
        return e.GetDof() * (this->fespace->GetVDim());
      }();
      auto F_tmp = mfem::Vector{};
      F.SetSize(ndofs);
      // although not required as F_tmp will be resized and initialized by all
      // behaviour integrators, this line avoids repeated allocations
      F_tmp.SetSize(ndofs);
      F = 0.;
      for (auto& bi : integrators) {
        update(F_tmp, *bi);
        F += F_tmp;
      }
    }  // end of AssembleRHSElementVect

   private:
#ifdef MFEM_USE_MPI
    //! \brief pointer to the finite element space for parallel resolutions
    const FiniteElementSpace<true>* const pfespace = nullptr;
#endif /* MFEM_USE_MPI */
    //! \brief pointer to the finite element space for sequential resolutions
    const FiniteElementSpace<false>* const fespace = nullptr;
    //! \brief unknowns
    const mfem::Vector& unknowns;
    //! \brief temporary vector containing the unknown of the current element
    mfem::Vector Ue;
    //! \brief list of behaviour integrators
    std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
        behaviour_integrators;
  };

  struct StiffnessMatrixIntegrator final : public BilinearFormIntegrator {
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] u: unknown
     * \param[in] fes: finite element space
     * \param[in] bis: behaviour integrators
     */
    StiffnessMatrixIntegrator(
        const FiniteElementSpace<true>& fes,
        const mfem::Vector& u,
        std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
            bis)
        : pfespace(&fes), unknowns(u), behaviour_integrators(bis) {}
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] u: unknown
     * \param[in] fes: finite element space
     * \param[in] bis: behaviour integrators
     */
    StiffnessMatrixIntegrator(
        const FiniteElementSpace<false>& fes,
        const mfem::Vector& u,
        std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
            bis)
        : fespace(&fes), unknowns(u), behaviour_integrators(bis) {}
    //
    void AssembleElementMatrix(const mfem::FiniteElement& e,
                               mfem::ElementTransformation& tr,
                               mfem::DenseMatrix& K) override {
      const auto m = tr.Attribute;
      checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                            "AssembleElementGrad", m);
      auto update = [this, &e, &tr](mfem::DenseMatrix& Ke,
                                    AbstractBehaviourIntegrator& bi) {
        if (bi.requiresCurrentSolutionForJacobianAssembly()) {
          auto vdofs = mfem::Array<int>{};
#ifdef MFEM_USE_MPI
          if (this->pfespace != nullptr) {
            this->pfespace->GetElementVDofs(tr.ElementNo, vdofs);
          } else {
            this->fespace->GetElementVDofs(tr.ElementNo, vdofs);
          }
#else  /* MFEM_USE_MPI */
          this->fespace->GetElementVDofs(tr.ElementNo, vdofs);
#endif /* MFEM_USE_MPI */
          this->unknowns.GetSubVector(vdofs, this->Ue);
        }
        bi.updateJacobian(Ke, e, tr, this->Ue);
      };
      //
      const auto& integrators = this->behaviour_integrators.at(m);
      if (integrators.size() == 1) {
        update(K, *(integrators.front()));
        return;
      }
      const auto ndofs = [this, &e] {
#ifdef MFEM_USE_MPI
        if (this->pfespace != nullptr) {
          return e.GetDof() * (this->pfespace->GetVDim());
        }
#endif /* MFEM_USE_MPI */
        return e.GetDof() * (this->fespace->GetVDim());
      }();
      auto K_tmp = mfem::DenseMatrix{};
      K.SetSize(ndofs, ndofs);
      // although not required as K_tmp will be resized and initialized by all
      // behaviour integrators, this line avoids repeated allocations
      K_tmp.SetSize(ndofs, ndofs);
      K = 0.;
      for (auto& bi : integrators) {
        update(K_tmp, *bi);
        K += K_tmp;
      }
    }  // end of AssembleElementMatrix

   private:
#ifdef MFEM_USE_MPI
    //! \brief pointer to the finite element space for parallel resolutions
    const FiniteElementSpace<true>* const pfespace = nullptr;
#endif /* MFEM_USE_MPI */
    //! \brief pointer to the finite element space for sequential resolutions
    const FiniteElementSpace<false>* const fespace = nullptr;
    //! \brief unknowns
    const mfem::Vector& unknowns;
    //! \brief temporary vector containing the unknown of the current element
    mfem::Vector Ue;
    //! \brief list of behaviour integrators
    std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
        behaviour_integrators;
  };  // end of StiffnessMatrixIntegrator

  MultiMaterialNonLinearIntegrator::MultiMaterialNonLinearIntegrator(
      std::shared_ptr<const FiniteElementDiscretization> fed,
      const Hypothesis h)
      : fe_discretization(fed), hypothesis(h) {
    if (this->fe_discretization->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto& mesh = this->fe_discretization->getMesh<true>();
      // shifting by one allows to directly use the material id to get
      // the behaviour integrators in all the other methods.
      // However, behaviour_integrators[0] will always be null.
      this->behaviour_integrators.resize(mesh.attributes.Max() + 1);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      const auto& mesh = this->fe_discretization->getMesh<false>();
      // shifting by one allows to directly use the material id to get
      // the behaviour integrators in all the other methods.
      // However, behaviour_integrators[0] will always be empty.
      this->behaviour_integrators.resize(mesh.attributes.Max() + 1);
    }
  }  // end of MultiMaterialNonLinearIntegrator

  bool MultiMaterialNonLinearIntegrator::integrate(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      const IntegrationType it) {
    const auto m = tr.Attribute;
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "integrate", m);
    for (const auto& bi : this->behaviour_integrators.at(m)) {
      if (!bi->integrate(e, tr, U, it)) {
        return false;
      }
    }
    return true;
  }  // end of integrate

  bool MultiMaterialNonLinearIntegrator::checkBatchIntegrationSupport(
      Context& ctx) const noexcept {
    for (size_type m = 0;
         m != static_cast<size_type>(this->behaviour_integrators.size()); ++m) {
      for (const auto& bi : this->behaviour_integrators[m]) {
        if (!bi->supportsBatchIntegration()) {
          return ctx.registerErrorMessage(
              "the behaviour integrator of material '" + std::to_string(m) +
              "' does not support the integration of the behaviours by "
              "batches");
        }
      }
    }
    return true;
  }  // end of checkBatchIntegrationSupport

  bool MultiMaterialNonLinearIntegrator::checkPointByPointIntegrationSupport(
      Context& ctx) const noexcept {
    auto ctx2 = Context{};
    for (size_type m = 0;
         m != static_cast<size_type>(this->behaviour_integrators.size()); ++m) {
      for (const auto& bi : this->behaviour_integrators[m]) {
        if (!bi->hasMaterial()) {
          continue;
        }
        const auto om =
            static_cast<const AbstractBehaviourIntegrator&>(*bi).getMaterial(
                ctx2);
        if ((!isInvalid(om)) && (om->b.batch_b != nullptr)) {
          return ctx.registerErrorMessage(
              "the behaviour of material '" + std::to_string(m) +
              "' is generated by the generic-parallel interface: it can only "
              "be integrated by batches (parameter BehaviourIntegration set "
              "to Batch)");
        }
      }
    }
    return true;
  }  // end of checkPointByPointIntegrationSupport

  bool MultiMaterialNonLinearIntegrator::computeGradients(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U) {
    const auto m = tr.Attribute;
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "computeGradients", m);
    for (const auto& bi : this->behaviour_integrators.at(m)) {
      if (!bi->computeGradients(e, tr, U)) {
        return false;
      }
    }
    return true;
  }  // end of computeGradients

  bool MultiMaterialNonLinearIntegrator::integrateBatch(
      const IntegrationType it) {
    // all the behaviour integrators are integrated, even if one of them
    // fails, as all the points of a batch
    auto success = true;
    const auto host_readers = (!this->inner_forces_computed_by_kernels) ||
                              (!this->stiffness_matrices_computed_by_kernels);
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        if (!bi->integrateBatch(it)) {
          success = false;
        }
        if ((host_readers) && (bi->integratesBatchesOnDevice())) {
          auto ctx = Context{};
          auto om = bi->getMaterial(ctx);
          if (!isInvalid(om)) {
            om->synchronizeOnHostForReading();
          }
        }
      }
    }
    return success;
  }  // end of integrateBatch

  bool MultiMaterialNonLinearIntegrator::integratesBatchesOnDevice()
      const noexcept {
    auto empty = true;
    for (const auto& bis : this->behaviour_integrators) {
      for (const auto& bi : bis) {
        if (!bi->integratesBatchesOnDevice()) {
          return false;
        }
        empty = false;
      }
    }
    return !empty;
  }  // end of integratesBatchesOnDevice

  bool MultiMaterialNonLinearIntegrator::checkKernelsSupport(
      Context& ctx) const noexcept {
    for (size_type m = 0;
         m != static_cast<size_type>(this->behaviour_integrators.size()); ++m) {
      for (const auto& bi : this->behaviour_integrators[m]) {
        if (!bi->supportsKernels()) {
          return ctx.registerErrorMessage(
              "the gradients, the inner forces and the stiffness matrices of "
              "the behaviour integrator of material '" +
              std::to_string(m) +
              "' can't be computed by kernels (supported: standard "
              "mechanical behaviours in 3D, plane strain and plane stress, "
              "small or finite strain, without FBar nor regularization, all "
              "the elements of the material using the same integration "
              "rule)");
        }
      }
    }
    try {
      const auto& fes = [this]() -> const mfem::FiniteElementSpace& {
#ifdef MFEM_USE_MPI
        if (this->fe_discretization->describesAParallelComputation()) {
          return this->fe_discretization->getFiniteElementSpace<true>();
        }
#endif /* MFEM_USE_MPI */
        return this->fe_discretization->getFiniteElementSpace<false>();
      }();
      const auto& mesh = *(fes.GetMesh());
      if (mesh.GetNumGeometries(mesh.Dimension()) > 1) {
        return ctx.registerErrorMessage(
            "the computations by kernels do not support meshes mixing several "
            "geometries");
      }
      if (mesh.SpaceDimension() != mesh.Dimension()) {
        return ctx.registerErrorMessage(
            "the computations by kernels require the dimension of the space "
            "to be the dimension of the mesh");
      }
      if (!mfem::QuadratureInterpolator::SupportsFESpace(fes)) {
        return ctx.registerErrorMessage(
            "the finite element space is not supported by the quadrature "
            "interpolator of MFEM, used by the computations by kernels");
      }
      // all the elements must have a behaviour integrator (with the
      // computations element by element, `AssembleElementVector` reports an
      // error otherwise)
      for (int i = 0; i != mesh.attributes.Size(); ++i) {
        const auto a = static_cast<size_type>(mesh.attributes[i]);
        if ((a >= static_cast<size_type>(this->behaviour_integrators.size())) ||
            (this->behaviour_integrators[a].empty())) {
          return ctx.registerErrorMessage(
              "the computations by kernels require a behaviour integrator on "
              "each material, and material '" +
              std::to_string(a) + "' has none");
        }
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(e.what());
    } catch (...) {
      return ctx.registerErrorMessage("unknown exception");
    }
    return true;
  }  // end of checkKernelsSupport

  bool MultiMaterialNonLinearIntegrator::computeGradientsByKernels(
      Context& ctx,
      const mfem::FiniteElementSpace& fes,
      const mfem::Vector& u) noexcept {
    // derivatives of the unknowns, computed once for each integration rule
    auto rules = std::vector<const mfem::IntegrationRule*>{};
    for (const auto& bis : this->behaviour_integrators) {
      for (const auto& bi : bis) {
        const auto* const r =
            bi->getPartialQuadratureSpace().getUniformIntegrationRule();
        if ((r != nullptr) &&
            (std::find(rules.begin(), rules.end(), r) == rules.end())) {
          rules.push_back(r);
        }
      }
    }
    // the derivatives of the rules which are no longer used are released
    std::erase_if(this->unknowns_derivatives, [&rules](const auto& d) {
      return std::find(rules.begin(), rules.end(), d.first) == rules.end();
    });
    for (const auto* const r : rules) {
      if (!batch::computeUnknownsDerivatives(ctx, this->unknowns_derivatives[r],
                                             fes, u, *r)) {
        return false;
      }
    }
    if (!this->computeElementsGeometry(ctx, fes)) {
      return false;
    }
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        if (!bi->computeGradientsByKernels(ctx, this->unknowns_derivatives,
                                           this->elements_geometry)) {
          return false;
        }
      }
    }
    return true;
  }  // end of computeGradientsByKernels

  bool MultiMaterialNonLinearIntegrator::addInnerForcesByKernels(
      Context& ctx,
      const mfem::FiniteElementSpace& fes,
      mfem::Vector& f) noexcept {
    if (!this->computeElementsGeometry(ctx, fes)) {
      return false;
    }
    try {
      const auto* const restriction =
          fes.GetElementRestriction(mfem::ElementDofOrdering::NATIVE);
      auto& fe = this->inner_forces_element_values;
      fe.UseDevice(true);
      fe.SetSize(restriction->Height());
      fe = real{0};
      for (auto& bis : this->behaviour_integrators) {
        for (auto& bi : bis) {
          if (!bi->addInnerForcesByKernels(ctx, fe, this->elements_geometry)) {
            return false;
          }
        }
      }
      restriction->AddMultTranspose(fe, f);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"addInnerForcesByKernels: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addInnerForcesByKernels: unknown exception");
    }
    return true;
  }  // end of addInnerForcesByKernels

  void MultiMaterialNonLinearIntegrator::setInnerForcesComputedByKernels(
      const bool b) noexcept {
    this->inner_forces_computed_by_kernels = b;
  }  // end of setInnerForcesComputedByKernels

  bool MultiMaterialNonLinearIntegrator::computeElementsGeometry(
      Context& ctx, const mfem::FiniteElementSpace& fes) noexcept {
    // computed once for each integration rule (the mesh does not change)
    for (const auto& bis : this->behaviour_integrators) {
      for (const auto& bi : bis) {
        const auto* const r =
            bi->getPartialQuadratureSpace().getUniformIntegrationRule();
        if ((r != nullptr) && (!this->elements_geometry.contains(r))) {
          if (!batch::computeElementsGeometry(ctx, this->elements_geometry[r],
                                              fes, *r)) {
            // no partial entry is kept
            this->elements_geometry.erase(r);
            return false;
          }
        }
      }
    }
    return true;
  }  // end of computeElementsGeometry

  bool MultiMaterialNonLinearIntegrator::computeStiffnessMatricesByKernels(
      Context& ctx, const mfem::FiniteElementSpace& fes) noexcept {
    if (!this->computeElementsGeometry(ctx, fes)) {
      return false;
    }
    try {
      const auto* const fe = fes.GetTypicalFE();
      const auto ndofs =
          static_cast<std::int64_t>(fe != nullptr ? fe->GetDof() : 0) *
          fes.GetVDim();
      const auto size = ndofs * ndofs * fes.GetNE();
      if (size > std::numeric_limits<int>::max()) {
        return ctx.registerErrorMessage(
            "computeStiffnessMatricesByKernels: too many values for MFEM");
      }
      auto& m = this->stiffness_matrices;
      m.UseDevice(true);
      m.SetSize(static_cast<int>(size));
      m = real{0};
      for (auto& bis : this->behaviour_integrators) {
        for (auto& bi : bis) {
          if (!bi->addStiffnessMatricesByKernels(ctx, m,
                                                 this->elements_geometry)) {
            return false;
          }
        }
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"computeStiffnessMatricesByKernels: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "computeStiffnessMatricesByKernels: unknown exception");
    }
    return true;
  }  // end of computeStiffnessMatricesByKernels

  void MultiMaterialNonLinearIntegrator::setStiffnessMatricesComputedByKernels(
      const bool b) noexcept {
    this->stiffness_matrices_computed_by_kernels = b;
  }  // end of setStiffnessMatricesComputedByKernels

  bool MultiMaterialNonLinearIntegrator::computeElementsDegreesOfFreedom(
      Context& ctx, const mfem::FiniteElementSpace& fes) noexcept {
    if ((this->element_dofs.Size() != 0) || (fes.GetNE() == 0)) {
      return true;
    }
    return batch::computeElementsDegreesOfFreedom(ctx, this->element_dofs, fes);
  }  // end of computeElementsDegreesOfFreedom

  std::unique_ptr<mfem::SparseMatrix>
  MultiMaterialNonLinearIntegrator::makeAssembledStiffnessMatrix(
      Context& ctx, const mfem::FiniteElementSpace& fes) noexcept {
    if (!this->computeElementsDegreesOfFreedom(ctx, fes)) {
      return {};
    }
    return batch::makeAssembledStiffnessMatrix(ctx, fes, this->element_dofs);
  }  // end of makeAssembledStiffnessMatrix

  bool MultiMaterialNonLinearIntegrator::assembleStiffnessMatrixByKernels(
      Context& ctx,
      const mfem::FiniteElementSpace& fes,
      mfem::SparseMatrix& m) noexcept {
    if ((!this->computeElementsGeometry(ctx, fes)) ||
        (!this->computeElementsDegreesOfFreedom(ctx, fes))) {
      return false;
    }
    if (!batch::zeroSparseMatrixValues(ctx, m)) {
      return false;
    }
    const auto a = batch::AssembledStiffnessMatrix{
        .matrix = &m, .element_dofs = &(this->element_dofs)};
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        if (!bi->addStiffnessMatricesByKernels(ctx, a,
                                               this->elements_geometry)) {
          return false;
        }
      }
    }
    return true;
  }  // end of assembleStiffnessMatrixByKernels

  void MultiMaterialNonLinearIntegrator::AssembleElementVector(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      mfem::Vector& F) {
    if (this->inner_forces_computed_by_kernels) {
      F.SetSize(e.GetDof() * getNumberOfComponents(*(this->fe_discretization)));
      F = 0.;
      return;
    }
    const auto m = tr.Attribute;
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "AssembleElementVector", m);
    const auto& integrators = this->behaviour_integrators.at(m);
    if (integrators.size() == 1) {
      if (usePETSc()) {
        MFEM_VERIFY(
            integrators.front()->integrate(
                e, tr, U,
                IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR),
            "behaviour integration failed on element " << tr.ElementNo
                                                       << " of material " << m);
      }
      integrators.front()->updateResidual(F, e, tr, U);
      return;
    }
    const auto ndofs =
        e.GetDof() * getNumberOfComponents(*(this->fe_discretization));
    auto F_tmp = mfem::Vector{};
    F.SetSize(ndofs);
    // although not required as F_tmp will be resized and initialized by all
    // behaviour integrators, this line avoids repeated allocations
    F_tmp.SetSize(ndofs);
    F = 0.;
    for (const auto& bi : integrators) {
      if (usePETSc()) {
        MFEM_VERIFY(
            bi->integrate(
                e, tr, U,
                IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR),
            "behaviour integration failed on element " << tr.ElementNo
                                                       << " of material " << m);
      }
      bi->updateResidual(F_tmp, e, tr, U);
      F += F_tmp;
    }
  }  // end of AssembleElementVector

  void MultiMaterialNonLinearIntegrator::AssembleElementGrad(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      mfem::DenseMatrix& K) {
    if (this->stiffness_matrices_computed_by_kernels) {
      // matrix computed by `computeStiffnessMatricesByKernels`
      const auto ndofs =
          e.GetDof() * getNumberOfComponents(*(this->fe_discretization));
      K.SetSize(ndofs);
      const auto* const m = this->stiffness_matrices.HostRead() +
                            static_cast<std::size_t>(ndofs) *
                                static_cast<std::size_t>(ndofs) *
                                static_cast<std::size_t>(tr.ElementNo);
      std::copy(m, m + ndofs * ndofs, K.GetData());
      return;
    }
    const auto m = tr.Attribute;
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "AssembleElementGrad", m);
    const auto& integrators = this->behaviour_integrators.at(m);
    if (integrators.size() == 1) {
      integrators.front()->updateJacobian(K, e, tr, U);
      return;
    }
    const auto ndofs =
        e.GetDof() * getNumberOfComponents(*(this->fe_discretization));
    auto K_tmp = mfem::DenseMatrix{};
    K.SetSize(ndofs, ndofs);
    // although not required as K_tmp will be resized and initialized by all
    // behaviour integrators, this line avoids repeated allocations
    K_tmp.SetSize(ndofs, ndofs);
    K = 0.;
    for (const auto& bi : integrators) {
      bi->updateJacobian(K_tmp, e, tr, U);
      K += K_tmp;
    }
  }  // end of AssembleElementGrad

  std::optional<size_type>
  MultiMaterialNonLinearIntegrator::addBehaviourIntegrator(
      Context& ctx,
      const std::string& n,
      const size_type m,
      const std::string& l,
      const std::string& b,
      const Parameters& params) noexcept {
    if ((m == 0) || (m >= this->behaviour_integrators.size())) {
      return ctx.registerErrorMessage("invalid material index '" +
                                      std::to_string(m) + "'");
    }
    const auto of = BehaviourIntegratorFactory::get(ctx, this->hypothesis);
    if (isInvalid(of)) {
      return {};
    }
    auto& bis = this->behaviour_integrators[m];
    const auto s = static_cast<size_type>(bis.size());
    auto bptr = mfem_mgis::load(ctx, l, b, this->hypothesis);
    if (isInvalid(bptr)) {
      return {};
    }
    auto bi = of->generate(ctx, n, *(this->fe_discretization), m,
                           std::move(bptr), params);
    if (isInvalid(bi)) {
      return {};
    }
    bis.push_back(std::move(bi));
    return s;
  }  // end of addBehaviourIntegrator

  OptionalReference<const Material>
  MultiMaterialNonLinearIntegrator::getMaterial(
      Context& ctx, const size_type m, const size_type b) const noexcept {
    const auto obi = this->getBehaviourIntegrator(ctx, m, b);
    if (isInvalid(obi)) {
      return {};
    }
    const auto& om = obi->getMaterial(ctx);
    if (isInvalid(om)) {
      return {};
    }
    return {&(*om)};
  }  // end of getMaterial

  OptionalReference<Material> MultiMaterialNonLinearIntegrator::getMaterial(
      Context& ctx, const size_type m, const size_type b) noexcept {
    auto obi = this->getBehaviourIntegrator(ctx, m, b);
    if (isInvalid(obi)) {
      return {};
    }
    const auto& om = obi->getMaterial(ctx);
    if (isInvalid(om)) {
      return {};
    }
    return {&(*om)};
  }  // end of getMaterial

  std::optional<size_type>
  MultiMaterialNonLinearIntegrator::getNumberOfBehaviourIntegrators(
      Context& ctx, const size_type m) const noexcept {
    if ((m == 0) || (m >= this->behaviour_integrators.size())) {
      return ctx.registerErrorMessage("invalid material index '" +
                                      std::to_string(m) + "'");
    }
    const auto& bis = this->behaviour_integrators.at(m);
    return bis.size();
  }  // end of getBehaviourIntegrator

  OptionalReference<const AbstractBehaviourIntegrator>
  MultiMaterialNonLinearIntegrator::getBehaviourIntegrator(
      Context& ctx, const size_type m, const size_type b) const noexcept {
    if ((m == 0) || (m >= this->behaviour_integrators.size())) {
      return ctx.registerErrorMessage("invalid material index '" +
                                      std::to_string(m) + "'");
    }
    const auto& bis = this->behaviour_integrators.at(m);
    if (b >= bis.size()) {
      return ctx.registerErrorMessage("invalid behaviour index '" +
                                      std::to_string(b) + "'");
    }
    return {bis.at(b).get()};
  }  // end of getBehaviourIntegrator

  OptionalReference<AbstractBehaviourIntegrator>
  MultiMaterialNonLinearIntegrator::getBehaviourIntegrator(
      Context& ctx, const size_type m, const size_type b) noexcept {
    if ((m == 0) || (m >= this->behaviour_integrators.size())) {
      return ctx.registerErrorMessage("invalid material index '" +
                                      std::to_string(m) + "'");
    }
    auto& bis = this->behaviour_integrators.at(m);
    if (b >= bis.size()) {
      return ctx.registerErrorMessage("invalid behaviour index '" +
                                      std::to_string(b) + "'");
    }
    return {bis.at(b).get()};
  }  // end of getBehaviourIntegrator

  real MultiMaterialNonLinearIntegrator::getTimeIncrement() const noexcept {
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        if (bi != nullptr) {
          return bi->getTimeIncrement();
        }
      }
    }
    return real{};
  }  // end of setTimeIncrement

  void MultiMaterialNonLinearIntegrator::setTimeIncrement(const real dt) {
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        if (bi != nullptr) {
          bi->setTimeIncrement(dt);
        }
      }
    }
  }  // end of setTimeIncrement

  void MultiMaterialNonLinearIntegrator::setMacroscopicGradients(
      std::span<const real> g) {
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        bi->setMacroscopicGradients(g);
      }
    }
  }  // end of setMacroscopicGradients

  bool MultiMaterialNonLinearIntegrator::setup(Context& ctx,
                                               const real t,
                                               const real dt) noexcept {
    // note: synchronization among MPI processes is done later
    // in the calling nonlinear evolution problem
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        if (!bi->setup(ctx, t, dt)) {
          return false;
        }
      }
    }
    return true;
  }  // end of setup

  void MultiMaterialNonLinearIntegrator::synchronizeOnHost() {
    auto ctx = Context{};
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        if (!bi->hasMaterial()) {
          continue;
        }
        auto om = bi->getMaterial(ctx);
        if (!isInvalid(om)) {
          om->synchronizeOnHost();
        }
      }
    }
  }  // end of synchronizeOnHost

  void MultiMaterialNonLinearIntegrator::synchronizeOnHostForReading() const {
    auto ctx = Context{};
    for (const auto& bis : this->behaviour_integrators) {
      for (const auto& bi : bis) {
        if (!bi->hasMaterial()) {
          continue;
        }
        const auto om =
            static_cast<const AbstractBehaviourIntegrator&>(*bi).getMaterial(
                ctx);
        if (!isInvalid(om)) {
          om->synchronizeOnHostForReading();
        }
      }
    }
  }  // end of synchronizeOnHostForReading

  void MultiMaterialNonLinearIntegrator::revert() {
    // the materials synchronize their arrays (`Material::revertMaterialStates`)
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        bi->revert();
      }
    }
  }  // end of revert

  void MultiMaterialNonLinearIntegrator::update() {
    // the materials synchronize their arrays (`Material::updateMaterialStates`)
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        bi->update();
      }
    }
  }  // end of update

  std::vector<size_type>
  MultiMaterialNonLinearIntegrator::getAssignedMaterialsIdentifiers() const {
    std::vector<size_type> mids;
    for (size_type i = 0; i != this->behaviour_integrators.size(); ++i) {
      if (!this->behaviour_integrators[i].empty()) {
        mids.push_back(i);
      }
    }
    return mids;
  }  // end of getAssignedMaterialsIdentifiers

  LinearizedOperators MultiMaterialNonLinearIntegrator::getLinearizedOperators(
      const mfem::Vector& u) {
    if (this->fe_discretization->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto& fespace =
          this->fe_discretization->getFiniteElementSpace<true>();
      return LinearizedOperators{
          .K = std::make_unique<StiffnessMatrixIntegrator>(
              fespace, u, this->behaviour_integrators),
          .mFi = std::make_unique<InternalForcesIntegrator>(
              fespace, u, this->behaviour_integrators)};
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    const auto& fespace =
        this->fe_discretization->getFiniteElementSpace<false>();
    return LinearizedOperators{
        .K = std::make_unique<StiffnessMatrixIntegrator>(
            fespace, u, this->behaviour_integrators),
        .mFi = std::make_unique<InternalForcesIntegrator>(
            fespace, u, this->behaviour_integrators)};
  }  // end of getLinearizedOperators

  MultiMaterialNonLinearIntegrator::~MultiMaterialNonLinearIntegrator() =
      default;

}  // end of namespace mfem_mgis
