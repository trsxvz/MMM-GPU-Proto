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
#include "mfem/general/device.hpp"
#include "mfem/fem/linearform.hpp"
#include "mfem/fem/restriction.hpp"
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
#include "ExecutionSettings.hxx"

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

  bool MultiMaterialNonLinearIntegrator::checkStiffnessMatricesKernelsSupport(
      Context& ctx, const bool gradients_kernels) const noexcept {
    if (gradients_kernels) {
      return true;
    }
    for (size_type m = 0;
         m != static_cast<size_type>(this->behaviour_integrators.size()); ++m) {
      for (const auto& bi : this->behaviour_integrators[m]) {
        if (bi->requiresCurrentSolutionForJacobianAssembly()) {
          return ctx.registerErrorMessage(
              "the stiffness matrices of the behaviour integrator of "
              "material '" +
              std::to_string(m) +
              "' depend on the current estimate of the solution (FBar "
              "formulation): they can only be computed by kernels from "
              "gradients computed by kernels (GradientsComputation set to "
              "Kernels)");
        }
      }
    }
    return true;
  }  // end of checkStiffnessMatricesKernelsSupport

  bool MultiMaterialNonLinearIntegrator::checkInnerForcesKernelsSupport(
      Context& ctx, const bool gradients_kernels) const noexcept {
    if (gradients_kernels) {
      return true;
    }
    for (size_type m = 0;
         m != static_cast<size_type>(this->behaviour_integrators.size()); ++m) {
      for (const auto& bi : this->behaviour_integrators[m]) {
        if (bi->requiresCurrentSolutionForResidualAssembly()) {
          return ctx.registerErrorMessage(
              "the inner forces of the behaviour integrator of material '" +
              std::to_string(m) +
              "' depend on the current estimate of the solution "
              "(regularization of Faltus et al. for instance): they can only "
              "be computed "
              "by kernels from gradients computed by kernels "
              "(GradientsComputation set to Kernels)");
        }
      }
    }
    return true;
  }  // end of checkInnerForcesKernelsSupport

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
              "mechanical behaviours, small or finite strain, stationary non "
              "linear heat transfer behaviours, transient term of the heat "
              "equation and micromorphic damage behaviours, in 3D, plane "
              "strain and plane stress, the regularizations FBar and of "
              "Faltus et al. in 3D and in plane strain only, all the "
              "elements of the material using the same integration rule)");
        }
      }
    }
    try {
      const auto& fes = getFiniteElementSpace(*(this->fe_discretization));
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
      // the kernels use the gradients of the shape functions on the
      // reference element (`mfem::DofToQuad`) and the element restriction of
      // MFEM
      const auto* const fe = (fes.GetNE() != 0) ? fes.GetTypicalFE() : nullptr;
      if ((fes.IsVariableOrder()) || (mesh.NURBSext != nullptr) ||
          ((fe != nullptr) &&
           (fe->GetMapType() != mfem::FiniteElement::MapType::VALUE))) {
        return ctx.registerErrorMessage(
            "the computations by kernels require a finite element space of "
            "uniform order whose shape functions are mapped by value (H1), on "
            "a mesh which is not a NURBS mesh");
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

  std::vector<const mfem::IntegrationRule*>
  MultiMaterialNonLinearIntegrator::getIntegrationRules() const {
    auto rules = std::vector<const mfem::IntegrationRule*>{};
    auto add = [&rules](const mfem::IntegrationRule* const r) {
      if ((r != nullptr) &&
          (std::find(rules.begin(), rules.end(), r) == rules.end())) {
        rules.push_back(r);
      }
    };
    for (const auto& bis : this->behaviour_integrators) {
      for (const auto& bi : bis) {
        add(bi->getPartialQuadratureSpace().getUniformIntegrationRule());
        for (const auto* const r :
             bi->getAdditionalIntegrationRulesForKernels()) {
          add(r);
        }
      }
    }
    return rules;
  }  // end of getIntegrationRules

  bool MultiMaterialNonLinearIntegrator::computeGradientsByKernels(
      Context& ctx,
      const mfem::FiniteElementSpace& fes,
      const mfem::Vector& u) noexcept {
    // derivatives of the unknowns, computed once for each integration rule
    const auto rules = this->getIntegrationRules();
    // the derivatives of the rules which are no longer used are released
    std::erase_if(this->unknowns_derivatives, [&rules](const auto& d) {
      return std::find(rules.begin(), rules.end(), d.first) == rules.end();
    });
    for (const auto* const r : rules) {
      // the values of the unknowns at the points are computed if a behaviour
      // integrator of this rule requires them
      auto values = false;
      for (const auto& bis : this->behaviour_integrators) {
        for (const auto& bi : bis) {
          values =
              values ||
              ((bi->getPartialQuadratureSpace().getUniformIntegrationRule() ==
                r) &&
               (bi->requiresUnknownsValuesForKernels()));
        }
      }
      if (!batch::computeUnknownsDerivatives(ctx, this->unknowns_derivatives[r],
                                             fes, u, *r, values)) {
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
    try {
      for (const auto* const r : this->getIntegrationRules()) {
        if (!this->elements_geometry.contains(r)) {
          if (!batch::computeElementsGeometry(ctx, this->elements_geometry[r],
                                              fes, *r)) {
            // no partial entry is kept
            this->elements_geometry.erase(r);
            return false;
          }
        }
      }
    } catch (...) {
      return registerExceptionInErrorBacktrace(ctx);
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

  void MultiMaterialNonLinearIntegrator::setOrderedStiffnessMatricesAssembly(
      const bool b) noexcept {
    this->ordered_stiffness_matrices_assembly = b;
  }  // end of setOrderedStiffnessMatricesAssembly

  size_type MultiMaterialNonLinearIntegrator::
      getNumberOfOrderedStiffnessMatricesAssemblies() const noexcept {
    return this->ordered_stiffness_matrices_assemblies;
  }  // end of getNumberOfOrderedStiffnessMatricesAssemblies

  bool MultiMaterialNonLinearIntegrator::addTangentOperatorActionByKernels(
      Context& ctx,
      const mfem::FiniteElementSpace& fes,
      const mfem::Vector& x,
      mfem::Vector& y) noexcept {
    if (!this->computeElementsGeometry(ctx, fes)) {
      return false;
    }
    try {
      const auto* const restriction =
          fes.GetElementRestriction(mfem::ElementDofOrdering::NATIVE);
      auto& xe = this->tangent_operator_element_arguments;
      xe.UseDevice(true);
      xe.SetSize(restriction->Height());
      restriction->Mult(x, xe);
      auto& ye = this->tangent_operator_element_values;
      ye.UseDevice(true);
      ye.SetSize(restriction->Height());
      ye = real{0};
      for (auto& bis : this->behaviour_integrators) {
        for (auto& bi : bis) {
          if (!bi->addTangentOperatorActionByKernels(ctx, ye, xe,
                                                     this->elements_geometry)) {
            return false;
          }
        }
      }
      restriction->AddMultTranspose(ye, y);
    } catch (...) {
      return registerExceptionInErrorBacktrace(ctx);
    }
    return true;
  }  // end of addTangentOperatorActionByKernels

  bool MultiMaterialNonLinearIntegrator::addTangentOperatorDiagonalByKernels(
      Context& ctx,
      const mfem::FiniteElementSpace& fes,
      mfem::Vector& d) noexcept {
    if (!this->computeElementsGeometry(ctx, fes)) {
      return false;
    }
    try {
      const auto* const restriction =
          fes.GetElementRestriction(mfem::ElementDofOrdering::NATIVE);
      auto& de = this->tangent_operator_element_values;
      de.UseDevice(true);
      de.SetSize(restriction->Height());
      de = real{0};
      for (auto& bis : this->behaviour_integrators) {
        for (auto& bi : bis) {
          if (!bi->addTangentOperatorDiagonalByKernels(
                  ctx, de, this->elements_geometry)) {
            return false;
          }
        }
      }
      restriction->AddMultTranspose(de, d);
    } catch (...) {
      return registerExceptionInErrorBacktrace(ctx);
    }
    return true;
  }  // end of addTangentOperatorDiagonalByKernels

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

  int MultiMaterialNonLinearIntegrator::
      getNumberOfStiffnessMatricesPositionsComputations(
          const bool sorted) const noexcept {
    return this->stiffness_matrices_positions_computations[sorted ? 1 : 0];
  }  // end of getNumberOfStiffnessMatricesPositionsComputations

  bool MultiMaterialNonLinearIntegrator::hasStiffnessMatricesPositions(
      const bool sorted) const noexcept {
    return this->stiffness_matrices_positions[sorted ? 1 : 0].Size() != 0;
  }  // end of hasStiffnessMatricesPositions

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
    // ordered assembly with a device (see
    // `setOrderedStiffnessMatricesAssembly`): stiffness matrices stored
    // element by element, then added value by value in a fixed order.
    // Without device, the stiffness matrices are added one element after
    // the other on the host, which is already a fixed order.
    if ((this->ordered_stiffness_matrices_assembly) &&
        (mfem::Device::IsEnabled())) {
      if (fes.GetNE() == 0) {
        ++(this->ordered_stiffness_matrices_assemblies);
        return true;
      }
      const auto n = this->element_dofs.Size() / fes.GetNE();
      if (std::int64_t{n} * n * fes.GetNE() > std::numeric_limits<int>::max()) {
        return ctx.registerErrorMessage(
            "assembleStiffnessMatrixByKernels: too many values to store the "
            "stiffness matrices of all the elements in the ordered assembly "
            "(see the TangentOperatorAssembly parameter)");
      }
      if (!this->computeStiffnessMatricesByKernels(ctx, fes)) {
        return false;
      }
      // contributions of the stiffness matrices of the elements to the
      // values of the matrix, kept for one matrix of each layout of the rows
      // (as the positions below), computed again for another matrix, or if
      // they do not match the matrix (another matrix at the same address)
      const auto layout = m.ColumnsAreSorted() ? 1 : 0;
      auto& contributions = this->stiffness_matrices_contributions[layout];
      auto& matrix = this->stiffness_matrices_contributions_matrices[layout];
      const auto compute_contributions = [&]() {
        ++(this->stiffness_matrices_contributions_computations[layout]);
        matrix = nullptr;
        if (!batch::computeEntriesContributions(ctx, contributions, m,
                                                this->element_dofs, n)) {
          return false;
        }
        matrix = m.GetJ();
        return true;
      };
      const auto a = batch::AssembledStiffnessMatrix{
          .matrix = &m, .element_dofs = &(this->element_dofs)};
      const auto computed =
          (matrix != m.GetJ()) ||
          (contributions.offsets.Size() != m.NumNonZeroElems() + 1);
      if (computed && (!compute_contributions())) {
        return false;
      }
      auto first_attempt = Context{};
      if (!batch::addStiffnessMatricesByEntries(computed ? ctx : first_attempt,
                                                a, contributions,
                                                this->stiffness_matrices)) {
        if (computed) {
          return false;
        }
        // the values were partially updated
        if ((!batch::zeroSparseMatrixValues(ctx, m)) ||
            (!compute_contributions()) ||
            (!batch::addStiffnessMatricesByEntries(ctx, a, contributions,
                                                   this->stiffness_matrices))) {
          return false;
        }
      }
      ++(this->ordered_stiffness_matrices_assemblies);
      return true;
    }
    // positions of the entries in the rows of the matrix, kept for each
    // layout of the rows (columns sorted, for example by a direct solver, or
    // diagonal first, see `batch::makeAssembledStiffnessMatrix`): computed
    // at the first assembly in this layout, and again after an assembly in
    // which some of them did not designate their columns (another matrix),
    // at most three times. They make the assembly 1.1 to 2.9 times faster on
    // the CPU, but gave no gain on a GPU, where they would use memory: they
    // are not used with a GPU or the debug device of MFEM. If their
    // computation fails, the entries are searched.
    constexpr auto maximum_number_of_positions_computations = 3;
    const auto layout = m.ColumnsAreSorted() ? 1 : 0;
    auto& positions = this->stiffness_matrices_positions[layout];
    auto& computations =
        this->stiffness_matrices_positions_computations[layout];
    const auto device = mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
    if ((!device) && (positions.Size() == 0) &&
        (computations < maximum_number_of_positions_computations)) {
      ++computations;
      auto positions_ctx = Context{};
      if (!batch::computeStiffnessMatricesPositions(
              positions_ctx, positions, fes, this->element_dofs, m)) {
        positions.DeleteAll();
        computations = maximum_number_of_positions_computations;
      }
    }
    const auto use_positions = (!device) && (positions.Size() != 0);
    auto searched_entries = std::size_t{};
    const auto a = batch::AssembledStiffnessMatrix{
        .matrix = &m,
        .element_dofs = &(this->element_dofs),
        .positions = use_positions ? &positions : nullptr,
        .searched_entries = &searched_entries};
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        if (!bi->addStiffnessMatricesByKernels(ctx, a,
                                               this->elements_geometry)) {
          return false;
        }
      }
    }
    if (searched_entries != 0) {
      // positions of another layout of the rows: computed again at the next
      // assembly, or abandoned
      if (computations < maximum_number_of_positions_computations) {
        positions.SetSize(0);
      } else {
        positions.DeleteAll();
      }
    }
    if (this->ordered_stiffness_matrices_assembly) {
      // no device: the elements were added one after the other
      ++(this->ordered_stiffness_matrices_assemblies);
    }
    return true;
  }  // end of assembleStiffnessMatrixByKernels

  int MultiMaterialNonLinearIntegrator::
      getNumberOfEntriesContributionsComputations(
          const bool sorted) const noexcept {
    return this->stiffness_matrices_contributions_computations[sorted ? 1 : 0];
  }  // end of getNumberOfEntriesContributionsComputations

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
    if (!internals::checkBehaviourExecutionLocation(ctx, *bptr)) {
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
    // the copies of the arrays of the materials on the host are valid after
    // the revert (see `BehaviourIntegratorBase::revert`)
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        bi->revert();
      }
    }
  }  // end of revert

  void MultiMaterialNonLinearIntegrator::update() {
    // the copies of the arrays of the materials on the host are valid after
    // the update (see `BehaviourIntegratorBase::update`)
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
