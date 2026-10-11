/*!
 * \file   L2Projection.cxx
 * \brief  This file implements the functions declared in
 * `MFEMMGIS/L2Projection.hxx`
 * \author Thomas Helfer
 * \date   14/01/2026
 */

#include <vector>
#include <exception>
#include <optional>
#include "mfem/mesh/mesh.hpp"
#include "mfem/fem/fespace.hpp"
#include "mfem/fem/fe/fe_base.hpp"
#include "mfem/fem/fe_coll.hpp"
#include "mfem/general/device.hpp"
#include "mfem/general/backends.hpp"
#include "mfem/mesh/submesh/submesh.hpp"
#include "mfem/fem/linearform.hpp"
#include "mfem/fem/bilinearform.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pfespace.hpp"
#include "mfem/mesh/submesh/psubmesh.hpp"
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/pbilinearform.hpp"
#endif

#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"
#include "MFEMMGIS/Batch/IntegrationPointsFunctions.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/L2Projection.hxx"
#include "PostProcessingsComputation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  [[nodiscard]] static bool checkConsistency_impl(
      Context& ctx,
      const FiniteElementDiscretization& fed0,
      const FiniteElementDiscretization& fed1) noexcept {
    //
    const auto& mesh0 = fed0.getMesh<parallel>();
    const auto& mesh1 = fed1.getMesh<parallel>();
    if (&mesh0 != &mesh1) {
      return ctx.registerErrorMessage(
          "inconsistent finite element discretization space: "
          "underling meshes are distinct");
    }
    //
    const auto& fes0 = fed0.getFiniteElementSpace<parallel>();
    const auto& fes1 = fed1.getFiniteElementSpace<parallel>();
    if ((fes0.FEColl() != fes1.FEColl()) ||
        (fes0.GetOrdering() != fes1.GetOrdering())) {
      return ctx.registerErrorMessage(
          "inconsistent finite element discretization space: "
          "underling finite element spaces are distinct");
    }
    return true;
  }  // end of checkConsistency_impl

  [[nodiscard]] static bool checkConsistency(
      Context& ctx,
      const FiniteElementDiscretization& fed0,
      const FiniteElementDiscretization& fed1) {
    const auto parallel = fed0.describesAParallelComputation();
    if (parallel) {
      if (!fed1.describesAParallelComputation()) {
        return ctx.registerErrorMessage(
            "inconsistent finite element discretization space: one is "
            "parallel, while the other is not");
      }
#ifdef MFEM_USE_MPI
      return checkConsistency_impl<true>(ctx, fed0, fed1);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    if (fed1.describesAParallelComputation()) {
      return ctx.registerErrorMessage(
          "inconsistent finite element discretization space: one is "
          "parallel, while the other is not");
    }
    return checkConsistency_impl<false>(ctx, fed0, fed1);
  }  // end of checkConsistency

  template <bool parallel>
  [[nodiscard]] static bool checkConsistency(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    if (fcts.empty()) {
      return true;
    }
    const auto& fed0 =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& mesh = fed0.getMesh<parallel>();
    const auto& elts_attributes = mesh.attributes;
    const auto on_boundaries =
        fcts.at(0).getPartialQuadratureSpace().isDefinedOnABoundary();
    auto ids = std::vector<size_type>();
    ids.reserve(fcts.size());
    for (const auto& f : fcts) {
      if (f.getNumberOfComponents() != fcts.front().getNumberOfComponents()) {
        return ctx.registerErrorMessage(
            "functions with different numbers of components given");
      }
      const auto& qspace = f.getPartialQuadratureSpace();
      const auto& fed = qspace.getFiniteElementDiscretization();
      if (!checkConsistency(ctx, fed0, fed)) {
        return false;
      }
      if (qspace.isDefinedOnABoundary() != on_boundaries) {
        return ctx.registerErrorMessage(
            "quadrature functions defined on materials and on boundaries can "
            "not be mixed");
      }
      const auto id = qspace.getId();
      if ((!on_boundaries) && (elts_attributes.Find(id) == -1)) {
        return ctx.registerErrorMessage("material id '" + std::to_string(id) +
                                        "' is not a element attribute");
      }
      if (std::find(ids.begin(), ids.end(), id) != ids.end()) {
        return ctx.registerErrorMessage(
            "several quadrature functions defined on '" +
            qspace.getLocationName() + "' given");
      }
      ids.push_back(id);
    }
    return true;
  }  // end of checkConsistency

  template <bool parallel>
  [[nodiscard]] static std::optional<bool> areDefinedOnWholeMesh_impl(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    if (fcts.empty()) {
      return ctx.registerErrorMessage(
          "empty list of partial quadrature functions given");
    }
    //
    if (!checkConsistency<parallel>(ctx, fcts)) {
      return {};
    }
    // functions defined on boundaries are never defined on the whole mesh
    if (fcts.at(0).getPartialQuadratureSpace().isDefinedOnABoundary()) {
      return false;
    }
    // thanks to checkConsistency we know that the element attributes contains
    // all the material ids, so we just need to check that the number of
    // attributes equals to number of functions
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& mesh = fed.getMesh<parallel>();
    const auto& elts_attributes = mesh.attributes;
    return elts_attributes.Size() == static_cast<size_type>(fcts.size());
  }  // end of areDefinedOnWholeMesh_impl

  [[nodiscard]] static std::optional<bool> areDefinedOnWholeMesh(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    if (fcts.empty()) {
      return ctx.registerErrorMessage(
          "empty list of partial quadrature functions given");
    }
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return areDefinedOnWholeMesh_impl<true>(ctx, fcts);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return areDefinedOnWholeMesh_impl<false>(ctx, fcts);
  }

  template <bool parallel>
  [[nodiscard]] static std::optional<L2ProjectionResult<parallel>>
  createL2ProjectionResult_impl(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    if (fcts.empty()) {
      return ctx.registerErrorMessage(
          "empty list of partial quadrature functions given");
    }
    if (!checkConsistency<parallel>(ctx, fcts)) {
      return {};
    }
    const auto onAll = areDefinedOnWholeMesh(ctx, fcts);
    if (isInvalid(onAll)) {
      return {};
    }
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& mesh = fed.getMesh<parallel>();
    auto r = L2ProjectionResult<parallel>{};
    if (*onAll) {
      auto oresult = makeGridFunction<parallel>(ctx, fcts, mesh);
      if (isInvalid(oresult)) {
        return {};
      }
      r.result.swap(oresult);
      return r;
    }
    // create sub mesh
    auto ids = std::set<size_type>{};
    for (const auto& f : fcts) {
      const auto& qspace = f.getPartialQuadratureSpace();
      const auto id = qspace.getId();
      if (!ids.insert(id).second) {
        return ctx.registerErrorMessage(
            "multiple quadrature functions defined on '" +
            qspace.getLocationName() + "'");
      }
    }
    // partial quadrature spaces defined on boundaries are built on the submesh
    // of the boundary
    const auto location =
        fcts.at(0).getPartialQuadratureSpace().isDefinedOnABoundary()
            ? MeshDiscretization::Location::ON_BOUNDARIES
            : MeshDiscretization::Location::ON_MATERIALS;
    r.submesh = fed.template getMutableSubMeshPointer<parallel>(
        ctx, Parameter::from(ids), location);
    if (isInvalid(r.submesh)) {
      return {};
    }
    //
    auto oresult = makeGridFunction<parallel>(ctx, fcts, *(r.submesh));
    if (isInvalid(oresult)) {
      return {};
    }
    r.result.swap(oresult);
    return r;
  }  // end of createL2ProjectionResult_impl

  template <>
  std::optional<L2ProjectionResult<true>> createL2ProjectionResult<true>(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
#ifdef MFEM_USE_MPI
    return createL2ProjectionResult_impl<true>(ctx, fcts);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of createL2ProjectionResult

  template <>
  std::optional<L2ProjectionResult<false>> createL2ProjectionResult<false>(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    return createL2ProjectionResult_impl<false>(ctx, fcts);
  }  // end of createL2ProjectionResult

  struct L2ProjectionFormIntegratorBase {
    L2ProjectionFormIntegratorBase(
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : functions(makeFunctionsMapping(throwing, fcts)) {}

   private:
    [[nodiscard]] static std::
        unordered_map<size_type, ImmutablePartialQuadratureFunctionView>
        makeFunctionsMapping(
            attributes::Throwing,
            const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
      auto functions =
          std::unordered_map<size_type,
                             ImmutablePartialQuadratureFunctionView>{};
      for (const auto& f : fcts) {
        const auto mid = f.getPartialQuadratureSpace().getId();
        if (!functions.insert({mid, f}).second) {
          raise("multiple functions defined for material '" +
                std::to_string(mid) + "'");
        }
      }
      return functions;
    }  // end of makeFunctionsMapping

   protected:
    // functions
    std::unordered_map<size_type, ImmutablePartialQuadratureFunctionView>
        functions;
    mfem::Vector shape;
  };

  struct ScalarL2ProjectionRHSFormIntegrator final
      : public L2ProjectionFormIntegratorBase,
        public mfem::LinearFormIntegrator {
    //
    using L2ProjectionFormIntegratorBase::L2ProjectionFormIntegratorBase;
    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& elvect) override {
      const auto m = tr.Attribute;
      const auto& f = this->functions.at(m);
      const int dof = e.GetDof();
      //
      shape.SetSize(dof);  // vector of size dof
      elvect.SetSize(dof);
      elvect = 0.0;
      //
      const auto& ir = f.getPartialQuadratureSpace().getIntegrationRule(e, tr);
      //
      for (int i = 0; i < ir.GetNPoints(); i++) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        const auto v =
            tr.Weight() * f.getIntegrationPointValue(tr.ElementNo, i);
        e.CalcPhysShape(tr, shape);
        add(elvect, ip.weight * v, shape, elvect);
      }
    }  // end of AssembleRHSElementVect
  };

  struct ScalarL2ProjectionRHSFormIntegratorII final
      : public L2ProjectionFormIntegratorBase,
        public mfem::LinearFormIntegrator {
    //
    ScalarL2ProjectionRHSFormIntegratorII(
        const mfem::Array<int>& m,
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : L2ProjectionFormIntegratorBase(fcts), elts_mapping(m) {}

    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& elvect) override {
      const auto m = tr.Attribute;
      const auto& f = this->functions.at(m);
      const int dof = e.GetDof();
      //
      shape.SetSize(dof);  // vector of size dof
      elvect.SetSize(dof);
      elvect = 0.0;
      //
      const auto& ir = f.getPartialQuadratureSpace().getIntegrationRule(e, tr);
      //
      for (int i = 0; i < ir.GetNPoints(); i++) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        const auto n = elts_mapping[tr.ElementNo];
        const auto v = tr.Weight() * f.getIntegrationPointValue(n, i);
        e.CalcPhysShape(tr, shape);
        add(elvect, ip.weight * v, shape, elvect);
      }
    }  // end of AssembleRHSElementVect

   private:
    const mfem::Array<int> elts_mapping;
  };

  struct ComponentL2ProjectionRHSFormIntegrator final
      : public L2ProjectionFormIntegratorBase,
        public mfem::LinearFormIntegrator {
    //
    ComponentL2ProjectionRHSFormIntegrator(
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
        const size_type c)
        : L2ProjectionFormIntegratorBase(fcts), component(c) {}
    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& elvect) override {
      const auto m = tr.Attribute;
      const auto& f = this->functions.at(m);
      const int dof = e.GetDof();
      //
      shape.SetSize(dof);  // vector of size dof
      elvect.SetSize(dof);
      elvect = 0.0;
      //
      const auto& ir = f.getPartialQuadratureSpace().getIntegrationRule(e, tr);
      //
      for (int i = 0; i < ir.GetNPoints(); i++) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        const auto fv =
            f.getIntegrationPointValues(tr.ElementNo, i)[this->component];
        const auto v = tr.Weight() * fv;
        e.CalcPhysShape(tr, shape);
        add(elvect, ip.weight * v, shape, elvect);
      }
    }  // end of AssembleRHSElementVect
   private:
    const size_type component;
  };

  struct ComponentL2ProjectionRHSFormIntegratorII final
      : public L2ProjectionFormIntegratorBase,
        public mfem::LinearFormIntegrator {
    //
    ComponentL2ProjectionRHSFormIntegratorII(
        const mfem::Array<int>& m,
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
        const size_type c)
        : L2ProjectionFormIntegratorBase(fcts), elts_mapping(m), component(c) {}

    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& elvect) override {
      const auto m = tr.Attribute;
      const auto& f = this->functions.at(m);
      const int dof = e.GetDof();
      //
      shape.SetSize(dof);  // vector of size dof
      elvect.SetSize(dof);
      elvect = 0.0;
      //
      const auto& ir = f.getPartialQuadratureSpace().getIntegrationRule(e, tr);
      //
      for (int i = 0; i < ir.GetNPoints(); i++) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        const auto n = elts_mapping[tr.ElementNo];
        const auto fv = f.getIntegrationPointValues(n, i)[this->component];
        const auto v = tr.Weight() * fv;
        e.CalcPhysShape(tr, shape);
        add(elvect, ip.weight * v, shape, elvect);
      }
    }  // end of AssembleRHSElementVect

   private:
    const mfem::Array<int> elts_mapping;
    const size_type component;
  };

  template <bool parallel>
  [[nodiscard]] static bool checkUpdateFunctionsArguments(
      Context& ctx,
      L2ProjectionResult<parallel>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    // checks
    if (l.linear_solver == nullptr) {
      return ctx.registerErrorMessage("uninitialized linear solver");
    }
    if (r.result.get() == nullptr) {
      return ctx.registerErrorMessage("uninitialized result");
    }
    // check consistency
    if (fcts.empty()) {
      return ctx.registerErrorMessage(
          "empty list of partial quadrature functions given");
    }
    if (!checkConsistency<parallel>(ctx, fcts)) {
      return {};
    }
    const auto& rfespace = *(
        static_cast<const FiniteElementSpace<parallel>*>(r.result->FESpace()));
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto fespaces_manager = fed.getFiniteElementSpacesManager();
    if (!fespaces_manager.manages(rfespace)) {
      return ctx.registerErrorMessage("inconsistent finite element spaces");
    }
    if (rfespace.GetVDim() != getNumberOfComponents(fcts.front())) {
      return ctx.registerErrorMessage("inconsistent number of components");
    }
    //
    const auto* const mesh = [&rfespace] {
      if constexpr (parallel) {
        return rfespace.GetParMesh();
      } else {
        return rfespace.GetMesh();
      }
    }();
    const auto& elts_attributes = mesh->attributes;
    for (int i = 0; i != elts_attributes.Size(); ++i) {
      const auto id = elts_attributes[i];
      const auto found = [&id, &fcts] {
        for (const auto& f : fcts) {
          const auto& qspace = f.getPartialQuadratureSpace();
          if (id == qspace.getId()) {
            return true;
          }
        }
        return false;
      }();
      if (!found) {
        return ctx.registerErrorMessage(
            "no partial quadrature function defined on material id '" +
            std::to_string(id) + "'");
      }
    }
    return true;
  }  // end of checkUpdateFunctionsArguments

  /*!
   * \brief data of the computation by kernels of the right-hand sides of the
   * L2 projections of functions defined at the integration points: geometric
   * data of the finite element space of the projection for the integration
   * rules of the functions, and values of the functions (see
   * `batch::gatherValues`)
   */
  struct L2ProjectionKernelsData {
    //! \brief geometric data, by integration rule
    batch::ElementsGeometrySet geometries;
    //! \brief values of the functions
    std::vector<mfem::Vector> values;
  };

  /*!
   * \return the data of the computation by kernels of the right-hand sides
   * of the L2 projections of functions, or an empty value if the
   * post-processings are not computed by kernels (see
   * `internals::usePostProcessingsKernels`) or if the kernels can't be used:
   * functions defined on boundaries or on a part of the mesh (submesh),
   * finite element space of another family than the Lagrange elements,
   * device of libCEED (whose assembly of the operators does not use the
   * element restrictions), elements of a function using different
   * integration rules, geometric data which can't be computed by kernels
   * (elements of different geometries). In parallel, the decision is made
   * by all the processes, which then assemble the operator in the same way.
   * \param[in] fes: scalar finite element space of the projection
   * \param[in] fcts: functions
   * \param[in] on_submesh: if true, the projection is made on a submesh
   */
  template <bool parallel>
  static std::optional<L2ProjectionKernelsData> makeL2ProjectionKernelsData(
      attributes::MayThrow,
      const FiniteElementSpace<parallel>& fes,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const bool on_submesh) {
    if ((on_submesh) || (!internals::usePostProcessingsKernels(may_throw))) {
      return {};
    }
    auto ctx = Context{};
    auto d = L2ProjectionKernelsData{};
    const auto ok = [&]() -> bool {
      const auto* const fec =
          dynamic_cast<const mfem::H1_FECollection*>(fes.FEColl());
      if ((fec == nullptr) ||
          (fec->GetBasisType() == mfem::BasisType::Positive) ||
          (mfem::Device::Allows(mfem::Backend::CEED_MASK))) {
        return false;
      }
      for (const auto& f : fcts) {
        const auto& qs = f.getPartialQuadratureSpace();
        auto& v = d.values.emplace_back();
        if (qs.getNumberOfElements() == 0) {
          continue;
        }
        const auto* const ir = qs.getUniformIntegrationRule();
        if (ir == nullptr) {
          return false;
        }
        auto& g = d.geometries[ir];
        if ((g.rule == nullptr) &&
            (!batch::computeElementsGeometry(ctx, g, fes, *ir))) {
          return false;
        }
        if (!batch::gatherValues(ctx, v, f)) {
          return false;
        }
      }
      return true;
    }();
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      auto lok = static_cast<int>(ok);
      auto gok = int{};
      MPI_Allreduce(&lok, &gok, 1, MPI_INT, MPI_MIN, fes.GetComm());
      if (gok == 0) {
        return {};
      }
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!ok) {
        return {};
      }
    }
    return d;
  }  // end of makeL2ProjectionKernelsData

  /*!
   * \brief computes by kernels the right-hand side of the L2 projection of a
   * component of functions (see `batch::addL2ProjectionRightHandSide`),
   * summed over the elements by the element restriction of MFEM
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] b: right-hand side
   * \param[in] fes: scalar finite element space of the projection
   * \param[in] fcts: functions
   * \param[in] d: data of the computation by kernels
   * \param[in] c: component
   */
  template <bool parallel>
  [[nodiscard]] static bool computeL2ProjectionRightHandSideByKernels(
      Context& ctx,
      mfem::Vector& b,
      const FiniteElementSpace<parallel>& fes,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const L2ProjectionKernelsData& d,
      const size_type c) {
    const auto* const r =
        fes.GetElementRestriction(mfem::ElementDofOrdering::NATIVE);
    auto e = mfem::Vector(r->Height());
    e.UseDevice(true);
    e = real{0};
    for (size_type i = 0; i != fcts.size(); ++i) {
      const auto& f = fcts[i];
      const auto& qs = f.getPartialQuadratureSpace();
      if (qs.getNumberOfElements() == 0) {
        continue;
      }
      const auto& g = d.geometries.at(qs.getUniformIntegrationRule());
      if (!batch::addL2ProjectionRightHandSide(
              ctx, e, qs, d.values[i], f.getNumberOfComponents(), c, g)) {
        return false;
      }
    }
    b.UseDevice(true);
    r->MultTranspose(e, b);
    return true;
  }  // end of computeL2ProjectionRightHandSideByKernels

  /*!
   * \brief integrator of the mass matrices of the elements computed by
   * kernels (see `batch::computeMassMatrices`), for the assembly of MFEM on
   * the device (`AssemblyLevel::FULL`), in the ordering of the element
   * matrices of MFEM (`mfem::GetEVectorOrdering`). The integration rule is
   * the one of `mfem::MassIntegrator`. Unlike the assembly of
   * `mfem::MassIntegrator` on the device, it supports the simplices, and it
   * does not give nodes to a mesh without nodes (the geometric factors of
   * MFEM are not used).
   */
  struct MassMatricesKernelsIntegrator final
      : public mfem::BilinearFormIntegrator {
    void AssembleEA(const mfem::FiniteElementSpace& fes,
                    mfem::Vector& m,
                    const bool add) override {
      if (fes.GetNE() == 0) {
        return;
      }
      const auto& fe = *(fes.GetFE(0));
      const auto& ir = mfem::MassIntegrator::GetRule(
          fe, fe, *(fes.GetElementTransformation(0)));
      auto ctx = Context{};
      auto g = batch::ElementsGeometry{};
      if ((!batch::computeElementsGeometry(ctx, g, fes, ir)) ||
          (!batch::computeMassMatrices(ctx, m, g, mfem::GetEVectorOrdering(fes),
                                       add))) {
        raise(ctx.getErrorMessage());
      }
    }  // end of AssembleEA
  };

  /*!
   * \brief computes the regularization of functions defined at the
   * integration points: L2 projection, operator `M + d K`, where `M` is the
   * mass matrix and `K` the matrix of the diffusion, if `d` is not null
   * (implicit gradient regularization)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] r: result
   * \param[in, out] l: linear solver
   * \param[in] fcts: functions
   * \param[in] d: coefficient of the diffusion
   */
  template <bool parallel>
  static bool updatePartialQuadratureRegularization_impl(
      Context& ctx,
      L2ProjectionResult<parallel>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const real d) {
    if (!checkUpdateFunctionsArguments(ctx, r, l, fcts)) {
      return false;
    }
    //
    *(r.result) = real{};
    //
    // Thus we introduce:
    // -  local_fespace and local_gridfunction which are only allocated
    // -  fespace and x which points with the finite element space and grid
    //    function used for the resolution
    //
    const auto fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& rfespace = *(
        static_cast<const FiniteElementSpace<parallel>*>(r.result->FESpace()));
    const auto* const mesh = [&rfespace] {
      if constexpr (parallel) {
        return rfespace.GetParMesh();
      } else {
        return rfespace.GetMesh();
      }
    }();
    auto fespaces_manager = fed.getFiniteElementSpacesManager();
    auto olocal_fespace =
        fespaces_manager.template getFiniteElementSpace<parallel>(ctx, *mesh,
                                                                  1);
    if (isInvalid(olocal_fespace)) {
      return {};
    }
    // partial quadrature spaces defined on boundaries are built on the submesh
    const auto on_boundaries = fed.isDefinedOnBoundaries(ctx, *mesh);
    if (isInvalid(on_boundaries)) {
      return {};
    }
    const auto use_parent_elements =
        (r.submesh != nullptr) && (!*on_boundaries);
    // Here begin the tricky part. In order to reduce the
    // computational cost, we will make the projection component by
    // component. But of course, if the functions are scalar, we do
    // want to avoid creating a temporary grid function and copying it
    // in *(r.result).
    //
    auto local_gridfunction = std::unique_ptr<GridFunction<parallel>>{};
    auto* const x = [&r, &rfespace, &local_gridfunction, &olocal_fespace] {
      if (rfespace.GetVDim() == 1) {
        return r.result.get();
      }
      local_gridfunction =
          std::make_unique<GridFunction<parallel>>(&(*olocal_fespace));
      return local_gridfunction.get();
    }();
    // computation by kernels: right-hand sides by kernels, operator
    // assembled by MFEM, on the device if a device is enabled
    const auto kernels_data = [&] {
      CatchTimeSection(ctx, "L2Projection::KernelsData");
      return makeL2ProjectionKernelsData<parallel>(may_throw, *olocal_fespace,
                                                   fcts, r.submesh != nullptr);
    }();
    // Regularization operator
    auto operator_section = std::optional<mgis::ProfilingSection>(
        std::in_place, ctx, "L2Projection::Operator", ctx.isProfilingEnabled());
    BilinearForm<parallel> a(&(*olocal_fespace));
    auto diffusion = mfem::ConstantCoefficient(d);
    // computation by kernels of the L2 projection: mass matrices computed by
    // kernels and assembled by MFEM, on the device if a device is enabled;
    // the operator of the implicit gradient regularization is assembled
    // element by element
    if ((kernels_data.has_value()) && (d == 0)) {
      a.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
      a.AddDomainIntegrator(new MassMatricesKernelsIntegrator);
    } else {
      a.AddDomainIntegrator(new mfem::MassIntegrator);
    }
    if (d != 0) {
      a.AddDomainIntegrator(new mfem::DiffusionIntegrator(diffusion));
    }
    a.Assemble();
    operator_section.reset();
    // resolution(s)
    for (size_type c = 0; c != rfespace.GetVDim(); ++c) {
      // initialize the solution, if needed
      if (rfespace.GetVDim() != 1) {
        *x = real{};
      }
      // right-hand side
      auto rhs_section = std::optional<mgis::ProfilingSection>(
          std::in_place, ctx, "L2Projection::RightHandSide",
          ctx.isProfilingEnabled());
      LinearForm<parallel> b(&(*olocal_fespace));
      if (kernels_data.has_value()) {
        if (!computeL2ProjectionRightHandSideByKernels<parallel>(
                ctx, b, *olocal_fespace, fcts, *kernels_data, c)) {
          return false;
        }
      } else {
        if (rfespace.GetVDim() == 1) {
          if (use_parent_elements) {
            b.AddDomainIntegrator(new ScalarL2ProjectionRHSFormIntegratorII(
                r.submesh->GetParentElementIDMap(), fcts));
          } else {
            b.AddDomainIntegrator(
                new ScalarL2ProjectionRHSFormIntegrator(fcts));
          }
        } else {
          if (use_parent_elements) {
            b.AddDomainIntegrator(new ComponentL2ProjectionRHSFormIntegratorII(
                r.submesh->GetParentElementIDMap(), fcts, c));
          } else {
            b.AddDomainIntegrator(
                new ComponentL2ProjectionRHSFormIntegrator(fcts, c));
          }
        }
        b.Assemble();
      }
      rhs_section.reset();
      // resolution
      CatchTimeSection(ctx, "L2Projection::Solve");
      mfem::Array<int> boundary_dofs;
      // the operator handle owns the matrix built by the full assembly
      auto A = mfem::OperatorHandle{};
      mfem::Vector B, X;
      a.FormLinearSystem(boundary_dofs, *x, b, A, X, B);
      l.linear_solver->SetOperator(*A);
      l.linear_solver->Mult(B, X);
      clearHypreErrors();
      auto* const isolver =
          dynamic_cast<IterativeSolver*>(l.linear_solver.get());
      if (isolver != nullptr) {
        if (!isolver->GetConverged()) {
          return ctx.registerErrorMessage("linear solver failed");
        }
      }
      a.RecoverFEMSolution(X, b, *x);
      // copy solution to r.result, if needed. The copy is done on the host,
      // where the solution, which may have been computed on the device, is
      // made valid, as the components already copied
      if (rfespace.GetVDim() != 1) {
        auto* const dest = r.result->HostReadWrite();
        const auto* const src = x->HostRead();
        const auto bynodes =
            (rfespace.GetOrdering() == mfem::Ordering::byNODES);
        if (bynodes) {
          const auto n = x->Size();
          const auto offset = x->Size() * c;
          for (size_type idx = 0; idx != n; ++idx) {
            dest[offset + idx] = src[idx];
          }
        } else {
          const auto n = rfespace.GetVDim();
          for (size_type idx = 0; idx != x->Size(); ++idx) {
            dest[idx * n + c] = src[idx];
          }
        }
      }
    }
    if (kernels_data.has_value()) {
      internals::countPostProcessingComputationByKernels();
    }
    return true;
  }  // end of updatePartialQuadratureRegularization_impl

  /*!
   * \brief calls `updatePartialQuadratureRegularization_impl`, the
   * exceptions being reported in the execution context
   */
  template <bool parallel>
  static bool updatePartialQuadratureRegularization(
      Context& ctx,
      L2ProjectionResult<parallel>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const real d) noexcept {
    try {
      return updatePartialQuadratureRegularization_impl<parallel>(ctx, r, l,
                                                                  fcts, d);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "updatePartialQuadratureRegularization: unknown exception");
    }
  }  // end of updatePartialQuadratureRegularization

  template <bool parallel>
  static bool updateL2Projection_impl(
      Context& ctx,
      L2ProjectionResult<parallel>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    return updatePartialQuadratureRegularization<parallel>(ctx, r, l, fcts, 0);
  }  // end of updateL2Projection_impl

  template <>
  bool updateL2Projection<true>(
      Context& ctx,
      L2ProjectionResult<true>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
#ifdef MFEM_USE_MPI
    return updateL2Projection_impl<true>(ctx, r, l, fcts);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of updateL2Projection

  template <>
  bool updateL2Projection<false>(
      Context& ctx,
      L2ProjectionResult<false>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    return updateL2Projection_impl<false>(ctx, r, l, fcts);
  }  // end of updateL2Projection

  template <bool parallel>
  static std::optional<L2ProjectionResult<parallel>> computeL2Projection_impl(
      Context& ctx,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    auto ores = createL2ProjectionResult<parallel>(ctx, fcts);
    if (isInvalid(ores)) {
      return {};
    }
    const auto ok = updateL2Projection<parallel>(ctx, *ores, l, fcts);
    if (!ok) {
      return {};
    }
    return ores;
  }  // end of computeL2Projection_impl

  template <>
  std::optional<L2ProjectionResult<true>> computeL2Projection<true>(
      Context& ctx,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
#ifdef MFEM_USE_MPI
    return computeL2Projection_impl<true>(ctx, l, fcts);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of computeL2Projection

  template <>
  std::optional<L2ProjectionResult<false>> computeL2Projection<false>(
      Context& ctx,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    return computeL2Projection_impl<false>(ctx, l, fcts);
  }  // end of computeL2Projection

  template <bool parallel>
  [[nodiscard]] static bool updateImplicitGradientRegularization_impl(
      Context& ctx,
      ImplicitGradientRegularizationResult<parallel>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const real lc) noexcept {
    return updatePartialQuadratureRegularization<parallel>(ctx, r, l, fcts,
                                                           lc * lc);
  }  // end of updateImplicitGradientRegularization_impl

  template <>
  bool updateImplicitGradientRegularization<true>(
      Context& ctx,
      ImplicitGradientRegularizationResult<true>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const real lc) noexcept {
#ifdef MFEM_USE_MPI
    return updateImplicitGradientRegularization_impl<true>(ctx, r, l, fcts, lc);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of updateImplicitGradientRegularization

  template <>
  bool updateImplicitGradientRegularization<false>(
      Context& ctx,
      ImplicitGradientRegularizationResult<false>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const real lc) noexcept {
    return updateImplicitGradientRegularization_impl<false>(ctx, r, l, fcts,
                                                            lc);
  }  // end of updateImplicitGradientRegularization

  template <bool parallel>
  static std::optional<ImplicitGradientRegularizationResult<parallel>>
  computeImplicitGradientRegularization_impl(
      Context& ctx,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const real lc) noexcept {
    auto ores = createL2ProjectionResult<parallel>(ctx, fcts);
    if (isInvalid(ores)) {
      return {};
    }
    const auto ok =
        updateImplicitGradientRegularization<parallel>(ctx, *ores, l, fcts, lc);
    if (!ok) {
      return {};
    }
    return ores;
  }  // end of computeImplicitGradientRegularization_impl

  template <>
  std::optional<ImplicitGradientRegularizationResult<true>>
  computeImplicitGradientRegularization<true>(
      Context& ctx,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const real lc) noexcept {
#ifdef MFEM_USE_MPI
    return computeImplicitGradientRegularization_impl<true>(ctx, l, fcts, lc);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of computeImplicitGradientRegularization

  template <>
  std::optional<ImplicitGradientRegularizationResult<false>>
  computeImplicitGradientRegularization<false>(
      Context& ctx,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const real lc) noexcept {
    return computeImplicitGradientRegularization_impl<false>(ctx, l, fcts, lc);
  }  // end of computeImplicitGradientRegularization

}  // end of namespace mfem_mgis
