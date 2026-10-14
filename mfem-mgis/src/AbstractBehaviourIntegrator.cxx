/*!
 * \file   src/AbstractBehaviourIntegrator.cxx
 * \brief  This file implements the `AbstractBehaviourIntegrator` class and
 * the `computeMeasure` and `computeIntegral` functions
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#include <cmath>
#include <optional>
#include <cstdlib>
#include <algorithm>
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"
#include "MFEMMGIS/Batch/IntegrationPointsFunctions.hxx"
#include "PostProcessingsComputation.hxx"

namespace mfem_mgis {

  AbstractBehaviourIntegrator::~AbstractBehaviourIntegrator() = default;

  bool AbstractBehaviourIntegrator::supportsBatchIntegration() const noexcept {
    return false;
  }  // end of supportsBatchIntegration

  bool AbstractBehaviourIntegrator::computeGradients(
      const mfem::FiniteElement &,
      mfem::ElementTransformation &,
      const mfem::Vector &) {
    return false;
  }  // end of computeGradients

  bool AbstractBehaviourIntegrator::integrateBatch(const IntegrationType) {
    return false;
  }  // end of integrateBatch

  bool AbstractBehaviourIntegrator::integratesBatchesOnDevice() const noexcept {
    return false;
  }  // end of integratesBatchesOnDevice

  bool AbstractBehaviourIntegrator::supportsKernels() const noexcept {
    return false;
  }  // end of supportsKernels

  bool AbstractBehaviourIntegrator::requiresUnknownsValuesForKernels()
      const noexcept {
    return false;
  }  // end of requiresUnknownsValuesForKernels

  std::vector<const mfem::IntegrationRule*>
  AbstractBehaviourIntegrator::getAdditionalIntegrationRulesForKernels() const {
    return {};
  }  // end of getAdditionalIntegrationRulesForKernels

  bool AbstractBehaviourIntegrator::computeGradientsByKernels(
      Context& ctx,
      const batch::UnknownsDerivativesSet&,
      const batch::ElementsGeometrySet&) noexcept {
    return ctx.registerErrorMessage(
        "AbstractBehaviourIntegrator::computeGradientsByKernels: the "
        "gradients of this behaviour integrator can't be computed by "
        "kernels");
  }  // end of computeGradientsByKernels

  bool AbstractBehaviourIntegrator::addInnerForcesByKernels(
      Context& ctx, mfem::Vector&, const batch::ElementsGeometrySet&) noexcept {
    return ctx.registerErrorMessage(
        "AbstractBehaviourIntegrator::addInnerForcesByKernels: the inner "
        "forces of this behaviour integrator can't be computed by kernels");
  }  // end of addInnerForcesByKernels

  bool AbstractBehaviourIntegrator::addStiffnessMatricesByKernels(
      Context& ctx, mfem::Vector&, const batch::ElementsGeometrySet&) noexcept {
    return ctx.registerErrorMessage(
        "AbstractBehaviourIntegrator::addStiffnessMatricesByKernels: the "
        "stiffness matrices of this behaviour integrator can't be computed by "
        "kernels");
  }  // end of addStiffnessMatricesByKernels

  bool AbstractBehaviourIntegrator::addTangentOperatorActionByKernels(
      Context& ctx,
      mfem::Vector&,
      const mfem::Vector&,
      const batch::ElementsGeometrySet&) noexcept {
    return ctx.registerErrorMessage(
        "AbstractBehaviourIntegrator::addTangentOperatorActionByKernels: the "
        "action of the tangent operator of this behaviour integrator can't be "
        "computed by kernels");
  }  // end of addTangentOperatorActionByKernels

  bool AbstractBehaviourIntegrator::addTangentOperatorDiagonalByKernels(
      Context& ctx, mfem::Vector&, const batch::ElementsGeometrySet&) noexcept {
    return ctx.registerErrorMessage(
        "AbstractBehaviourIntegrator::addTangentOperatorDiagonalByKernels: the "
        "diagonal of the tangent operator of this behaviour integrator can't "
        "be computed by kernels");
  }  // end of addTangentOperatorDiagonalByKernels

  bool AbstractBehaviourIntegrator::addStiffnessMatricesByKernels(
      Context& ctx,
      const batch::AssembledStiffnessMatrix&,
      const batch::ElementsGeometrySet&) noexcept {
    return ctx.registerErrorMessage(
        "AbstractBehaviourIntegrator::addStiffnessMatricesByKernels: the "
        "stiffness matrices of this behaviour integrator can't be computed by "
        "kernels");
  }  // end of addStiffnessMatricesByKernels

  template <bool parallel, typename Functor>
  static void performLoopOverIntegrationPoints(
      Functor& f, const AbstractBehaviourIntegrator& bi) {
    const auto& qspace = bi.getPartialQuadratureSpace();
    const auto& fespace = qspace.getFiniteElementDiscretization()
                              .getFiniteElementSpace<parallel>();
    for (const auto [e, eo] : qspace.getOffsets()) {
      const auto& fe = *(fespace.GetFE(e));
      auto& tr = *(fespace.GetElementTransformation(e));
      const auto& ir = bi.getIntegrationRule(fe, tr);
      for (size_type i = 0; i != ir.GetNPoints(); ++i) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        f(fe, tr, ip, eo + i);
      }
    }
  }  // end of performLoopOverIntegrationPoints

  /*!
   * \return the integral of a scalar function, or the measure if no function
   * is given, over the elements of the quadrature space of a behaviour
   * integrator owned by this process, computed by kernels with the weight
   * `w det J` of the behaviour integrators of mfem-mgis (zero if this
   * process owns no element of the space), or an empty value if the
   * post-processings are not computed by kernels (see
   * `internals::usePostProcessingsKernels`), if the behaviour integrator
   * does not support kernels, or if the geometric data of the elements
   * can't be computed by kernels (elements of different geometries). Throws
   * if the values of the function can't be gathered or integrated.
   * \param[in] bi: behaviour integrator
   * \param[in] f: function, or null
   */
  template <bool parallel>
  static std::optional<real> computeLocalIntegralByKernels(
      attributes::MayThrow,
      const AbstractBehaviourIntegrator& bi,
      const ImmutablePartialQuadratureFunctionView* const f) {
    if ((!internals::usePostProcessingsKernels(may_throw)) ||
        (!bi.supportsKernels())) {
      return {};
    }
    const auto& qspace = bi.getPartialQuadratureSpace();
    const auto* const ir = qspace.getUniformIntegrationRule();
    if (ir == nullptr) {
      // no element owned by this process
      internals::countPostProcessingComputationByKernels();
      return real{0};
    }
    const auto& fespace = qspace.getFiniteElementDiscretization()
                              .getFiniteElementSpace<parallel>();
    auto ctx = Context{};
    // geometric data of the elements of the space only
    auto g = batch::ElementsGeometry{};
    if (!batch::computeElementsGeometry(ctx, g, fespace, *ir,
                                        qspace.getElements())) {
      return {};
    }
    auto values = mfem::Vector{};
    if ((f != nullptr) && (!batch::gatherValues(ctx, values, *f))) {
      raise(ctx.getErrorMessage());
    }
    const auto r = batch::computeIntegral(ctx, qspace,
                                          f != nullptr ? &values : nullptr, g);
    if (!r.has_value()) {
      raise(ctx.getErrorMessage());
    }
    internals::countPostProcessingComputationByKernels();
    return r;
  }  // end of computeLocalIntegralByKernels

  template <bool parallel>
  static real BehaviourIntegrator_computeMeasure(
      const AbstractBehaviourIntegrator& bi) {
    const auto rk =
        computeLocalIntegralByKernels<parallel>(may_throw, bi, nullptr);
    if (rk.has_value()) {
      return *rk;
    }
    // contributions accumulated in double, as by the kernels
    auto r = double{};
    auto integrate =
        [&r, &bi](const mfem::FiniteElement&, mfem::ElementTransformation& tr,
                  const mfem::IntegrationPoint& ip, const size_type) {
          r += bi.getIntegrationPointWeight(tr, ip);
        };
    performLoopOverIntegrationPoints<parallel>(integrate, bi);
    return static_cast<real>(r);
  }  // end of BehaviourIntegrator_computeMeasure

  real computeMeasure(const AbstractBehaviourIntegrator& bi) {
    const auto& qspace = bi.getPartialQuadratureSpace();
    const auto& fed = qspace.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto lv = BehaviourIntegrator_computeMeasure<true>(bi);
      auto r = real{};
      MPI_Allreduce(&lv, &r, 1, mpi_type<real>, MPI_SUM,
                    getMPICommunicator(fed));
      return r;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return BehaviourIntegrator_computeMeasure<false>(bi);
  }  // end of computeIntegral

  template <bool parallel>
  static real BehaviourIntegrator_computeScalarIntegral(
      const AbstractBehaviourIntegrator& bi,
      const ImmutablePartialQuadratureFunctionView& f) {
    const auto rk = computeLocalIntegralByKernels<parallel>(may_throw, bi, &f);
    if (rk.has_value()) {
      return *rk;
    }
    const auto* const values = f.getValues().data() + f.getDataOffset();
    // contributions accumulated in double, as by the kernels
    auto r = double{};
    if (f.getDataStride() == 1) {
      auto integrate = [&r, &bi, values](const mfem::FiniteElement&,
                                         mfem::ElementTransformation& tr,
                                         const mfem::IntegrationPoint& ip,
                                         const size_type o) {
        const auto v = values[o];
        r += v * bi.getIntegrationPointWeight(tr, ip);
      };
      performLoopOverIntegrationPoints<parallel>(integrate, bi);
    } else {
      const auto s = f.getDataStride();
      auto integrate = [&r, &bi, &s, values](const mfem::FiniteElement&,
                                             mfem::ElementTransformation& tr,
                                             const mfem::IntegrationPoint& ip,
                                             const size_type o) {
        const auto v = values[o * s];
        r += v * bi.getIntegrationPointWeight(tr, ip);
      };
      performLoopOverIntegrationPoints<parallel>(integrate, bi);
    }
    return static_cast<real>(r);
  }  // end of BehaviourIntegrator_computeScalarIntegral

  template <>
  real computeIntegral(const AbstractBehaviourIntegrator& bi,
                       const ImmutablePartialQuadratureFunctionView& f) {
    const auto& qspace = f.getPartialQuadratureSpace();
    if (&qspace != &(bi.getPartialQuadratureSpace())) {
      raise(
          "computeIntegral: the partial quadrature space of the behaviour "
          "integrator don't match the one of the partial quadrature function");
    }
    if (f.getNumberOfComponents() != 1u) {
      raise(
          "computeIntegral: the partial quadrature function is not scalar "
          "valuated");
    }
    const auto& fed = qspace.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto lv = BehaviourIntegrator_computeScalarIntegral<true>(bi, f);
      auto r = real{};
      MPI_Allreduce(&lv, &r, 1, mpi_type<real>, MPI_SUM,
                    getMPICommunicator(fed));
      return r;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return BehaviourIntegrator_computeScalarIntegral<false>(bi, f);
  }  // end of computeIntegral

}  // end of namespace mfem_mgis
