/*!
 * \file   src/NonLinearEvolutionProblemImplementation.cxx
 * \brief  This file implements the `NonLinearEvolutionProblemImplementation`
 * class
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include <tuple>
#include <optional>
#include <algorithm>
#include <string>
#include "mfem/fem/datacollection.hpp"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/fem/linearform.hpp"
#include "mfem/fem/bilinearform.hpp"
#include "mfem/fem/bilininteg.hpp"
#include "mfem/fem/lininteg.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/pbilinearform.hpp"
#endif /* MFEM_USE_MPI */

#include "MGIS/Raise.hxx"
#include "MGIS/Profiling.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/PostProcessingFactory.hxx"
#include "MFEMMGIS/AbstractBoundaryCondition.hxx"
#include "MFEMMGIS/AbstractDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/NonLinearSolvers/NewtonSolver.hxx"
#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverFactory.hxx"
#include "MFEMMGIS/PostProcessing/NonLinearEvolutionProblemPostProcessingBase.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  /*!
   * \return true if all the behaviours are integrated on the device, without
   * gradients computed on the host: the arrays of the materials do not need
   * to be synchronized on the host
   * \param[in] mi: multi-material integrator
   * \param[in] batch: if true, the behaviours are integrated by batches
   * \param[in] kernels: if true, the gradients are computed by kernels
   * \param[in] it: integration type
   */
  static bool integratesOnDevice(const MultiMaterialNonLinearIntegrator& mi,
                                 const bool batch,
                                 const bool kernels,
                                 const IntegrationType it) noexcept {
    const auto prediction =
        (it == IntegrationType::PREDICTION_TANGENT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_SECANT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_ELASTIC_OPERATOR);
    return batch && (kernels || prediction) && mi.integratesBatchesOnDevice();
  }  // end of integratesOnDevice

  /*!
   * \return if the given multi-material integrator is the only integrator of
   * a non linear form
   * \param[in] f: non linear form
   * \param[in] mi: multi-material integrator
   */
  static bool isTheOnlyIntegrator(
      const mfem::NonlinearForm& f,
      const MultiMaterialNonLinearIntegrator* const mi) noexcept {
    const auto& dnfi = *(f.GetDNFI());
    return (mi != nullptr) && (dnfi.Size() == 1) && (dnfi[0] == mi) &&
           (f.GetBNFI()->Size() == 0) &&
           (f.GetInteriorFaceIntegrators().Size() == 0) &&
           (f.GetBdrFaceIntegrators().Size() == 0);
  }  // end of isTheOnlyIntegrator

  namespace {

    /*!
     * \brief gives access to the markers of the integrators of a non linear
     * form and to its extension, protected members of
     * `mfem::NonlinearForm`: the pointers to these members are formed in
     * this derived class, which is never instantiated, and applied to any
     * non linear form
     */
    struct NonlinearFormMembers : mfem::NonlinearForm {
      //! \return the markers of the domain integrators
      static const mfem::Array<mfem::Array<int>*>& getDomainMarkers(
          const mfem::NonlinearForm& f) noexcept {
        return f.*(&NonlinearFormMembers::dnfi_marker);
      }
      //! \return the markers of the boundary integrators
      static const mfem::Array<mfem::Array<int>*>& getBoundaryMarkers(
          const mfem::NonlinearForm& f) noexcept {
        return f.*(&NonlinearFormMembers::bnfi_marker);
      }
      //! \return the markers of the boundary face integrators
      static const mfem::Array<mfem::Array<int>*>& getBoundaryFaceMarkers(
          const mfem::NonlinearForm& f) noexcept {
        return f.*(&NonlinearFormMembers::bfnfi_marker);
      }
      //! \return true if the form has an extension (partial assembly)
      static bool hasExtension(const mfem::NonlinearForm& f) noexcept {
        return f.*(&NonlinearFormMembers::ext) != nullptr;
      }
    };

  }  // end of namespace

  /*!
   * \return true if the integrators of a non linear form other than the
   * multi-material integrator can be evaluated by
   * `addOtherIntegratorsContributions`, the
   * multi-material integrator being one of the domain integrators: the form
   * has no interior face integrator nor extension
   * \param[in] f: non linear form
   * \param[in] mi: multi-material integrator
   */
  static bool supportsOtherIntegrators(
      const mfem::NonlinearForm& f,
      const MultiMaterialNonLinearIntegrator* const mi) noexcept {
    const auto& dnfi = *(f.GetDNFI());
    return (mi != nullptr) &&
           (std::find(dnfi.begin(), dnfi.end(), mi) != dnfi.end()) &&
           (f.GetInteriorFaceIntegrators().Size() == 0) &&
           (!NonlinearFormMembers::hasExtension(f));
  }  // end of supportsOtherIntegrators

  /*!
   * \return the attributes of the elements (or of the boundary elements) on
   * which at least one of the given integrators, other than the excluded
   * one, is evaluated, as by the loops of `mfem::NonlinearForm`
   * \param[in] integrators: integrators
   * \param[in] markers: markers of the integrators (null for all the
   * attributes)
   * \param[in] n: number of attributes
   * \param[in] excluded: integrator not evaluated
   * \throws std::runtime_error if a marker has less values than the number
   * of attributes
   */
  static mfem::Array<int> getEvaluatedAttributes(
      const mfem::Array<mfem::NonlinearFormIntegrator*>& integrators,
      const mfem::Array<mfem::Array<int>*>& markers,
      const int n,
      const mfem::NonlinearFormIntegrator* const excluded) {
    auto attributes = mfem::Array<int>(n);
    attributes = 0;
    auto all = false;
    for (int k = 0; k != integrators.Size(); ++k) {
      if (integrators[k] == excluded) {
        continue;
      }
      if (markers[k] == nullptr) {
        all = true;
        continue;
      }
      if (markers[k]->Size() < n) {
        raise("invalid marker of an integrator of the non linear form: " +
              std::to_string(markers[k]->Size()) + " attributes instead of " +
              std::to_string(n));
      }
      for (int i = 0; i != n; ++i) {
        attributes[i] |= (*markers[k])[i];
      }
    }
    if (all) {
      attributes = 1;
    }
    return attributes;
  }  // end of getEvaluatedAttributes

  /*!
   * \brief adds the contributions of the integrators of a non linear form
   * other than the multi-material integrator (domain, boundary and boundary
   * face integrators) to a residual or to a matrix, by the loops of
   * `mfem::NonlinearForm::Mult` and `mfem::NonlinearForm::GetGradient` (see
   * `supportsOtherIntegrators`), on the host
   * \tparam Residual: if true, the contributions to a residual are added,
   * otherwise the contributions to a matrix
   * \param[in] f: non linear form
   * \param[in] mi: multi-material integrator, not evaluated
   * \param[in] u: unknowns (local vector of the finite element space)
   * \param[in] add: function adding the vector (residual) or the matrix of an
   * element, given the degrees of freedom
   */
  template <bool Residual, typename AddFunction>
  static void addOtherIntegratorsContributions(
      const mfem::NonlinearForm& f,
      const MultiMaterialNonLinearIntegrator* const mi,
      const mfem::Vector& u,
      const AddFunction& add) {
    const auto& fes = *(f.FESpace());
    const auto& mesh = *(fes.GetMesh());
    const auto& dnfi = *(f.GetDNFI());
    const auto& bnfi = *(f.GetBNFI());
    const auto& bfnfi = f.GetBdrFaceIntegrators();
    const auto& dmarkers = NonlinearFormMembers::getDomainMarkers(f);
    const auto& bmarkers = NonlinearFormMembers::getBoundaryMarkers(f);
    const auto& bfmarkers = NonlinearFormMembers::getBoundaryFaceMarkers(f);
    const auto* const excluded =
        static_cast<const mfem::NonlinearFormIntegrator*>(mi);
    const auto nattributes =
        mesh.attributes.Size() != 0 ? mesh.attributes.Max() : 0;
    const auto nbattributes =
        mesh.bdr_attributes.Size() != 0 ? mesh.bdr_attributes.Max() : 0;
    auto vdofs = mfem::Array<int>{};
    auto ue = mfem::Vector{};
    auto ve = mfem::Vector{};
    auto me = mfem::DenseMatrix{};
    auto doftrans = mfem::DofTransformation{};
    // evaluates an integrator on an element and adds its contribution
    const auto evaluate = [&](mfem::NonlinearFormIntegrator& i,
                              const mfem::FiniteElement& e,
                              mfem::ElementTransformation& tr) {
      if constexpr (Residual) {
        i.AssembleElementVector(e, tr, ue, ve);
        doftrans.TransformDual(ve);
        add(vdofs, ve);
      } else {
        i.AssembleElementGrad(e, tr, ue, me);
        doftrans.TransformDual(me);
        add(vdofs, me);
      }
    };
    // domain integrators
    const auto domain =
        getEvaluatedAttributes(dnfi, dmarkers, nattributes, excluded);
    if (std::find(domain.begin(), domain.end(), 1) != domain.end()) {
      for (int i = 0; i != fes.GetNE(); ++i) {
        const auto attribute = mesh.GetAttribute(i);
        if (domain[attribute - 1] == 0) {
          continue;
        }
        const auto& e = *(fes.GetFE(i));
        fes.GetElementVDofs(i, vdofs, doftrans);
        auto& tr = *(fes.GetElementTransformation(i));
        u.GetSubVector(vdofs, ue);
        doftrans.InvTransformPrimal(ue);
        for (int k = 0; k != dnfi.Size(); ++k) {
          if ((dnfi[k] == excluded) || ((dmarkers[k] != nullptr) &&
                                        ((*dmarkers[k])[attribute - 1] == 0))) {
            continue;
          }
          evaluate(*dnfi[k], e, tr);
        }
      }
    }
    // boundary integrators
    const auto boundary =
        getEvaluatedAttributes(bnfi, bmarkers, nbattributes, excluded);
    if (std::find(boundary.begin(), boundary.end(), 1) != boundary.end()) {
      for (int i = 0; i != fes.GetNBE(); ++i) {
        const auto attribute = mesh.GetBdrAttribute(i);
        if (boundary[attribute - 1] == 0) {
          continue;
        }
        const auto& e = *(fes.GetBE(i));
        fes.GetBdrElementVDofs(i, vdofs, doftrans);
        auto& tr = *(fes.GetBdrElementTransformation(i));
        u.GetSubVector(vdofs, ue);
        doftrans.InvTransformPrimal(ue);
        for (int k = 0; k != bnfi.Size(); ++k) {
          if ((bmarkers[k] != nullptr) &&
              ((*bmarkers[k])[attribute - 1] == 0)) {
            continue;
          }
          evaluate(*bnfi[k], e, tr);
        }
      }
    }
    // boundary face integrators, evaluated with the degrees of freedom of
    // the element of the face
    const auto faces =
        getEvaluatedAttributes(bfnfi, bfmarkers, nbattributes, excluded);
    if (std::find(faces.begin(), faces.end(), 1) != faces.end()) {
      auto tr = mfem::FaceElementTransformations{};
      auto tr1 = mfem::IsoparametricTransformation{};
      auto tr2 = mfem::IsoparametricTransformation{};
      for (int i = 0; i != fes.GetNBE(); ++i) {
        const auto attribute = mesh.GetBdrAttribute(i);
        if (faces[attribute - 1] == 0) {
          continue;
        }
        // interior, shared or non conforming faces are skipped, as by MFEM
        mesh.GetBdrFaceTransformations(i, tr, tr1, tr2);
        if (tr.GetGeometryType() == mfem::Geometry::INVALID) {
          continue;
        }
        fes.GetElementVDofs(tr.Elem1No, vdofs);
        u.GetSubVector(vdofs, ue);
        const auto& e = *(fes.GetFE(tr.Elem1No));
        for (int k = 0; k != bfnfi.Size(); ++k) {
          if ((bfmarkers[k] != nullptr) &&
              ((*bfmarkers[k])[attribute - 1] == 0)) {
            continue;
          }
          if constexpr (Residual) {
            bfnfi[k]->AssembleFaceVector(e, e, tr, ue, ve);
            add(vdofs, ve);
          } else {
            bfnfi[k]->AssembleFaceGrad(e, e, tr, ue, me);
            add(vdofs, me);
          }
        }
      }
    }
  }  // end of addOtherIntegratorsContributions

  /*!
   * \brief adds the contributions of the integrators of a non linear form
   * other than the multi-material integrator to the local matrix assembled
   * by kernels (their entries are in the pattern of the elements), on the
   * host. A failure on one process is reported on all of them.
   * \param[in] fed: finite element discretization
   * \param[in] f: non linear form
   * \param[in] mi: multi-material integrator, not evaluated
   * \param[in] u: unknowns (local vector of the finite element space),
   * prolongated by the caller on all the processes
   * \param[in, out] K: local matrix
   */
  static void addOtherIntegratorsContributionsToTheLocalMatrix(
      const FiniteElementDiscretization& fed,
      const mfem::NonlinearForm& f,
      const MultiMaterialNonLinearIntegrator* const mi,
      const mfem::Vector& u,
      mfem::SparseMatrix& K) {
    auto ctx = Context{};
    auto success = false;
    try {
      static_cast<void>(K.HostReadWriteData());
      addOtherIntegratorsContributions<false>(
          f, mi, u,
          [&K](const mfem::Array<int>& dofs, const mfem::DenseMatrix& m) {
            K.AddSubMatrix(dofs, dofs, m, 0);
          });
      success = true;
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    raiseOnFailureOnAnyProcess(fed, ctx, success,
                               "computation of the tangent operator");
  }  // end of addOtherIntegratorsContributionsToTheLocalMatrix

  /*!
   * \brief checks that the tangent operator can be applied without assembly
   * (see `TangentOperatorComputation`): the non linear form must only
   * contain the integrator of mfem-mgis, the contributions of the other
   * integrators to the tangent operator being otherwise lost
   * \param[in] f: non linear form
   * \param[in] mi: multi-material integrator
   */
  static void checkMatrixFreeTangentOperatorSupport(
      const mfem::NonlinearForm& f,
      const MultiMaterialNonLinearIntegrator* const mi) {
    // this check gives the same result on all the processes
    if (!isTheOnlyIntegrator(f, mi)) {
      throw(CollectiveFailure(
          "the tangent operator applied without assembly (parameter "
          "TangentOperatorComputation set to MatrixFree) requires that the "
          "non linear form only contains the integrator of mfem-mgis (with "
          "the inner forces computed by kernels, the boundary conditions "
          "independent of the unknowns, such as an imposed pressure, are not "
          "integrators of the non linear form)"));
    }
  }  // end of checkMatrixFreeTangentOperatorSupport

  /*!
   * \brief adds the inner forces computed by kernels and the contributions
   * of the boundary conditions independent of the unknowns to the residual
   * (true degrees of freedom), and cancels the residual at the essential
   * degrees of freedom
   * \param[in] fed: finite element discretization
   * \param[in, out] mi: multi-material integrator, null if the inner forces
   * are not computed by kernels
   * \param[in] bcs: boundary conditions independent of the unknowns
   * \param[in] fes: finite element space
   * \param[in] essential_dofs: essential true degrees of freedom
   * \param[out] f: inner forces, contributions of the boundary conditions
   * and of the other integrators of the non linear form, if any (local
   * vector of the finite element space)
   * \param[in, out] r: residual
   * \param[in] addOtherContributions: function adding the contributions of
   * the other integrators of the non linear form to `f`, if not empty, on
   * the host, called if the previous computations succeeded
   */
  static void addInnerForcesComputedByKernels(
      const FiniteElementDiscretization& fed,
      MultiMaterialNonLinearIntegrator* const mi,
      const std::vector<AbstractBoundaryCondition*>& bcs,
      const mfem::FiniteElementSpace& fes,
      const mfem::Array<int>& essential_dofs,
      mfem::Vector& f,
      mfem::Vector& r,
      const std::function<void(mfem::Vector&)>& addOtherContributions = {}) {
    auto ctx = Context{};
    auto success = false;
    try {
      f.UseDevice(true);
      f.SetSize(fes.GetVSize());
      f = real{0};
      success = (mi == nullptr) || mi->addInnerForcesByKernels(ctx, fes, f);
      for (auto* const bc : bcs) {
        success = success && bc->addContributionToTheResidual(ctx, f);
      }
      if (success && addOtherContributions) {
        static_cast<void>(f.HostReadWrite());
        addOtherContributions(f);
      }
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    raiseOnFailureOnAnyProcess(fed, ctx, success,
                               "computation of the residual");
    // from the local vector of the finite element space to the true degrees
    // of freedom
    if (const auto* const P = fes.GetProlongationMatrix(); P != nullptr) {
      P->AddMultTranspose(f, r);
    } else {
      r += f;
    }
    r.SetSubVector(essential_dofs, real{0});
    // with the debug device of MFEM, any access on the device, even a
    // reading, protects the copy on the host (`MemoryManager::GetDevicePtr`):
    // the essential degrees of freedom, read on the device by
    // `SetSubVector`, are made readable on the host, where `GetGradient`
    // reads them without declaring it
    static_cast<void>(essential_dofs.HostRead());
  }  // end of addInnerForcesComputedByKernels

  /*!
   * \brief computes the stiffness matrices of the elements by kernels, given
   * to the assembly of the tangent operator by MFEM
   * \param[in] fed: finite element discretization
   * \param[in, out] mi: multi-material integrator
   * \param[in] fes: finite element space
   */
  static void computeStiffnessMatricesByKernels(
      const FiniteElementDiscretization& fed,
      MultiMaterialNonLinearIntegrator& mi,
      const mfem::FiniteElementSpace& fes) {
    auto ctx = Context{};
    const auto success = mi.computeStiffnessMatricesByKernels(ctx, fes);
    raiseOnFailureOnAnyProcess(fed, ctx, success,
                               "computation of the tangent operator");
  }  // end of computeStiffnessMatricesByKernels

  /*!
   * \return the local matrix in which the stiffness matrices of the
   * elements are assembled by kernels (see
   * `MultiMaterialNonLinearIntegrator::makeAssembledStiffnessMatrix`)
   * \param[in] fed: finite element discretization
   * \param[in, out] mi: multi-material integrator
   * \param[in] fes: finite element space
   */
  static std::unique_ptr<mfem::SparseMatrix> makeAssembledStiffnessMatrix(
      const FiniteElementDiscretization& fed,
      MultiMaterialNonLinearIntegrator& mi,
      const mfem::FiniteElementSpace& fes) {
    auto ctx = Context{};
    auto m = mi.makeAssembledStiffnessMatrix(ctx, fes);
    raiseOnFailureOnAnyProcess(fed, ctx, m != nullptr,
                               "computation of the tangent operator");
    return m;
  }  // end of makeAssembledStiffnessMatrix

  /*!
   * \brief assembles the stiffness matrices of the elements computed by
   * kernels in the local matrix of the tangent operator, created at the
   * first call (see
   * `MultiMaterialNonLinearIntegrator::makeAssembledStiffnessMatrix`)
   * \param[in] fed: finite element discretization
   * \param[in, out] mi: multi-material integrator
   * \param[in] fes: finite element space
   * \param[in, out] m: local matrix, owned by the caller
   */
  static void assembleStiffnessMatrixByKernels(
      const FiniteElementDiscretization& fed,
      MultiMaterialNonLinearIntegrator& mi,
      const mfem::FiniteElementSpace& fes,
      mfem::SparseMatrix*& m) {
    if (m == nullptr) {
      m = makeAssembledStiffnessMatrix(fed, mi, fes).release();
    }
    auto ctx = Context{};
    const auto success = mi.assembleStiffnessMatrixByKernels(ctx, fes, *m);
    raiseOnFailureOnAnyProcess(fed, ctx, success,
                               "computation of the tangent operator");
  }  // end of assembleStiffnessMatrixByKernels

  /*!
   * \brief integrate the behaviours by batches: gradients computed element
   * by element or by kernels (no computation of the gradients for a
   * prediction), then one call per behaviour integrator
   * \return true on success
   * \param[in] mi: multi-material integrator
   * \param[in] fespace: finite element space
   * \param[in] u: unknowns (local degrees of freedom)
   * \param[in] it: integration type
   * \param[in] kernels: if true, the gradients are computed by kernels
   */
  template <typename FiniteElementSpaceType>
  static bool integrateByBatches(MultiMaterialNonLinearIntegrator& mi,
                                 const FiniteElementSpaceType& fespace,
                                 const mfem::Vector& u,
                                 const IntegrationType it,
                                 const bool kernels) {
    const auto prediction =
        (it == IntegrationType::PREDICTION_TANGENT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_SECANT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_ELASTIC_OPERATOR);
    if ((!prediction) && (kernels)) {
      auto ctx = Context{};
      if (!mi.computeGradientsByKernels(ctx, fespace, u)) {
        getErrorStream() << ctx.getErrorMessage() << '\n';
        return false;
      }
      // the behaviours integrated on the host read the gradients there
      if (!mi.integratesBatchesOnDevice()) {
        try {
          mi.synchronizeOnHost();
        } catch (std::exception& e) {
          getErrorStream() << "integrate: " << e.what() << '\n';
          return false;
        }
      }
    } else if (!prediction) {
      mfem::Array<int> vdofs;
      mfem::Vector ue;
      for (size_type i = 0; i != fespace.GetNE(); ++i) {
        const auto& e = *(fespace.GetFE(i));
        auto& tr = *(fespace.GetElementTransformation(i));
        fespace.GetElementVDofs(i, vdofs);
        u.GetSubVector(vdofs, ue);
        if (!mi.computeGradients(e, tr, ue)) {
          return false;
        }
      }
    }
    return mi.integrateBatch(it);
  }  // end of integrateByBatches

  // template <bool parallel>
  // static void export_prediction(
  //     mfem_mgis::NonLinearEvolutionProblemImplementation<parallel>& p,
  //     mfem_mgis::GridFunction<parallel>& mdu,
  //     const real te) {
  //   static auto cycle = size_type{1};
  //   auto& fed = p.getFiniteElementDiscretization();
  //   auto exporter = mfem::ParaViewDataCollection{
  //       parallel ? "result-prediction-test-parallel"
  //                : "result-prediction-test"};
  //   exporter.SetMesh(&(fed.template getMesh<parallel>()));
  //   exporter.SetDataFormat(mfem::VTKFormat::BINARY);
  //   exporter.RegisterField("OppositeOfDisplacementIncrementPrediction",
  //   &mdu); exporter.SetCycle(cycle); exporter.SetTime(te); exporter.Save();
  //   ++cycle;
  // }

  template <bool parallel>
  struct PredictionResult {
    //! \brief prediction of the opposite of the increment of the unknowns
    std::unique_ptr<mfem_mgis::GridFunction<parallel>> mdu;
    //! \brief norm of the initial residual
    const real initial_residual_norm;
  };

  [[nodiscard]] static IntegrationType convertToIntegrationType(
      const PredictionOperator o) noexcept {
    if (o == PredictionOperator::ELASTIC) {
      return IntegrationType::PREDICTION_ELASTIC_OPERATOR;
    } else if (o == PredictionOperator::SECANT) {
      return IntegrationType::PREDICTION_SECANT_OPERATOR;
    } else if (o == PredictionOperator::TANGENT) {
      return IntegrationType::PREDICTION_TANGENT_OPERATOR;
    }
    abort(
        "convertToIntegrationType: internal error, "
        "unsupported prediction operator");
  }  // end of convertToIntegrationType

  [[nodiscard]] static IntegrationType convertToIntegrationType(
      const IntegrationOperator o) noexcept {
    if (o == IntegrationOperator::ELASTIC) {
      return IntegrationType::INTEGRATION_ELASTIC_OPERATOR;
    } else if (o == IntegrationOperator::SECANT) {
      return IntegrationType::INTEGRATION_SECANT_OPERATOR;
    } else if (o == IntegrationOperator::TANGENT) {
      return IntegrationType::INTEGRATION_TANGENT_OPERATOR;
    } else if (o == IntegrationOperator::CONSISTENT_TANGENT) {
      return IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR;
    }
    abort(
        "convertToIntegrationType: internal error, "
        "unsupported integration operator");
  }  // end of convertToIntegrationType

  /*!
   * \brief data of the computation of the prediction by kernels: the
   * operator of the prediction and the inner forces are computed by kernels
   * instead of the element loops of MFEM
   */
  struct PredictionKernels {
    //! \brief multi-material integrator
    MultiMaterialNonLinearIntegrator* integrator = nullptr;
    /*!
     * \brief local matrix of the operator of the prediction, created at the
     * first call (see
     * `MultiMaterialNonLinearIntegrator::makeAssembledStiffnessMatrix`)
     */
    std::unique_ptr<mfem::SparseMatrix>* matrix = nullptr;
    //! \brief inner forces (local vector, workspace)
    mfem::Vector* inner_forces = nullptr;
  };

  /*!
   * \return true if a bilinear form has integrators
   * \param[in] a: bilinear form
   */
  template <typename BilinearFormType>
  static bool hasIntegrators(BilinearFormType& a) noexcept {
    return (a.GetDBFI()->Size() != 0) || (a.GetBBFI()->Size() != 0) ||
           (a.GetFBFI()->Size() != 0) || (a.GetBFBFI()->Size() != 0);
  }  // end of hasIntegrators

  template <bool parallel>
  [[nodiscard]] std::optional<PredictionResult<parallel>> computePrediction(
      mfem_mgis::Context& ctx,
      mfem_mgis::NonLinearEvolutionProblemImplementation<parallel>& p,
      LinearSolver& ls,
      const real t,
      const real dt,
      const bool discardLinearSolverFailure,
      const PredictionKernels& kernels) noexcept {
    auto& fed = p.getFiniteElementDiscretization();
    auto& fespace = fed.template getFiniteElementSpace<parallel>();
    const auto& u0 = p.getUnknowns(mfem_mgis::bts);
    // the phases of the prediction are timed separately (sections `PRED::`)
    auto success = [&ctx, &p, &u0, dt] {
      CatchTimeSection(ctx, "PRED::integrate");
      const auto policy = p.getPredictionPolicy();
      if (policy.strategy ==
          PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION) {
        if (policy.prediction_operator ==
            PredictionOperator::LAST_ITERATE_OPERATOR) {
          return p.integrate(
              u0, convertToIntegrationType(PredictionOperator::ELASTIC), {});
        } else {
          return p.integrate(
              u0, convertToIntegrationType(policy.prediction_operator), {});
        }
      }
      const auto ndt = [&policy, dt] {
        if (policy.null_time_increment) {
          return real{};
        }
        return dt;
      }();
      return p.integrate(
          u0, convertToIntegrationType(policy.integration_operator), ndt);
    }();
#ifdef MFEM_USE_MPI
    if constexpr (parallel) {
      success = isTrueOnAllProcesses(fed, success);
    }
#endif /* MFEM_USE_MPI */
    if (!success) {
      return ctx.registerErrorMessage("integration failure");
    }
    //
    auto profiler_forms = std::optional<mgis::ProfilingSection>{};
    profiler_forms.emplace(ctx, "PRED::forms", ctx.isProfilingEnabled());
    auto u0_local = mfem_mgis::GridFunction<parallel>(&fespace);
    if constexpr (parallel) {
      u0_local.Distribute(u0);
    } else {
      u0_local = u0;
    }
    // contributions of the boundary conditions
    auto mdu = std::make_unique<mfem_mgis::GridFunction<parallel>>(&fespace);
    *mdu = 0.0;
    auto a = mfem_mgis::BilinearForm<parallel>(&fespace);
    auto b = mfem_mgis::LinearForm<parallel>(&fespace);
    for (const auto& bc : p.getBoundaryConditions()) {
      if (!bc->addLinearFormIntegrators(ctx, a, b, u0_local, t, dt)) {
        return {};
      }
    }
    // the operator and the inner forces are computed by kernels if the
    // boundary conditions do not contribute to the operator, and, in
    // sequential, if there is no prolongation (non conforming mesh)
    const auto use_kernels =
        (kernels.integrator != nullptr) && (!hasIntegrators(a)) &&
        (parallel || (fespace.GetProlongationMatrix() == nullptr));
    if (!use_kernels) {
      auto operators = p.getLinearizedOperators(ctx, u0_local);
      if (isInvalid(operators)) {
        return {};
      }
      a.AddDomainIntegrator(operators->K.release());
      if (operators->mFi.get() != nullptr) {
        b.AddDomainIntegrator(operators->mFi.release());
      }
      a.Assemble();
    }
    b.Assemble();
    //
    auto A = [] {
#ifdef MFEM_USE_MPI
      if constexpr (parallel) {
        return mfem::HypreParMatrix{};
      } else {
        return mfem::SparseMatrix{};
      }
#else
      return mfem::SparseMatrix{};
#endif /* MFEM_USE_MPI */
    }();
    // operator built by kernels
    auto pK = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
    const mfem::Operator* op = &A;
    mfem::Vector B, X;
    //
    auto essential_dofs = p.getEssentialDegreesOfFreedom();
    // copy allocated by MFEM, readable on the device (an array wrapping the
    // memory of the vector can't be read on the debug device of MFEM)
    auto edofs_list = mfem::Array<mfem_mgis::size_type>(
        static_cast<int>(essential_dofs.size()));
    std::copy(essential_dofs.begin(), essential_dofs.end(),
              edofs_list.HostWrite());
    auto mdu_values = mfem::Vector(fespace.GetTrueVSize());
    mdu_values = 0.0;
    if constexpr (parallel) {
      //
      for (const auto& bc : p.getDirichletBoundaryConditions()) {
        // we impose the opposite of the displacement increment
        bc->setImposedValuesIncrements(mdu_values, t, t + dt, -1);
      }
      mdu->Distribute(mdu_values);
    } else {
      for (const auto& bc : p.getDirichletBoundaryConditions()) {
        // we impose the opposite of the displacement increment
        bc->setImposedValuesIncrements(mdu_values, t, t + dt, -1);
      }
      mdu->SetFromTrueDofs(mdu_values);
    }
    //
    profiler_forms.reset();
    if (use_kernels) {
      CatchTimeSection(ctx, "PRED::assemble");
      // as `FormLinearSystem` of MFEM, the elimination of the essential
      // degrees of freedom of the operator being made in place: B = Pᵀ b - K
      // X, X being null except at the essential degrees of freedom, then
      // B[e] = X[e] (diagonal of the operator set to one)
      try {
        auto& mi = *(kernels.integrator);
        if (kernels.matrix->get() == nullptr) {
          *(kernels.matrix) = makeAssembledStiffnessMatrix(fed, mi, fespace);
        }
        auto* m = kernels.matrix->get();
        assembleStiffnessMatrixByKernels(fed, mi, fespace, m);
        X = mdu_values;
        B.SetSize(fespace.GetTrueVSize());
        B.UseDevice(true);
        if (const auto* const P = fespace.GetProlongationMatrix();
            P != nullptr) {
          P->MultTranspose(b, B);
        } else {
          B = b;
        }
        // the boundary conditions contribute to `b` by their linear forms
        addInnerForcesComputedByKernels(fed, &mi, {}, fespace, edofs_list,
                                        *(kernels.inner_forces), B);
        // elimination in the local matrix, without prolongation
        const auto eliminateInLocalMatrix = [&X, &B, &edofs_list, &op, m] {
          m->AddMult(X, B, real{-1});
          m->EliminateBC(edofs_list, mfem::Operator::DIAG_ONE);
          // the linear solvers may read the matrix on the host without
          // declaring it
          static_cast<void>(m->HostReadI());
          static_cast<void>(m->HostReadJ());
          static_cast<void>(m->HostReadData());
          op = m;
        };
#ifdef MFEM_USE_MPI
        if constexpr (parallel) {
          mfem::OperatorHandle dA(mfem::Operator::Hypre_ParCSR);
          mfem::OperatorHandle Ph(mfem::Operator::Hypre_ParCSR);
          dA.MakeSquareBlockDiag(fespace.GetComm(), fespace.GlobalVSize(),
                                 fespace.GetDofOffsets(), m);
          Ph.ConvertFrom(fespace.Dof_TrueDof_Matrix());
          pK.MakePtAP(dA, Ph);
          auto& K = *(pK.As<mfem::HypreParMatrix>());
          K.Mult(real{-1}, X, real{1}, B);
          K.EliminateBC(edofs_list, mfem::Operator::DIAG_ONE);
          op = &K;
        } else {
          eliminateInLocalMatrix();
        }
#else  /* MFEM_USE_MPI */
        eliminateInLocalMatrix();
#endif /* MFEM_USE_MPI */
        auto values = mfem::Vector{};
        X.GetSubVector(edofs_list, values);
        B.SetSubVector(edofs_list, values);
      } catch (...) {
        return registerExceptionInErrorBacktrace(ctx);
      }
    } else {
      // diagonal set to one at the essential degrees of freedom, as in
      // parallel and with the kernels (the sequential bilinear forms of MFEM
      // keep the diagonal of the stiffness by default): an iterative linear
      // solver preconditioned by the elastic operator, whose diagonal is set
      // to one, converges in one iteration
      a.SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
      a.FormLinearSystem(edofs_list, *mdu, b, A, X, B);
    }
    // the essential degrees of freedom are excluded from the norm: the
    // residual of the non linear solver is null for them
    const auto norm = [&B, &fespace, &essential_dofs] {
      auto r = B;
      auto* const rv = r.HostReadWrite();
      for (const auto i : essential_dofs) {
        rv[i] = real{0};
      }
      if constexpr (parallel) {
        return std::sqrt(mfem::InnerProduct(fespace.GetComm(), r, r));
      } else {
        return std::sqrt(mfem::InnerProduct(r, r));
      }
    }();
    {
      CatchTimeSection(ctx, "PRED::solve");
      ls.SetOperator(*op);
      ls.Mult(B, X);
      clearHypreErrors();
    }
    // check for convergence
    if ((!discardLinearSolverFailure) && (!hasConverged(ls))) {
      return ctx.registerErrorMessage("linear solver did not converge");
    }
    // an iterative linear solver only approximates the imposed increments,
    // which are known exactly
    {
      auto* const x = X.HostReadWrite();
      const auto* const mv = mdu_values.HostRead();
      for (const auto i : essential_dofs) {
        x[i] = mv[i];
      }
    }
    if (use_kernels) {
      if constexpr (parallel) {
        mdu->Distribute(X);
      } else {
        mdu->SetFromTrueDofs(X);
      }
    } else {
      a.RecoverFEMSolution(X, b, *mdu);
    }
    //
    return PredictionResult<parallel>{.mdu = std::move(mdu),
                                      .initial_residual_norm = norm};
  }  // end of computePrediction

  /*!
   * \brief post-processing defined by an std::function
   * \tparam parallel: boolean stating if the computations are performed in
   * parallel.
   */
  template <bool parallel>
  struct StdFunctionPostProcessing final
      : NonLinearEvolutionProblemPostProcessingBase<parallel> {
    /*!
     * \brief constructor
     * \param[in] fct: function executing the postprocessing
     */
    StdFunctionPostProcessing(
        const std::function<void(const real, const real)>& fct)
        : f(fct) {
      if (!f) {
        raise("invalid function");
      }
    }  // end of StdFunctionPostProcessing
    //
    [[nodiscard]] bool execute(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const real t,
        const real dt) noexcept override {
      auto success = true;
      auto fed = p.getFiniteElementDiscretization();
      try {
        this->f(t, dt);
      } catch (...) {
        std::ignore = registerExceptionInErrorBacktrace(ctx);
        success = false;
      }
      return isTrueOnAllProcesses(fed, success);
    }  // end of execute
    //! \brief destructor
    ~StdFunctionPostProcessing() override = default;

   private:
    //! \brief function executing the post-processing
    std::function<void(const real, const real)> f;
  };  // end of struct StdFunctionPostProcessing

#ifdef MFEM_USE_MPI

  NonLinearEvolutionProblemImplementation<true>::
      NonLinearEvolutionProblemImplementation(
          Context& ctx,
          std::shared_ptr<FiniteElementDiscretization> fed,
          const Hypothesis h,
          const Parameters& p)
      : NonLinearEvolutionProblemImplementationBase(
            ctx,
            fed,
            h,
            extract(throwing,
                    p,
                    NonLinearEvolutionProblemImplementationBase::
                        getParametersList())),
        mfem::ParNonlinearForm(&(fed->getFiniteElementSpace<true>())) {
    //
    auto or_raise = ctx.getThrowingFailureHandler();
    auto validator = ParametersValidator{}.add(
        NonLinearEvolutionProblemImplementationBase::getParametersList());
    validator.validate(ctx, p) | or_raise;
    //
    if (this->fe_discretization->getMesh<true>().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      raise(
          "NonLinearEvolutionProblemImplementationBase::"
          "NonLinearEvolutionProblemImplementationBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    auto set_solver = [&ctx, this, &p, &or_raise] {
      if (contains(
              p,
              NonLinearEvolutionProblemImplementationBase::NonLinearSolver)) {
        const auto [n, params] = extractFactoryArgument(
            throwing,
            get<Parameters>(
                throwing, p,
                NonLinearEvolutionProblemImplementationBase::NonLinearSolver));
        auto& f = NonLinearSolverFactory::get();
        this->solver = f.generate(ctx, n, *this, params) | or_raise;
      } else {
        this->solver = std::make_unique<NewtonSolver>(*this);
      }
    };
#ifdef MFEM_USE_PETSC
    if (usePETSc()) {
      this->petsc_solver = std::make_unique<mfem::PetscNonlinearSolver>(
          this->getFiniteElementSpace().GetComm(), *this);
      this->petsc_solver->iterative_mode = true;
    } else {
      set_solver();
    }
#else  /* MFEM_USE_PETSC */
    set_solver();
#endif /* MFEM_USE_PETSC */
    if (this->mgis_integrator != nullptr) {
      this->AddDomainIntegrator(this->mgis_integrator);
    }
  }  // end of NonLinearEvolutionProblemImplementation

  std::optional<real>
  NonLinearEvolutionProblemImplementation<true>::computePrediction(
      Context& ctx, const real t, const real dt) noexcept {
    if (this->linear_solver.get() == nullptr) {
      return ctx.registerErrorMessage("linear solver not initialized");
    }
    const auto discardLinearSolverFailure = [this] {
      if (this->solver.get() != nullptr) {
        return this->solver->isLinearSolverFailureDiscarded();
      }
      return false;
    }();
    auto kernels = PredictionKernels{};
    if ((this->batch_integration) && (this->inner_forces_kernels) &&
        this->computesStiffnessMatricesByKernels() &&
        (this->mgis_integrator != nullptr)) {
      kernels =
          PredictionKernels{.integrator = this->mgis_integrator,
                            .matrix = &(this->prediction_matrix),
                            .inner_forces = &(this->inner_forces_by_kernels)};
    }
    // the non linear solver may provide the linear solver of an elastic
    // prediction at the beginning of the time step
    const auto& policy = this->getPredictionPolicy();
    const auto elastic_prediction =
        (policy.strategy ==
         PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION) &&
        (policy.prediction_operator == PredictionOperator::ELASTIC);
    auto* const ls = ((this->solver.get() != nullptr) && elastic_prediction)
                         ? this->solver->getElasticPredictionLinearSolver()
                         : nullptr;
    const auto oresult = ::mfem_mgis::computePrediction<true>(
        ctx, *this, (ls != nullptr) ? *ls : *(this->linear_solver), t, dt,
        discardLinearSolverFailure, kernels);
#pragma message("required ?")
    this->solver->setLinearSolver(*(this->linear_solver));
    if (isInvalid(oresult)) {
      return {};
    }
    this->u1 = this->u0;
    this->u1 -= oresult->mdu->GetTrueVector();
    return oresult->initial_residual_norm;
  }  // end of computePrediction

  void NonLinearEvolutionProblemImplementation<true>::Mult(
      const mfem::Vector& u, mfem::Vector& r) const {
    const auto kernels =
        (this->inner_forces_kernels) && (this->mgis_integrator != nullptr);
    if ((!kernels) && (this->residual_boundary_conditions.empty())) {
      mfem_mgis::NonlinearForm<true>::Mult(u, r);
      return;
    }
    // with the inner forces computed by kernels, the integrator of mfem-mgis
    // gives null contributions to the residual assembled by MFEM: the
    // assembly of MFEM is skipped if the form has no other integrator
    const auto only_mfem_mgis =
        kernels && isTheOnlyIntegrator(*this, this->mgis_integrator);
    // the contributions of the other integrators of the form are added to
    // the local vector of the inner forces computed by kernels, if possible
    const auto others = kernels && (!only_mfem_mgis) &&
                        supportsOtherIntegrators(*this, this->mgis_integrator);
    if (only_mfem_mgis || others) {
      r = real{0};
    } else {
      mfem_mgis::NonlinearForm<true>::Mult(u, r);
    }
    auto addOtherContributions = std::function<void(mfem::Vector&)>{};
    if (others) {
      ++(this->other_integrators_residual_evaluations);
      // the unknowns are prolongated on all the processes (exchanges
      // between processes) before the computations of the residual, which
      // may fail on one process only
      const auto& pu = this->Prolongate(u);
      addOtherContributions = [this, &pu](mfem::Vector& f) {
        addOtherIntegratorsContributions<true>(
            *this, this->mgis_integrator, pu,
            [&f](const mfem::Array<int>& dofs, const mfem::Vector& v) {
              f.AddElementVector(dofs, v);
            });
      };
    }
    addInnerForcesComputedByKernels(
        this->getFiniteElementDiscretization(),
        kernels ? this->mgis_integrator : nullptr,
        this->residual_boundary_conditions, this->getFiniteElementSpace(),
        this->ess_tdof_list, this->inner_forces_by_kernels, r,
        addOtherContributions);
  }  // end of Mult

  mfem::Operator& NonLinearEvolutionProblemImplementation<true>::
      assembleTangentOperatorByKernels(const mfem::Vector* const u) const {
    if ((u == nullptr) &&
        (!isTheOnlyIntegrator(*this, this->mgis_integrator))) {
      raise(
          "assembleTangentOperatorByKernels: the unknowns are required to "
          "add the contributions of the other integrators of the non linear "
          "form");
    }
    // the stiffness matrices are assembled by kernels in the local matrix,
    // then the global matrix is built as by
    // `mfem::ParNonlinearForm::GetGradient`
    ++(*(this->assembled_tangent_operator_revision));
    assembleStiffnessMatrixByKernels(this->getFiniteElementDiscretization(),
                                     *(this->mgis_integrator),
                                     this->getFiniteElementSpace(), this->Grad);
    if (u != nullptr) {
      // the unknowns are prolongated on all the processes (exchanges
      // between processes) before the evaluation of the other integrators
      addOtherIntegratorsContributionsToTheLocalMatrix(
          this->getFiniteElementDiscretization(), *this, this->mgis_integrator,
          this->Prolongate(*u), *(this->Grad));
    }
    auto* const pfes = this->ParFESpace();
    if (this->pGrad.Type() == mfem::Operator::Hypre_ParCSR) {
      // global matrix whose pattern and exchanges between processes are
      // computed once, if the finite element space allows it (checked once,
      // collectively)
      if (!this->global_tangent_operator_checked) {
        this->global_tangent_operator_checked = true;
        if (batch::GlobalStiffnessMatrix::isSupported(*pfes)) {
          this->global_tangent_operator =
              std::make_unique<batch::GlobalStiffnessMatrix>(*pfes);
        }
      }
      if (this->global_tangent_operator != nullptr) {
        auto ctx = Context{};
        auto* const K =
            this->global_tangent_operator->assemble(ctx, *(this->Grad));
        raiseOnFailureOnAnyProcess(this->getFiniteElementDiscretization(), ctx,
                                   K != nullptr,
                                   "computation of the tangent operator");
        // the matrix is owned by `global_tangent_operator`
        this->pGrad.Reset(K, false);
        K->EliminateBC(this->ess_tdof_list, mfem::Operator::DIAG_ONE);
        return *K;
      }
    }
    this->pGrad.Clear();
    mfem::OperatorHandle dA(this->pGrad.Type());
    mfem::OperatorHandle Ph(this->pGrad.Type());
    dA.MakeSquareBlockDiag(pfes->GetComm(), pfes->GlobalVSize(),
                           pfes->GetDofOffsets(), this->Grad);
    Ph.ConvertFrom(pfes->Dof_TrueDof_Matrix());
    this->pGrad.MakePtAP(dA, Ph);
    if (this->pGrad.Type() == mfem::Operator::Hypre_ParCSR) {
      // same result as `EliminateRowsCols`, the pattern of the matrix
      // being symmetric, but in place (on the device if hypre uses it),
      // without computing the eliminated part
      this->pGrad.As<mfem::HypreParMatrix>()->EliminateBC(
          this->ess_tdof_list, mfem::Operator::DIAG_ONE);
    } else {
      mfem::OperatorHandle pGrad_e;
      pGrad_e.EliminateRowsCols(this->pGrad, this->ess_tdof_list);
    }
    return *(this->pGrad.Ptr());
  }  // end of assembleTangentOperatorByKernels

  mfem::Operator& NonLinearEvolutionProblemImplementation<true>::GetGradient(
      const mfem::Vector& u) const {
    if (this->appliesTangentOperatorWithoutAssembly()) {
      checkMatrixFreeTangentOperatorSupport(*this, this->mgis_integrator);
      return this->getMatrixFreeTangentOperator(
          this->ess_tdof_list, [this]() -> mfem::Operator& {
            return this->assembleTangentOperatorByKernels();
          });
    }
    if ((!this->computesStiffnessMatricesByKernels()) ||
        (this->mgis_integrator == nullptr)) {
      return mfem_mgis::NonlinearForm<true>::GetGradient(u);
    }
    if (isTheOnlyIntegrator(*this, this->mgis_integrator)) {
      return this->assembleTangentOperatorByKernels();
    }
    if (supportsOtherIntegrators(*this, this->mgis_integrator)) {
      return this->assembleTangentOperatorByKernels(&u);
    }
    // other integrators in the form not supported by
    // `addOtherIntegratorsContributions`: the stiffness matrices computed by
    // kernels are copied by the element loop of MFEM
    computeStiffnessMatricesByKernels(this->getFiniteElementDiscretization(),
                                      *(this->mgis_integrator),
                                      this->getFiniteElementSpace());
    return mfem_mgis::NonlinearForm<true>::GetGradient(u);
  }  // end of GetGradient

  bool NonLinearEvolutionProblemImplementation<true>::addBoundaryCondition(
      Context& ctx,
      std::unique_ptr<AbstractDirichletBoundaryCondition> bc) noexcept {
    if (bc.get() == nullptr) {
      return ctx.registerErrorMessage("invalid boundary condition");
    }
    this->dirichlet_boundary_conditions.push_back(std::move(bc));
    return true;
  }  // end of addBoundaryCondition

  bool NonLinearEvolutionProblemImplementation<true>::addBoundaryCondition(
      Context& ctx, std::unique_ptr<AbstractBoundaryCondition> f) noexcept {
    if (f.get() == nullptr) {
      return ctx.registerErrorMessage("invalid boundary condition");
    }
    if ((this->inner_forces_kernels) && (f->isIndependentOfTheUnknowns())) {
      // the contribution of the boundary condition is added to the inner
      // forces computed by kernels (see `Mult`): the non linear form only
      // contains the integrator of mfem-mgis, and the tangent operator can
      // be assembled by kernels or applied without assembly
      this->residual_boundary_conditions.push_back(f.get());
    } else if (!f->addNonlinearFormIntegrator(ctx, *this, this->u1)) {
      return false;
    }
    this->boundary_conditions.push_back(std::move(f));
    return true;
  }  // end of addBoundaryCondition

  bool NonLinearEvolutionProblemImplementation<true>::addPostProcessing(
      Context& ctx,
      std::unique_ptr<AbstractNonLinearEvolutionProblemPostProcessing<true>>
          p) noexcept {
    if (p.get() == nullptr) {
      return ctx.registerErrorMessage("invalid post-processing");
    }
    this->postprocessings.push_back(std::move(p));
    return true;
  }  // end of addPostProcessing

  bool NonLinearEvolutionProblemImplementation<true>::addPostProcessing(
      Context& ctx, std::string_view n, const Parameters& p) noexcept {
    const auto& f = PostProcessingFactory<true>::getFactory();
    return this->addPostProcessing(ctx, f.generate(ctx, n, *this, p));
  }  // end of addPostProcessing

  bool NonLinearEvolutionProblemImplementation<true>::setLinearSolver(
      Context& ctx, LinearSolverHandler s) noexcept {
    if (isInvalid(s)) {
      return ctx.registerErrorMessage("invalid linear solver");
    }
    this->updateLinearSolver(ctx, std::move(s));
    return true;
  }

  bool NonLinearEvolutionProblemImplementation<true>::setLinearSolver(
      Context& ctx, std::string_view n, const Parameters& p) noexcept {
    const auto& f = LinearSolverFactory<true>::getFactory();
    auto& fespace = this->getFiniteElementSpace();
    auto s = f.generate(ctx, n, fespace, p);
    if (isInvalid(s)) {
      return false;
    }
    this->updateLinearSolver(ctx, std::move(s));
    return true;
  }  // end of setLinearSolver

  bool NonLinearEvolutionProblemImplementation<true>::addPostProcessing(
      Context& ctx,
      const std::function<void(const real, const real)>& p) noexcept {
    if (!p) {
      return ctx.registerErrorMessage("invalid post-processing function");
    }
    auto ptr = make_unique<StdFunctionPostProcessing<true>>(ctx, p);
    if (isInvalid(ptr)) {
      return false;
    }
    return this->addPostProcessing(ctx, std::move(ptr));
  }  // end of addPostProcessing

  bool
  NonLinearEvolutionProblemImplementation<true>::executeInitialPostProcessings(
      Context& ctx, const real t) noexcept {
    // the post-processings read the materials on the host
    if (this->mgis_integrator != nullptr) {
      try {
        this->mgis_integrator->synchronizeOnHostForReading();
      } catch (...) {
        return registerExceptionInErrorBacktrace(ctx);
      }
    }
    for (auto& p : this->postprocessings) {
      if (!p->hasExecuteInitialPostProcessingAlreadyBeenCalled()) {
        if (!p->executeInitialPostProcessing(ctx, *this, t)) {
          return false;
        }
      }
    }
    return true;
  }  // end of executeInitialPostProcessings

  bool NonLinearEvolutionProblemImplementation<true>::executePostProcessings(
      Context& ctx, const real t, const real dt) noexcept {
    CatchTimeSection(ctx, "NLEP::executePostProcessings");
    // the post-processings read the materials on the host
    if (this->mgis_integrator != nullptr) {
      CatchTimeSection(ctx, "NLEP::PP::synchronizeOnHost");
      try {
        this->mgis_integrator->synchronizeOnHostForReading();
      } catch (...) {
        return registerExceptionInErrorBacktrace(ctx);
      }
    }
    auto success = true;
    for (auto& p : this->postprocessings) {
      CatchTimeSection(ctx, "NLEP::PP::execute");
      if (!p->execute(ctx, *this, t, dt)) {
        success = false;
      }
    }
    return success;
  }  // end of executePostProcessings

  Mesh<true>& NonLinearEvolutionProblemImplementation<true>::getMesh() {
    return this->fe_discretization->getMesh<true>();
  }  // end of getMesh

  const Mesh<true>& NonLinearEvolutionProblemImplementation<true>::getMesh()
      const {
    return this->fe_discretization->getMesh<true>();
  }  // end of getMesh

  FiniteElementSpace<true>&
  NonLinearEvolutionProblemImplementation<true>::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace<true>();
  }  // end of getFiniteElementSpace

  const FiniteElementSpace<true>&
  NonLinearEvolutionProblemImplementation<true>::getFiniteElementSpace() const {
    return this->fe_discretization->getFiniteElementSpace<true>();
  }  // end of getFiniteElementSpace

  bool NonLinearEvolutionProblemImplementation<true>::integrate(
      const mfem::Vector& u,
      const IntegrationType it,
      const std::optional<real> odt) {
    if (this->mgis_integrator == nullptr) {
      return true;
    }
    // the behaviours are integrated on the host, unless they are all
    // integrated on the device from gradients computed by kernels (or
    // without gradients, for a prediction)
    if (!integratesOnDevice(*(this->mgis_integrator), this->batch_integration,
                            this->gradients_kernels, it)) {
      try {
        this->mgis_integrator->synchronizeOnHost();
      } catch (std::exception& e) {
        getErrorStream() << "integrate: " << e.what() << '\n';
        return false;
      } catch (...) {
        getErrorStream() << "integrate: unknown exception\n";
        return false;
      }
    }
    const auto dt = this->mgis_integrator->getTimeIncrement();
    if (odt.has_value()) {
      this->mgis_integrator->setTimeIncrement(*odt);
    }
    const auto& pu = this->Prolongate(u);
    const auto& fespace = this->getFiniteElementSpace();
    mfem::Array<int> vdofs;
    mfem::Vector ue;
    bool success = true;
    if (this->batch_integration) {
      success = integrateByBatches(*(this->mgis_integrator), fespace, pu, it,
                                   this->gradients_kernels);
    } else {
      for (size_type i = 0; success && (i != fespace.GetNE()); ++i) {
        const auto& e = *(fespace.GetFE(i));
        auto& tr = *(fespace.GetElementTransformation(i));
        fespace.GetElementVDofs(i, vdofs);
        pu.GetSubVector(vdofs, ue);
        success = this->mgis_integrator->integrate(e, tr, ue, it);
      }
    }
    success = isTrueOnAllProcesses(*(this->fe_discretization), success);
    if ((success) && (it != IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR)) {
      this->hasStiffnessOperatorsBeenComputed = true;
    }
    this->mgis_integrator->setTimeIncrement(dt);
    return success;
  }  // end of integrate

  void NonLinearEvolutionProblemImplementation<true>::
      markDegreesOfFreedomHandledByDirichletBoundaryConditions(
          std::vector<size_type> dofs) {
    auto tmp = mfem::Array<size_type>(dofs.data(), dofs.size());
    this->SetEssentialTrueDofs(tmp);
  }  // end of markDegreesOfFreedomHandledByDirichletBoundaryConditions

  NonLinearEvolutionProblemImplementation<
      true>::~NonLinearEvolutionProblemImplementation() = default;

#endif /* MFEM_USE_MPI */

  NonLinearEvolutionProblemImplementation<false>::
      NonLinearEvolutionProblemImplementation(
          Context& ctx,
          std::shared_ptr<FiniteElementDiscretization> fed,
          const Hypothesis h,
          const Parameters& p)
      : NonLinearEvolutionProblemImplementationBase(
            ctx,
            fed,
            h,
            extract(throwing,
                    p,
                    NonLinearEvolutionProblemImplementationBase::
                        getParametersList())),
        mfem::NonlinearForm(&(fed->getFiniteElementSpace<false>())) {
    //
    auto or_raise = ctx.getThrowingFailureHandler();
    auto validator = ParametersValidator{}.add(
        NonLinearEvolutionProblemImplementationBase::getParametersList());
    validator.validate(ctx, p) | or_raise;
    //
    if (this->fe_discretization->getMesh<false>().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      raise(
          "NonLinearEvolutionProblemImplementationBase::"
          "NonLinearEvolutionProblemImplementationBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    if (this->mgis_integrator != nullptr) {
      this->AddDomainIntegrator(this->mgis_integrator);
    }
    //
    if (contains(
            p, NonLinearEvolutionProblemImplementationBase::NonLinearSolver)) {
      const auto [n, params] = extractFactoryArgument(
          throwing,
          get<Parameters>(
              throwing, p,
              NonLinearEvolutionProblemImplementationBase::NonLinearSolver));
      auto& f = NonLinearSolverFactory::get();
      this->solver = f.generate(ctx, n, *this, params) | or_raise;
    } else {
      this->solver = std::make_unique<NewtonSolver>(*this);
    }
  }  // end of NonLinearEvolutionProblemImplementation

  std::optional<real>
  NonLinearEvolutionProblemImplementation<false>::computePrediction(
      Context& ctx, const real t, const real dt) noexcept {
    if (this->linear_solver.get() == nullptr) {
      return ctx.registerErrorMessage("linear solver not initialized");
    }
    const auto discardLinearSolverFailure = [this] {
      if (this->solver.get() != nullptr) {
        return this->solver->isLinearSolverFailureDiscarded();
      }
      return false;
    }();
    auto kernels = PredictionKernels{};
    if ((this->batch_integration) && (this->inner_forces_kernels) &&
        this->computesStiffnessMatricesByKernels() &&
        (this->mgis_integrator != nullptr)) {
      kernels =
          PredictionKernels{.integrator = this->mgis_integrator,
                            .matrix = &(this->prediction_matrix),
                            .inner_forces = &(this->inner_forces_by_kernels)};
    }
    // the non linear solver may provide the linear solver of an elastic
    // prediction at the beginning of the time step
    const auto& policy = this->getPredictionPolicy();
    const auto elastic_prediction =
        (policy.strategy ==
         PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION) &&
        (policy.prediction_operator == PredictionOperator::ELASTIC);
    auto* const ls = ((this->solver.get() != nullptr) && elastic_prediction)
                         ? this->solver->getElasticPredictionLinearSolver()
                         : nullptr;
    const auto oresult = ::mfem_mgis::computePrediction<false>(
        ctx, *this, (ls != nullptr) ? *ls : *(this->linear_solver), t, dt,
        discardLinearSolverFailure, kernels);
    this->solver->setLinearSolver(*(this->linear_solver));
    if (isInvalid(oresult)) {
      return {};
    }
    this->u1 = this->u0;
    this->u1 -= oresult->mdu->GetTrueVector();
    return oresult->initial_residual_norm;
  }  // end of computePrediction

  void NonLinearEvolutionProblemImplementation<false>::Mult(
      const mfem::Vector& u, mfem::Vector& r) const {
    const auto kernels =
        (this->inner_forces_kernels) && (this->mgis_integrator != nullptr);
    if ((!kernels) && (this->residual_boundary_conditions.empty())) {
      mfem_mgis::NonlinearForm<false>::Mult(u, r);
      return;
    }
    const auto only_mfem_mgis =
        kernels && isTheOnlyIntegrator(*this, this->mgis_integrator);
    // the contributions of the other integrators of the form are added to
    // the local vector of the inner forces computed by kernels, if possible
    const auto others = kernels && (!only_mfem_mgis) &&
                        supportsOtherIntegrators(*this, this->mgis_integrator);
    if (only_mfem_mgis || others) {
      r = real{0};
    } else {
      mfem_mgis::NonlinearForm<false>::Mult(u, r);
    }
    auto addOtherContributions = std::function<void(mfem::Vector&)>{};
    if (others) {
      ++(this->other_integrators_residual_evaluations);
      // the unknowns are prolongated (non conforming meshes) before the
      // computations of the residual
      const auto& pu = this->Prolongate(u);
      addOtherContributions = [this, &pu](mfem::Vector& f) {
        addOtherIntegratorsContributions<true>(
            *this, this->mgis_integrator, pu,
            [&f](const mfem::Array<int>& dofs, const mfem::Vector& v) {
              f.AddElementVector(dofs, v);
            });
      };
    }
    addInnerForcesComputedByKernels(
        this->getFiniteElementDiscretization(),
        kernels ? this->mgis_integrator : nullptr,
        this->residual_boundary_conditions, this->getFiniteElementSpace(),
        this->ess_tdof_list, this->inner_forces_by_kernels, r,
        addOtherContributions);
  }  // end of Mult

  mfem::Operator& NonLinearEvolutionProblemImplementation<false>::
      assembleTangentOperatorByKernels(const mfem::Vector* const u) const {
    if ((u == nullptr) &&
        (!isTheOnlyIntegrator(*this, this->mgis_integrator))) {
      raise(
          "assembleTangentOperatorByKernels: the unknowns are required to "
          "add the contributions of the other integrators of the non linear "
          "form");
    }
    ++(*(this->assembled_tangent_operator_revision));
    assembleStiffnessMatrixByKernels(this->getFiniteElementDiscretization(),
                                     *(this->mgis_integrator),
                                     this->getFiniteElementSpace(), this->Grad);
    if (u != nullptr) {
      addOtherIntegratorsContributionsToTheLocalMatrix(
          this->getFiniteElementDiscretization(), *this, this->mgis_integrator,
          this->Prolongate(*u), *(this->Grad));
    }
    // as `mfem::NonlinearForm::GetGradient`, the elimination of the
    // essential degrees of freedom giving the same result as
    // `EliminateRowCol`, the pattern of the matrix being symmetric
    auto* mGrad = this->Grad;
    if (this->cP != nullptr) {
      // `RAP` reads the matrix on the host without declaring it
      static_cast<void>(this->Grad->HostReadI());
      static_cast<void>(this->Grad->HostReadJ());
      static_cast<void>(this->Grad->HostReadData());
      delete this->cGrad;
      this->cGrad = mfem::RAP(*(this->cP), *(this->Grad), *(this->cP));
      mGrad = this->cGrad;
    }
    mGrad->EliminateBC(this->ess_tdof_list, mfem::Operator::DIAG_ONE);
    // the linear solvers may read the matrix on the host without declaring
    // it
    static_cast<void>(mGrad->HostReadI());
    static_cast<void>(mGrad->HostReadJ());
    static_cast<void>(mGrad->HostReadData());
    return *mGrad;
  }  // end of assembleTangentOperatorByKernels

  mfem::Operator& NonLinearEvolutionProblemImplementation<false>::GetGradient(
      const mfem::Vector& u) const {
    if (this->appliesTangentOperatorWithoutAssembly()) {
      checkMatrixFreeTangentOperatorSupport(*this, this->mgis_integrator);
      return this->getMatrixFreeTangentOperator(
          this->ess_tdof_list, [this]() -> mfem::Operator& {
            return this->assembleTangentOperatorByKernels();
          });
    }
    if ((!this->computesStiffnessMatricesByKernels()) ||
        (this->mgis_integrator == nullptr)) {
      return mfem_mgis::NonlinearForm<false>::GetGradient(u);
    }
    if (isTheOnlyIntegrator(*this, this->mgis_integrator)) {
      return this->assembleTangentOperatorByKernels();
    }
    if (supportsOtherIntegrators(*this, this->mgis_integrator)) {
      return this->assembleTangentOperatorByKernels(&u);
    }
    computeStiffnessMatricesByKernels(this->getFiniteElementDiscretization(),
                                      *(this->mgis_integrator),
                                      this->getFiniteElementSpace());
    return mfem_mgis::NonlinearForm<false>::GetGradient(u);
  }  // end of GetGradient

  bool NonLinearEvolutionProblemImplementation<false>::addBoundaryCondition(
      Context& ctx,
      std::unique_ptr<AbstractDirichletBoundaryCondition> bc) noexcept {
    if (bc.get() == nullptr) {
      return ctx.registerErrorMessage("invalid boundary condition");
    }
    this->dirichlet_boundary_conditions.push_back(std::move(bc));
    return true;
  }  // end of addBoundaryCondition

  bool NonLinearEvolutionProblemImplementation<false>::addBoundaryCondition(
      Context& ctx, std::unique_ptr<AbstractBoundaryCondition> f) noexcept {
    if (f.get() == nullptr) {
      return ctx.registerErrorMessage("invalid boundary condition");
    }
    if ((this->inner_forces_kernels) && (f->isIndependentOfTheUnknowns())) {
      // the contribution of the boundary condition is added to the inner
      // forces computed by kernels (see `Mult`): the non linear form only
      // contains the integrator of mfem-mgis, and the tangent operator can
      // be assembled by kernels or applied without assembly
      this->residual_boundary_conditions.push_back(f.get());
    } else if (!f->addNonlinearFormIntegrator(ctx, *this, this->u1)) {
      return false;
    }
    this->boundary_conditions.push_back(std::move(f));
    return true;
  }  // end of addBoundaryCondition

  bool NonLinearEvolutionProblemImplementation<false>::addPostProcessing(
      Context& ctx,
      std::unique_ptr<AbstractNonLinearEvolutionProblemPostProcessing<false>>
          p) noexcept {
    if (p.get() == nullptr) {
      return ctx.registerErrorMessage("invalid post-processing");
    }
    this->postprocessings.push_back(std::move(p));
    return true;
  }  // end of addPostProcessing

  bool NonLinearEvolutionProblemImplementation<false>::addPostProcessing(
      Context& ctx,
      const std::function<void(const real, const real)>& p) noexcept {
    if (!p) {
      return ctx.registerErrorMessage("invalid post-processing function");
    }
    auto ptr = make_unique<StdFunctionPostProcessing<false>>(ctx, p);
    if (isInvalid(ptr)) {
      return false;
    }
    return this->addPostProcessing(ctx, std::move(ptr));
  }  // end of addPostProcessing

  bool NonLinearEvolutionProblemImplementation<false>::addPostProcessing(
      Context& ctx, std::string_view n, const Parameters& p) noexcept {
    const auto& f = PostProcessingFactory<false>::getFactory();
    return this->addPostProcessing(ctx, f.generate(ctx, n, *this, p));
  }  // end of addPostProcessing

  bool
  NonLinearEvolutionProblemImplementation<false>::executeInitialPostProcessings(
      Context& ctx, const real t) noexcept {
    // the post-processings read the materials on the host
    if (this->mgis_integrator != nullptr) {
      try {
        this->mgis_integrator->synchronizeOnHostForReading();
      } catch (...) {
        return registerExceptionInErrorBacktrace(ctx);
      }
    }
    for (auto& p : this->postprocessings) {
      if (!p->hasExecuteInitialPostProcessingAlreadyBeenCalled()) {
        if (!p->executeInitialPostProcessing(ctx, *this, t)) {
          return false;
        }
      }
    }
    return true;
  }  // end of executeInitialPostProcessings

  bool NonLinearEvolutionProblemImplementation<false>::executePostProcessings(
      Context& ctx, const real t, const real dt) noexcept {
    CatchTimeSection(ctx, "NLEP::executePostProcessings");
    // the post-processings read the materials on the host
    if (this->mgis_integrator != nullptr) {
      CatchTimeSection(ctx, "NLEP::PP::synchronizeOnHost");
      try {
        this->mgis_integrator->synchronizeOnHostForReading();
      } catch (...) {
        return registerExceptionInErrorBacktrace(ctx);
      }
    }
    auto success = true;
    for (auto& p : this->postprocessings) {
      CatchTimeSection(ctx, "NLEP::PP::execute");
      if (!p->execute(ctx, *this, t, dt)) {
        success = false;
      }
    }
    return success;
  }  // end of executePostProcessings

  Mesh<false>&
  NonLinearEvolutionProblemImplementation<false>::getMesh() noexcept {
    return this->fe_discretization->getMesh<false>();
  }  // end of getMesh

  const Mesh<false>& NonLinearEvolutionProblemImplementation<false>::getMesh()
      const noexcept {
    return this->fe_discretization->getMesh<false>();
  }  // end of getMesh

  const FiniteElementSpace<false>& NonLinearEvolutionProblemImplementation<
      false>::getFiniteElementSpace() const noexcept {
    return this->fe_discretization->getFiniteElementSpace<false>();
  }  // end of getFiniteElementSpace

  FiniteElementSpace<false>& NonLinearEvolutionProblemImplementation<
      false>::getFiniteElementSpace() noexcept {
    return this->fe_discretization->getFiniteElementSpace<false>();
  }  // end of getFiniteElementSpace

  bool NonLinearEvolutionProblemImplementation<false>::setLinearSolver(
      Context& ctx, LinearSolverHandler s) noexcept {
    if (isInvalid(s)) {
      return ctx.registerErrorMessage("invalid linear solver");
    }
    this->updateLinearSolver(ctx, std::move(s));
    return true;
  }

  bool NonLinearEvolutionProblemImplementation<false>::setLinearSolver(
      Context& ctx, std::string_view n, const Parameters& p) noexcept {
    const auto& f = LinearSolverFactory<false>::getFactory();
    auto& fespace = this->getFiniteElementSpace();
    auto s = f.generate(ctx, n, fespace, p);
    if (isInvalid(s)) {
      return false;
    }
    this->updateLinearSolver(ctx, std::move(s));
    return true;
  }  // end of setLinearSolver

  bool NonLinearEvolutionProblemImplementation<false>::integrate(
      const mfem::Vector& u,
      const IntegrationType it,
      const std::optional<real> odt) {
    if (this->mgis_integrator == nullptr) {
      return true;
    }
    // the behaviours are integrated on the host, unless they are all
    // integrated on the device from gradients computed by kernels (or
    // without gradients, for a prediction)
    if (!integratesOnDevice(*(this->mgis_integrator), this->batch_integration,
                            this->gradients_kernels, it)) {
      try {
        this->mgis_integrator->synchronizeOnHost();
      } catch (std::exception& e) {
        getErrorStream() << "integrate: " << e.what() << '\n';
        return false;
      } catch (...) {
        getErrorStream() << "integrate: unknown exception\n";
        return false;
      }
    }
    const auto dt = this->mgis_integrator->getTimeIncrement();
    if (odt.has_value()) {
      this->mgis_integrator->setTimeIncrement(*odt);
    }
    const auto& pu = this->Prolongate(u);
    const auto& fespace = this->getFiniteElementSpace();
    mfem::Array<int> vdofs;
    mfem::Vector ue;
    if (this->batch_integration) {
      if (!integrateByBatches(*(this->mgis_integrator), fespace, pu, it,
                              this->gradients_kernels)) {
        this->mgis_integrator->setTimeIncrement(dt);
        return false;
      }
    } else {
      for (size_type i = 0; i != fespace.GetNE(); ++i) {
        const auto& e = *(fespace.GetFE(i));
        auto& tr = *(fespace.GetElementTransformation(i));
        fespace.GetElementVDofs(i, vdofs);
        pu.GetSubVector(vdofs, ue);
        if (!this->mgis_integrator->integrate(e, tr, ue, it)) {
          this->mgis_integrator->setTimeIncrement(dt);
          return false;
        }
      }
    }
    if (it != IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR) {
      this->hasStiffnessOperatorsBeenComputed = true;
    }
    this->mgis_integrator->setTimeIncrement(dt);
    return true;
  }  // end of integrate

  void NonLinearEvolutionProblemImplementation<false>::
      markDegreesOfFreedomHandledByDirichletBoundaryConditions(
          std::vector<size_type> dofs) {
    auto tmp = mfem::Array<size_type>(dofs.data(), dofs.size());
    this->SetEssentialTrueDofs(tmp);
  }  // end of markDegreesOfFreedomHandledByDirichletBoundaryConditions

  NonLinearEvolutionProblemImplementation<
      false>::~NonLinearEvolutionProblemImplementation() = default;

}  // end of namespace mfem_mgis
