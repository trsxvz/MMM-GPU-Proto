/*!
 * \file   src/MatrixFreeTangentOperator.cxx
 * \brief  This file implements the `MatrixFreeTangentOperator` class
 * \date   10/10/2026
 */

#include <string>
#include <utility>
#include "mfem/fem/fespace.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/MatrixFreeTangentOperator.hxx"

namespace mfem_mgis {

  /*!
   * \return the finite element space of a finite element discretization
   * \param[in] fed: finite element discretization
   */
  static const mfem::FiniteElementSpace& getFiniteElementSpace(
      const FiniteElementDiscretization& fed) {
#ifdef MFEM_USE_MPI
    if (fed.describesAParallelComputation()) {
      return fed.getFiniteElementSpace<true>();
    }
#endif /* MFEM_USE_MPI */
    return fed.getFiniteElementSpace<false>();
  }  // end of getFiniteElementSpace

  /*!
   * \brief raises an exception on all the processes if a computation failed
   * on one of them, before the collective operations following it
   * \param[in] fed: finite element discretization
   * \param[in, out] ctx: execution context of the computation
   * \param[in] success: result of the computation on this process
   * \param[in] what: description of the computation
   */
  static void raiseOnFailureOnAnyProcess(const FiniteElementDiscretization& fed,
                                         Context& ctx,
                                         const bool success,
                                         const std::string& what) {
    if (isTrueOnAllProcesses(fed, success)) {
      return;
    }
    if (!success) {
      raise(what + ": " + ctx.getErrorMessage());
    }
    raise(what + ": failure on another process");
  }  // end of raiseOnFailureOnAnyProcess

  MatrixFreeTangentOperator::MatrixFreeTangentOperator(
      const FiniteElementDiscretization& fed,
      MultiMaterialNonLinearIntegrator& mi)
      : mfem::Operator(getFiniteElementSpace(fed).GetTrueVSize()),
        fe_discretization(fed),
        fe_space(getFiniteElementSpace(fed)),
        integrator(mi) {
    this->xl.UseDevice(true);
    this->yl.UseDevice(true);
  }  // end of MatrixFreeTangentOperator

  void MatrixFreeTangentOperator::Mult(const mfem::Vector& x,
                                       mfem::Vector& y) const {
    const auto* const P = this->fe_space.GetProlongationMatrix();
    if (P != nullptr) {
      this->xl.SetSize(P->Height());
      P->Mult(x, this->xl);
    }
    const auto& xl_ = (P != nullptr) ? this->xl : x;
    auto& yl_ = (P != nullptr) ? this->yl : y;
    yl_.SetSize(this->fe_space.GetVSize());
    yl_ = real{0};
    auto ctx = Context{};
    const auto success = this->integrator.addTangentOperatorActionByKernels(
        ctx, this->fe_space, xl_, yl_);
    raiseOnFailureOnAnyProcess(this->fe_discretization, ctx, success,
                               "action of the tangent operator");
    if (P != nullptr) {
      P->MultTranspose(yl_, y);
    }
  }  // end of Mult

  void MatrixFreeTangentOperator::AssembleDiagonal(mfem::Vector& d) const {
    const auto* const P = this->fe_space.GetProlongationMatrix();
    auto& dl = (P != nullptr) ? this->yl : d;
    dl.SetSize(this->fe_space.GetVSize());
    dl = real{0};
    auto ctx = Context{};
    const auto success = this->integrator.addTangentOperatorDiagonalByKernels(
        ctx, this->fe_space, dl);
    raiseOnFailureOnAnyProcess(this->fe_discretization, ctx, success,
                               "diagonal of the tangent operator");
    if (P != nullptr) {
      // the prolongation of a conforming space is boolean: the diagonal of
      // Pᵀ A P is Pᵀ diag(A) for the true degrees of freedom
      d.SetSize(this->Height());
      P->MultTranspose(dl, d);
    }
  }  // end of AssembleDiagonal

  MatrixFreeTangentOperator::~MatrixFreeTangentOperator() = default;

  AssembledOperatorProvider::~AssembledOperatorProvider() = default;

  ConstrainedMatrixFreeTangentOperator::ConstrainedMatrixFreeTangentOperator(
      MatrixFreeTangentOperator& op,
      const mfem::Array<int>& essential_dofs,
      std::function<mfem::Operator&()> f)
      : mfem::ConstrainedOperator(
            &op, essential_dofs, false, mfem::Operator::DIAG_ONE),
        assemble(std::move(f)) {
  }  // end of ConstrainedMatrixFreeTangentOperator

  mfem::Operator& ConstrainedMatrixFreeTangentOperator::getAssembledOperator()
      const {
    return this->assemble();
  }  // end of getAssembledOperator

  ConstrainedMatrixFreeTangentOperator::
      ~ConstrainedMatrixFreeTangentOperator() = default;

}  // end of namespace mfem_mgis
