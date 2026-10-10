/*!
 * \file   src/MatrixFreeTangentOperator.cxx
 * \brief  This file implements the `MatrixFreeTangentOperator` class
 * \date   10/10/2026
 */

#include <utility>
#include "mfem/fem/fespace.hpp"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/MatrixFreeTangentOperator.hxx"

namespace mfem_mgis {

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
    if (P == nullptr) {
      return;
    }
    d.SetSize(this->Height());
    if (this->fe_space.Conforming()) {
      // the prolongation of a conforming space is boolean: the diagonal of
      // Pᵀ A P is Pᵀ diag(A)
      P->MultTranspose(dl, d);
    } else {
      // non-conforming mesh: |Pᵀ| diag(A), where |Pᵀ| has the absolute
      // values of the entries of Pᵀ, which approximates the diagonal of
      // Pᵀ A P, as `mfem::BilinearForm::AssembleDiagonal`
      P->AbsMultTranspose(dl, d);
    }
  }  // end of AssembleDiagonal

  MatrixFreeTangentOperator::~MatrixFreeTangentOperator() = default;

  AssembledOperatorProvider::~AssembledOperatorProvider() = default;

  ConstrainedMatrixFreeTangentOperator::ConstrainedMatrixFreeTangentOperator(
      MatrixFreeTangentOperator& op,
      const mfem::Array<int>& essential_dofs,
      std::function<mfem::Operator&()> f,
      std::shared_ptr<const size_type> r)
      : mfem::ConstrainedOperator(
            &op, essential_dofs, false, mfem::Operator::DIAG_ONE),
        assemble(std::move(f)),
        revision(std::move(r)) {
  }  // end of ConstrainedMatrixFreeTangentOperator

  mfem::Operator& ConstrainedMatrixFreeTangentOperator::getAssembledOperator()
      const {
    return this->assemble();
  }  // end of getAssembledOperator

  std::shared_ptr<const size_type>
  ConstrainedMatrixFreeTangentOperator::getAssembledOperatorRevision() const {
    return this->revision;
  }  // end of getAssembledOperatorRevision

  ConstrainedMatrixFreeTangentOperator::
      ~ConstrainedMatrixFreeTangentOperator() = default;

}  // end of namespace mfem_mgis
