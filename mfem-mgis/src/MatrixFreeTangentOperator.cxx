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
      : mfem::Operator(mfem_mgis::getFiniteElementSpace(fed).GetTrueVSize()),
        fe_discretization(fed),
        fe_space(mfem_mgis::getFiniteElementSpace(fed)),
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

  void MatrixFreeTangentOperator::AssembleNodalBlocks(mfem::Vector& b) const {
    const auto vdim = this->fe_space.GetVDim();
    const auto nl = this->fe_space.GetVSize();
    const auto* const P = this->fe_space.GetProlongationMatrix();
    auto& bl_ = (P != nullptr) ? this->bl : b;
    bl_.UseDevice(true);
    bl_.SetSize(vdim * nl);
    bl_ = real{0};
    auto ctx = Context{};
    const auto success =
        this->integrator.addTangentOperatorNodalBlocksByKernels(
            ctx, this->fe_space, bl_);
    raiseOnFailureOnAnyProcess(this->fe_discretization, ctx, success,
                               "nodal blocks of the tangent operator");
    if (P == nullptr) {
      return;
    }
    // one column of the blocks after the other (see `AssembleDiagonal`)
    const auto nt = this->Height();
    b.UseDevice(true);
    b.SetSize(vdim * nt);
    for (int c = 0; c != vdim; ++c) {
      auto bc = mfem::Vector(b, c * nt, nt);
      const auto blc = mfem::Vector(bl_, c * nl, nl);
      if (this->fe_space.Conforming()) {
        P->MultTranspose(blc, bc);
      } else {
        P->AbsMultTranspose(blc, bc);
      }
    }
  }  // end of AssembleNodalBlocks

  const mfem::FiniteElementSpace&
  MatrixFreeTangentOperator::getFiniteElementSpace() const noexcept {
    return this->fe_space;
  }  // end of getFiniteElementSpace

  MatrixFreeTangentOperator::~MatrixFreeTangentOperator() = default;

  NodalBlocksProvider::~NodalBlocksProvider() = default;

  AssembledOperatorProvider::~AssembledOperatorProvider() = default;

  ConstrainedMatrixFreeTangentOperator::ConstrainedMatrixFreeTangentOperator(
      MatrixFreeTangentOperator& op,
      const mfem::Array<int>& essential_dofs,
      std::function<mfem::Operator&()> f,
      std::shared_ptr<const size_type> r)
      : mfem::ConstrainedOperator(
            &op, essential_dofs, false, mfem::Operator::DIAG_ONE),
        assemble(std::move(f)),
        revision(std::move(r)),
        tangent_operator(op) {
  }  // end of ConstrainedMatrixFreeTangentOperator

  void ConstrainedMatrixFreeTangentOperator::AssembleNodalBlocks(
      mfem::Vector& b) const {
    this->tangent_operator.AssembleNodalBlocks(b);
    const auto ne = this->constraint_list.Size();
    if (ne == 0) {
      return;
    }
    const auto& fes = this->tangent_operator.getFiniteElementSpace();
    const auto vdim = fes.GetVDim();
    if (this->unit_entries.Size() != ne) {
      // entries of the rows and of the columns of the essential degrees of
      // freedom, computed once (the essential degrees of freedom do not
      // change)
      const auto nt = this->Height();
      const auto nb = nt / vdim;
      const auto by_nodes = fes.GetOrdering() == mfem::Ordering::byNODES;
      const auto unknown = [nb, vdim, by_nodes](const int node, const int c) {
        return by_nodes ? c * nb + node : node * vdim + c;
      };
      const auto* const dofs = this->constraint_list.HostRead();
      this->null_entries.SetSize(2 * ne * vdim);
      this->unit_entries.SetSize(ne);
      auto* const null_entries_ = this->null_entries.HostWrite();
      auto* const unit_entries_ = this->unit_entries.HostWrite();
      for (int i = 0; i != ne; ++i) {
        const auto t = dofs[i];
        const auto node = by_nodes ? t % nb : t / vdim;
        const auto c = by_nodes ? t / nb : t % vdim;
        for (int d = 0; d != vdim; ++d) {
          // row of the degree of freedom, column d of the block
          null_entries_[2 * (vdim * i + d)] = d * nt + t;
          // column of the degree of freedom, row d of the block
          null_entries_[2 * (vdim * i + d) + 1] = c * nt + unknown(node, d);
        }
        unit_entries_[i] = c * nt + t;
      }
      this->null_entries.UseDevice(true);
      this->unit_entries.UseDevice(true);
    }
    b.SetSubVector(this->null_entries, real{0});
    b.SetSubVector(this->unit_entries, real{1});
  }  // end of AssembleNodalBlocks

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
