/*!
 * \file   LaplaceSystem.hxx
 * \brief  Laplace equation on a cube, used by the tests of the linear
 * solvers
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_TESTS_LAPLACESYSTEM_HXX
#define LIB_MFEM_MGIS_TESTS_LAPLACESYSTEM_HXX

#include <cmath>
#include <memory>
#include <type_traits>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::tests {

  //! \brief Laplace equation on a cube, sequential or parallel
  template <bool is_parallel>
  struct LaplaceSystem {
    //! \brief type of the mesh
    using Mesh = std::conditional_t<is_parallel, mfem::ParMesh, mfem::Mesh>;
    //! \brief type of the bilinear form
    using BilinearForm = std::
        conditional_t<is_parallel, mfem::ParBilinearForm, mfem::BilinearForm>;
    //! \brief type of the linear form
    using LinearForm =
        std::conditional_t<is_parallel, mfem::ParLinearForm, mfem::LinearForm>;
    //! \brief type of the grid function
    using GridFunction = std::
        conditional_t<is_parallel, mfem::ParGridFunction, mfem::GridFunction>;
    //! \brief type of the matrix
    using Matrix = std::
        conditional_t<is_parallel, mfem::HypreParMatrix, mfem::SparseMatrix>;
    /*!
     * \param[in] smesh: sequential mesh (partitioned in parallel)
     * \param[in] fec: finite element collection
     * \param[in] dirichlet: if true, the unknown is null on the boundary and
     * the source is uniform; otherwise, the flux is null on the boundary (the
     * kernel of the operator is the constant vector) and the source, whose
     * integral is null, is `x - 0.5` (on the unit cube)
     */
    LaplaceSystem(mfem::Mesh& smesh,
                  mfem::FiniteElementCollection& fec,
                  const bool dirichlet = true)
        : mesh(makeMesh(smesh)), fespace(mesh.get(), &fec), a(&fespace) {
      auto ess_bdr = mfem::Array<int>(this->mesh->bdr_attributes.Max());
      ess_bdr = dirichlet ? 1 : 0;
      auto ess_tdofs = mfem::Array<int>{};
      this->fespace.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
      this->a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
      this->a.SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
      this->a.Assemble();
      auto one = mfem::ConstantCoefficient(1.);
      auto source = mfem::FunctionCoefficient(
          [](const mfem::Vector& p) { return p[0] - mfem::real_t{0.5}; });
      auto b = LinearForm(&(this->fespace));
      if (dirichlet) {
        b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
      } else {
        b.AddDomainIntegrator(new mfem::DomainLFIntegrator(source));
      }
      b.Assemble();
      auto x = GridFunction(&(this->fespace));
      x = 0.;
      // the right hand side given by MFEM may refer to the memory of the
      // linear form: it is copied
      auto X = mfem::Vector{};
      auto rhs = mfem::Vector{};
      this->a.FormLinearSystem(ess_tdofs, x, b, this->A, X, rhs);
      this->B = rhs;
    }
    //! \return the norm of the residual relative to the right hand side
    double getRelativeResidual(const mfem::Vector& x) {
      return this->getRelativeResidual(this->A, x);
    }
    /*!
     * \return the norm of the residual of the system of another matrix, with
     * the same right hand side, relative to the right hand side
     * \param[in] m: matrix
     * \param[in] x: solution
     */
    double getRelativeResidual(const Matrix& m, const mfem::Vector& x) {
      auto r = mfem::Vector(this->B.Size());
      r.UseDevice(true);
      m.Mult(x, r);
      r -= this->B;
      return std::sqrt(this->dot(r, r) / this->dot(this->B, this->B));
    }
    /*!
     * \return the mean of a vector relative to its norm, |sum_i x_i| /
     * (sqrt(n) ||x||), null for a vector orthogonal to the constant vector
     * \param[in] x: vector
     */
    double getRelativeMean(const mfem::Vector& x) {
      auto one = mfem::Vector(x.Size());
      one.UseDevice(true);
      one = 1.;
      return std::abs(this->dot(one, x)) /
             std::sqrt(this->dot(one, one) * this->dot(x, x));
    }
    //! \brief mesh (copy of the sequential mesh, or partition of it)
    std::unique_ptr<Mesh> mesh;
    //! \brief finite element space
    FiniteElementSpace<is_parallel> fespace;
    //! \brief bilinear form
    BilinearForm a;
    //! \brief matrix
    Matrix A;
    //! \brief right hand side
    mfem::Vector B;

   private:
    //! \return the mesh of the system
    static std::unique_ptr<Mesh> makeMesh(mfem::Mesh& smesh) {
      if constexpr (is_parallel) {
        return std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, smesh);
      } else {
        return std::make_unique<mfem::Mesh>(smesh);
      }
    }
    //! \return the scalar product of two vectors
    double dot(const mfem::Vector& u, const mfem::Vector& v) {
      if constexpr (is_parallel) {
        return mfem::InnerProduct(this->fespace.GetComm(), u, v);
      } else {
        return mfem::InnerProduct(u, v);
      }
    }
  };

}  // end of namespace mfem_mgis::tests

#endif /* LIB_MFEM_MGIS_TESTS_LAPLACESYSTEM_HXX */
