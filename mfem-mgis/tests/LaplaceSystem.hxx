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
     */
    LaplaceSystem(mfem::Mesh& smesh, mfem::FiniteElementCollection& fec)
        : mesh(makeMesh(smesh)), fespace(mesh.get(), &fec), a(&fespace) {
      auto ess_bdr = mfem::Array<int>(this->mesh->bdr_attributes.Max());
      ess_bdr = 1;
      auto ess_tdofs = mfem::Array<int>{};
      this->fespace.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
      this->a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
      this->a.SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
      this->a.Assemble();
      auto one = mfem::ConstantCoefficient(1.);
      auto b = LinearForm(&(this->fespace));
      b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
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
      auto r = mfem::Vector(this->B.Size());
      r.UseDevice(true);
      this->A.Mult(x, r);
      r -= this->B;
      return std::sqrt(this->dot(r, r) / this->dot(this->B, this->B));
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
