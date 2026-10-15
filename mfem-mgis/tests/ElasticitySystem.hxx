/*!
 * \file   ElasticitySystem.hxx
 * \brief  Cantilever beam in linear elasticity, used by the tests of the
 * linear solvers
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_TESTS_ELASTICITYSYSTEM_HXX
#define LIB_MFEM_MGIS_TESTS_ELASTICITYSYSTEM_HXX

#include <cmath>
#include <memory>
#include <type_traits>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::tests {

  /*!
   * \brief cantilever beam in linear elasticity (Lamé coefficients 1.5 and
   * 1), clamped on the face `x = 0` (boundary attribute 5 of the meshes
   * built by `mfem::Mesh::MakeCartesian3D`), loaded by a uniform vertical
   * force, sequential or parallel. The coefficients of the bilinear form
   * are members, which it keeps pointers to: the form may be assembled
   * again.
   */
  template <bool is_parallel>
  struct ElasticitySystem {
    //! \brief type of the mesh
    using Mesh = std::conditional_t<is_parallel, mfem::ParMesh, mfem::Mesh>;
    //! \brief type of the bilinear form
    using BilinearForm = std::
        conditional_t<is_parallel, mfem::ParBilinearForm, mfem::BilinearForm>;
    //! \brief type of the linear form
    using LinearForm =
        std::conditional_t<is_parallel, mfem::ParLinearForm, mfem::LinearForm>;
    //! \brief type of the grid functions
    using GridFunction = std::
        conditional_t<is_parallel, mfem::ParGridFunction, mfem::GridFunction>;
    //! \brief type of the matrix of the system
    using Matrix = std::
        conditional_t<is_parallel, mfem::HypreParMatrix, mfem::SparseMatrix>;
    /*!
     * \param[in] smesh: sequential mesh (partitioned in parallel)
     * \param[in] fec: finite element collection
     * \param[in] o: ordering of the unknowns
     * \param[in] partitioning: process of each element in parallel (a
     * partitioning of METIS if null)
     */
    ElasticitySystem(mfem::Mesh& smesh,
                     mfem::FiniteElementCollection& fec,
                     const mfem::Ordering::Type o = mfem::Ordering::byNODES,
                     const int* const partitioning = nullptr)
        : mesh(makeMesh(smesh, partitioning)),
          fespace(mesh.get(), &fec, 3, o),
          a(&fespace) {
      this->a.AddDomainIntegrator(
          new mfem::ElasticityIntegrator(this->lambda, this->mu));
      this->a.Assemble();
      auto f = mfem::Vector(3);
      f = 0.;
      f(2) = -1;
      auto fc = mfem::VectorConstantCoefficient(f);
      auto b = LinearForm(&(this->fespace));
      b.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(fc));
      b.Assemble();
      auto x = GridFunction(&(this->fespace));
      x = 0.;
      auto boundaries = mfem::Array<int>(this->mesh->bdr_attributes.Max());
      boundaries = 0;
      boundaries[4] = 1;
      auto ess_tdofs = mfem::Array<int>{};
      this->fespace.GetEssentialTrueDofs(boundaries, ess_tdofs);
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
    //! \brief first Lamé coefficient
    mfem::ConstantCoefficient lambda = mfem::ConstantCoefficient(1.5);
    //! \brief second Lamé coefficient
    mfem::ConstantCoefficient mu = mfem::ConstantCoefficient(1);
    //! \brief mesh
    std::unique_ptr<Mesh> mesh;
    //! \brief finite element space
    FiniteElementSpace<is_parallel> fespace;
    //! \brief bilinear form
    BilinearForm a;
    //! \brief matrix of the system
    Matrix A;
    //! \brief right hand side
    mfem::Vector B;

   private:
    /*!
     * \return the mesh, partitioned in parallel
     * \param[in] smesh: sequential mesh
     * \param[in] partitioning: process of each element, if not null
     */
    static std::unique_ptr<Mesh> makeMesh(mfem::Mesh& smesh,
                                          const int* const partitioning) {
      if constexpr (is_parallel) {
        return std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, smesh,
                                               partitioning);
      } else {
        static_cast<void>(partitioning);
        return std::make_unique<mfem::Mesh>(smesh);
      }
    }
    //! \return the scalar product of two vectors, on all the processes
    double dot(const mfem::Vector& u, const mfem::Vector& v) {
      if constexpr (is_parallel) {
        return mfem::InnerProduct(this->fespace.GetComm(), u, v);
      } else {
        return mfem::InnerProduct(u, v);
      }
    }
  };  // end of struct ElasticitySystem

}  // end of namespace mfem_mgis::tests

#endif /* LIB_MFEM_MGIS_TESTS_ELASTICITYSYSTEM_HXX */
