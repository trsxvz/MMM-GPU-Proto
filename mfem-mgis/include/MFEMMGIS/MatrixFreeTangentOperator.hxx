/*!
 * \file   include/MFEMMGIS/MatrixFreeTangentOperator.hxx
 * \brief  This file declares the `MatrixFreeTangentOperator` class
 * \date   10/10/2026
 */

#ifndef LIB_MFEM_MGIS_MATRIXFREETANGENTOPERATOR_HXX
#define LIB_MFEM_MGIS_MATRIXFREETANGENTOPERATOR_HXX

#include "mfem/linalg/operator.hpp"
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem {
  // forward declaration
  class FiniteElementSpace;
}  // end of namespace mfem

namespace mfem_mgis {

  // forward declarations
  struct FiniteElementDiscretization;
  struct MultiMaterialNonLinearIntegrator;

  /*!
   * \brief tangent operator of the behaviour integrators of mfem-mgis,
   * applied without assembly by kernels on the true degrees of freedom:
   * `y = Pᵀ (Σ Bᵀ K B w det J) P x`, `P` being the prolongation of the finite
   * element space (see
   * `MultiMaterialNonLinearIntegrator::addTangentOperatorActionByKernels`).
   *
   * Its diagonal is computed by kernels (`AssembleDiagonal`), as required
   * by the smoothers of MFEM (`mfem::OperatorJacobiSmoother`, ...).
   *
   * The essential degrees of freedom are not taken into account: see
   * `mfem::ConstrainedOperator`. The tangent operator blocks are the ones
   * of the last integration of the behaviours. A failure of the
   * computations by kernels on a process is reported by an exception on
   * all the processes.
   */
  struct MFEM_MGIS_EXPORT MatrixFreeTangentOperator final : mfem::Operator {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] mi: multi-material integrator
     */
    MatrixFreeTangentOperator(const FiniteElementDiscretization&,
                              MultiMaterialNonLinearIntegrator&);
    void Mult(const mfem::Vector&, mfem::Vector&) const override;
    void AssembleDiagonal(mfem::Vector&) const override;
    //! \brief destructor
    ~MatrixFreeTangentOperator() override;

   private:
    //! \brief finite element discretization
    const FiniteElementDiscretization& fe_discretization;
    //! \brief finite element space
    const mfem::FiniteElementSpace& fe_space;
    //! \brief multi-material integrator
    MultiMaterialNonLinearIntegrator& integrator;
    //! \brief local vector of the argument (workspace)
    mutable mfem::Vector xl;
    //! \brief local vector of the result (workspace)
    mutable mfem::Vector yl;
  };  // end of struct MatrixFreeTangentOperator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MATRIXFREETANGENTOPERATOR_HXX */
