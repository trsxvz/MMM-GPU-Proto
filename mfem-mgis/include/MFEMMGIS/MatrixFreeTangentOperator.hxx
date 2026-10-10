/*!
 * \file   include/MFEMMGIS/MatrixFreeTangentOperator.hxx
 * \brief  This file declares the `MatrixFreeTangentOperator` class
 * \date   10/10/2026
 */

#ifndef LIB_MFEM_MGIS_MATRIXFREETANGENTOPERATOR_HXX
#define LIB_MFEM_MGIS_MATRIXFREETANGENTOPERATOR_HXX

#include <memory>
#include <functional>
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
   * by the smoothers of MFEM (`mfem::OperatorJacobiSmoother`, ...). On a
   * non-conforming mesh, it is approximated by `|Pᵀ| diag(Σ Bᵀ K B w det J)`,
   * as by `mfem::BilinearForm::AssembleDiagonal`.
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

  /*!
   * \brief interface of the operators applied without assembly which can
   * give an assembled version of themselves, used by the preconditioners
   * requiring a matrix (see the `AssembledOperator` preconditioner of the
   * linear solvers)
   */
  struct MFEM_MGIS_EXPORT AssembledOperatorProvider {
    /*!
     * \return the assembled operator, assembled at each call. It stays valid
     * until the next assembly, by this provider or by another one of the
     * same problem: see `getAssembledOperatorRevision`.
     */
    [[nodiscard]] virtual mfem::Operator& getAssembledOperator() const = 0;
    /*!
     * \return the number of assemblies of the assembled operator of the
     * problem, shared by all its providers and incremented at each assembly:
     * an assembled operator stays valid as long as this number keeps the
     * value it had just after its assembly. The counter outlives the
     * providers.
     */
    [[nodiscard]] virtual std::shared_ptr<const size_type>
    getAssembledOperatorRevision() const = 0;
    //! \brief destructor
    virtual ~AssembledOperatorProvider();
  };  // end of struct AssembledOperatorProvider

  /*!
   * \brief tangent operator applied without assembly whose essential degrees
   * of freedom are constrained (see `mfem::ConstrainedOperator`), which can
   * give its assembled version
   */
  struct MFEM_MGIS_EXPORT ConstrainedMatrixFreeTangentOperator final
      : mfem::ConstrainedOperator,
        AssembledOperatorProvider {
    /*!
     * \brief constructor
     * \param[in] op: tangent operator applied without assembly
     * \param[in] essential_dofs: essential true degrees of freedom
     * \param[in] f: function assembling the tangent operator, the essential
     * degrees of freedom being eliminated
     * \param[in] r: number of assemblies of the assembled tangent operator
     * (see `getAssembledOperatorRevision`), incremented by `f`
     */
    ConstrainedMatrixFreeTangentOperator(MatrixFreeTangentOperator&,
                                         const mfem::Array<int>&,
                                         std::function<mfem::Operator&()>,
                                         std::shared_ptr<const size_type>);
    [[nodiscard]] mfem::Operator& getAssembledOperator() const override;
    [[nodiscard]] std::shared_ptr<const size_type>
    getAssembledOperatorRevision() const override;
    //! \brief destructor
    ~ConstrainedMatrixFreeTangentOperator() override;

   private:
    //! \brief function assembling the tangent operator
    std::function<mfem::Operator&()> assemble;
    //! \brief number of assemblies of the assembled tangent operator
    std::shared_ptr<const size_type> revision;
  };  // end of struct ConstrainedMatrixFreeTangentOperator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MATRIXFREETANGENTOPERATOR_HXX */
