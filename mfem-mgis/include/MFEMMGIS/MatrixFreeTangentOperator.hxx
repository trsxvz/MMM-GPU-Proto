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
   * \brief interface of the operators applied without assembly which can
   * give the nodal blocks of their diagonal, the blocks which couple the
   * components of the unknowns at a node (used by the Jacobi preconditioner
   * by nodal blocks)
   */
  struct MFEM_MGIS_EXPORT NodalBlocksProvider {
    /*!
     * \brief compute the nodal blocks
     * \param[out] b: columns of the blocks, one vector of the true degrees of
     * freedom after the other: the entry `(c1, c2)` of the block of a node
     * is the entry of the component `c1` of this node in the vector of the
     * column `c2`
     */
    virtual void AssembleNodalBlocks(mfem::Vector& b) const = 0;
    //! \brief destructor
    virtual ~NodalBlocksProvider();
  };  // end of struct NodalBlocksProvider

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
  struct MFEM_MGIS_EXPORT MatrixFreeTangentOperator final
      : mfem::Operator,
        NodalBlocksProvider {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] mi: multi-material integrator
     */
    MatrixFreeTangentOperator(const FiniteElementDiscretization& fed,
                              MultiMaterialNonLinearIntegrator& mi);
    /*!
     * \brief apply the tangent operator
     * \param[in] x: true degrees of freedom
     * \param[out] y: action of the tangent operator on `x`
     */
    void Mult(const mfem::Vector& x, mfem::Vector& y) const override;
    /*!
     * \brief compute the diagonal of the tangent operator (see the
     * description of the class for a non-conforming mesh)
     * \param[out] d: diagonal, on the true degrees of freedom
     */
    void AssembleDiagonal(mfem::Vector& d) const override;
    /*!
     * \brief compute the nodal blocks of the tangent operator (see
     * `NodalBlocksProvider`), by kernels. As the diagonal, they are
     * approximated on a non-conforming mesh by `|Pᵀ|` applied to the blocks
     * of the local degrees of freedom.
     */
    void AssembleNodalBlocks(mfem::Vector& b) const override;
    //! \return the finite element space
    [[nodiscard]] const mfem::FiniteElementSpace& getFiniteElementSpace()
        const noexcept;
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
    //! \brief local columns of the nodal blocks (workspace)
    mutable mfem::Vector bl;
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
   * of freedom are constrained (see `mfem::ConstrainedOperator`, diagonal
   * policy `DIAG_ONE`: the rows and columns of the essential degrees of
   * freedom are the ones of the identity), which can give its assembled
   * version
   */
  struct MFEM_MGIS_EXPORT ConstrainedMatrixFreeTangentOperator final
      : mfem::ConstrainedOperator,
        AssembledOperatorProvider,
        NodalBlocksProvider {
    /*!
     * \brief constructor
     * \param[in] op: tangent operator applied without assembly, referenced:
     * it must outlive this operator
     * \param[in] essential_dofs: essential true degrees of freedom,
     * referenced (not copied) by `mfem::ConstrainedOperator`: the array must
     * outlive this operator and must not be modified
     * \param[in] f: function assembling the tangent operator, the essential
     * degrees of freedom being eliminated
     * \param[in] r: number of assemblies of the assembled tangent operator
     * (see `getAssembledOperatorRevision`), incremented by `f`
     */
    ConstrainedMatrixFreeTangentOperator(MatrixFreeTangentOperator& op,
                                         const mfem::Array<int>& essential_dofs,
                                         std::function<mfem::Operator&()> f,
                                         std::shared_ptr<const size_type> r);
    [[nodiscard]] mfem::Operator& getAssembledOperator() const override;
    [[nodiscard]] std::shared_ptr<const size_type>
    getAssembledOperatorRevision() const override;
    /*!
     * \brief compute the nodal blocks of the constrained operator (see
     * `NodalBlocksProvider`): the rows and the columns of the essential
     * degrees of freedom are the ones of the identity
     */
    void AssembleNodalBlocks(mfem::Vector& b) const override;
    //! \brief destructor
    ~ConstrainedMatrixFreeTangentOperator() override;

   private:
    //! \brief function assembling the tangent operator
    std::function<mfem::Operator&()> assemble;
    //! \brief number of assemblies of the assembled tangent operator
    std::shared_ptr<const size_type> revision;
    //! \brief tangent operator applied without assembly
    const MatrixFreeTangentOperator& tangent_operator;
    /*!
     * \brief entries of the columns of the nodal blocks set to zero, in the
     * rows and the columns of the essential degrees of freedom
     */
    mutable mfem::Array<int> null_entries;
    /*!
     * \brief entries of the columns of the nodal blocks set to one, on the
     * diagonal of the essential degrees of freedom
     */
    mutable mfem::Array<int> unit_entries;
  };  // end of struct ConstrainedMatrixFreeTangentOperator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MATRIXFREETANGENTOPERATOR_HXX */
