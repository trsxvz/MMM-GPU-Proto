/*!
 * \file   MFEMMGIS/AbstractBoundaryCondition.hxx
 * \brief  This file declares the `AbstractBoundaryCondition` class
 * \author Thomas Helfer
 * \date   27/09/2024
 */

#ifndef LIB_MFEM_MGIS_ABSTRACTBOUNDARYCONDITION_HXX
#define LIB_MFEM_MGIS_ABSTRACTBOUNDARYCONDITION_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief base class for non essential boundary conditions
   */
  struct MFEM_MGIS_EXPORT AbstractBoundaryCondition {
#ifdef MFEM_USE_MPI
    /*!
     * \brief add the nonlinear form integrator describing the boundary
     * condition
     *
     * \param[in, out] ctx: execution context
     * \param[in, out] f: form
     * \param[in] u: current estimate of the solution at the end of the time
     * step
     * \return true on success
     */
    [[nodiscard]] virtual bool addNonlinearFormIntegrator(
        Context& ctx,
        NonlinearForm<true>& f,
        const mfem::Vector& u) noexcept = 0;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief add the nonlinear form integrator describing the boundary
     * condition
     *
     * \param[in, out] ctx: execution context
     * \param[in, out] f: form
     * \param[in] u: current estimate of the solution at the end of the time
     * step
     * \return true on success
     */
    [[nodiscard]] virtual bool addNonlinearFormIntegrator(
        Context& ctx,
        NonlinearForm<false>& f,
        const mfem::Vector& u) noexcept = 0;
#ifdef MFEM_USE_MPI
    /*!
     * \brief add the bilinear and linear form integrators
     * describing the boundary condition
     *
     * \param[in, out] ctx: execution context
     * \param[in, out] a: form computing the jacobian matrix
     * \param[in, out] b: form computing the right hand side
     * \param[in] u: solution at the beginning of the time step
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] virtual bool addLinearFormIntegrators(
        Context& ctx,
        BilinearForm<true>& a,
        LinearForm<true>& b,
        const mfem::Vector& u,
        const real t,
        const real dt) noexcept = 0;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief add the bilinear and linear form integrators
     * describing the boundary condition
     *
     * \param[in, out] ctx: execution context
     * \param[in, out] a: form computing the jacobian matrix
     * \param[in, out] b: form computing the right hand side
     * \param[in] u: solution at the beginning of the time step
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] virtual bool addLinearFormIntegrators(
        Context& ctx,
        BilinearForm<false>& a,
        LinearForm<false>& b,
        const mfem::Vector& u,
        const real t,
        const real dt) noexcept = 0;
    /*!
     * \brief method called at the beginning of each resolution
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] virtual bool setup(Context& ctx,
                                     const real t,
                                     const real dt) noexcept = 0;
    /*!
     * \return if the contribution of the boundary condition to the residual
     * is a vector independent of the unknowns multiplied by a scalar which
     * only depends on time, its contribution to the tangent operator being
     * null. With the inner forces computed by kernels (parameter
     * `InnerForcesComputation`), such a boundary condition is applied by
     * `addContributionToTheResidual` instead of an integrator of the non
     * linear form: if the non linear form contains no other integrator than
     * the one of mfem-mgis, the tangent operator can then be applied without
     * assembly. The default implementation returns false.
     */
    [[nodiscard]] virtual bool isIndependentOfTheUnknowns() const noexcept;
    /*!
     * \brief adds the contribution of a boundary condition independent of
     * the unknowns (see `isIndependentOfTheUnknowns`) to the residual, at the
     * end of the time step given to `setup`. The default implementation
     * reports an error.
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in, out] r: residual, in the layout of the local vectors of the
     * finite element space (one value per degree of freedom of the
     * process), on the device if a device is enabled
     */
    [[nodiscard]] virtual bool addContributionToTheResidual(
        Context& ctx, mfem::Vector& r) noexcept;
    //! \brief destructor
    virtual ~AbstractBoundaryCondition();
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACTBOUNDARYCONDITION_HXX */
