/*!
 * \file   src/LinearSolverModules.hxx
 * \brief  This file declares the modules registering the linear solvers
 * and the preconditioners of each library in the factories, and the
 * functions shared by these modules (internal header)
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_LINEARSOLVERMODULES_HXX
#define LIB_MFEM_MGIS_SRC_LINEARSOLVERMODULES_HXX

#include <memory>
#include <string>
#include <utility>
#include <optional>
#include <string_view>
#include "MGIS/Contract.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/LinearSolverDescription.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/PreconditionerFactory.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief Each library providing linear solvers or preconditioners has a
   * module which registers them in the factories. A module is a pair of
   * functions, called by the constructors of the factories when the library
   * is available:
   *
   * - `declare<Library>LinearSolvers(LinearSolverFactory<parallel>&)`;
   * - `declare<Library>Preconditioners(PreconditionerFactory<parallel>&)`.
   *
   * Adding a library consists of writing its module and calling it in
   * `declareLinearSolvers` and `declarePreconditioners`.
   */

  /*!
   * \brief register the linear solvers of all the available libraries
   * \param[in, out] f: factory
   */
  template <bool parallel>
  void declareLinearSolvers(LinearSolverFactory<parallel>& f);
  /*!
   * \brief register the preconditioners of all the available libraries
   * \param[in, out] f: factory
   */
  template <bool parallel>
  void declarePreconditioners(PreconditionerFactory<parallel>& f);

  /*!
   * \brief register the solvers of MFEM (Krylov solvers, direct solvers
   * UMFPACK and KLU of SuiteSparse in sequential)
   * \param[in, out] f: factory
   */
  template <bool parallel>
  void declareMFEMLinearSolvers(LinearSolverFactory<parallel>& f);
  /*!
   * \brief register the preconditioners of MFEM (Jacobi and Chebyshev
   * smoothers, smoothers and block incomplete LU factorization of the sparse
   * matrices, direct solvers UMFPACK and KLU of SuiteSparse, orthogonalizing
   * wrapper) and of mfem-mgis (`AssembledOperator`)
   * \param[in, out] f: factory
   */
  template <bool parallel>
  void declareMFEMPreconditioners(PreconditionerFactory<parallel>& f);
#ifdef MFEM_USE_MPI
  /*!
   * \brief register the solvers of hypre
   * \param[in, out] f: factory
   */
  void declareHypreLinearSolvers(LinearSolverFactory<true>& f);
  /*!
   * \brief register the preconditioners of hypre
   * \param[in, out] f: factory
   */
  void declareHyprePreconditioners(PreconditionerFactory<true>& f);
#endif /* MFEM_USE_MPI */
#if defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS)
  /*!
   * \brief register the solver of MUMPS
   * \param[in, out] f: factory
   */
  void declareMUMPSLinearSolvers(LinearSolverFactory<true>& f);
  /*!
   * \brief register the preconditioner of MUMPS (factorization applied once)
   * \param[in, out] f: factory
   */
  void declareMUMPSPreconditioners(PreconditionerFactory<true>& f);
#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS) */
#if defined(MFEM_USE_MPI) && defined(MFEM_USE_CUDSS)
  /*!
   * \brief register the solver of cuDSS
   * \param[in, out] f: factory
   */
  void declareCuDSSLinearSolvers(LinearSolverFactory<true>& f);
#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_CUDSS) */
#ifdef MFEM_USE_MKL_PARDISO
  /*!
   * \brief register the PARDISO direct solver of oneMKL (sequential)
   * \param[in, out] f: factory
   */
  void declareMKLPardisoLinearSolvers(LinearSolverFactory<false>& f);
#endif /* MFEM_USE_MKL_PARDISO */
#if defined(MFEM_USE_MPI) && defined(MFEM_USE_MKL_CPARDISO)
  /*!
   * \brief register the Cluster Sparse Solver of oneMKL (parallel)
   * \param[in, out] f: factory
   */
  void declareMKLClusterSparseLinearSolvers(LinearSolverFactory<true>& f);
#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_MKL_CPARDISO) */
#ifdef MFEM_MGIS_HAVE_MAGMA_SPARSE
  /*!
   * \brief register the Krylov solvers of MAGMA (sequential, one GPU)
   * \param[in, out] f: factory
   */
  void declareMAGMALinearSolvers(LinearSolverFactory<false>& f);
  /*!
   * \brief register the preconditioners of MAGMA (applied to the diagonal
   * block of each process in parallel)
   * \param[in, out] f: factory
   */
  template <bool parallel>
  void declareMAGMAPreconditioners(PreconditionerFactory<parallel>& f);
#endif /* MFEM_MGIS_HAVE_MAGMA_SPARSE */
  /*!
   * \brief register the Jacobi preconditioner by nodal blocks, inverted by
   * the batched dense linear algebra of MFEM
   * \param[in, out] f: factory
   */
  template <bool parallel>
  void declareNodalBlockJacobiPreconditioner(
      PreconditionerFactory<parallel>& f);
#if defined(MFEM_USE_MPI) && defined(MFEM_USE_SUPERLU)
  /*!
   * \brief register the direct solver of SuperLU_DIST
   * \param[in, out] f: factory
   */
  void declareSuperLULinearSolvers(LinearSolverFactory<true>& f);
#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_SUPERLU) */
#if defined(MFEM_USE_MPI) && defined(MFEM_USE_STRUMPACK)
  /*!
   * \brief register the direct solvers of STRUMPACK
   * \param[in, out] f: factory
   */
  void declareSTRUMPACKLinearSolvers(LinearSolverFactory<true>& f);
  /*!
   * \brief register the preconditioner of STRUMPACK
   * \param[in, out] f: factory
   */
  void declareSTRUMPACKPreconditioners(PreconditionerFactory<true>& f);
#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_STRUMPACK) */
#if defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC)
  /*!
   * \brief register the Krylov solver of PETSc
   * \param[in, out] f: factory
   */
  void declarePETScLinearSolvers(LinearSolverFactory<true>& f);
  /*!
   * \brief register the preconditioner of PETSc
   * \param[in, out] f: factory
   */
  void declarePETScPreconditioners(PreconditionerFactory<true>& f);
#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC) */
#ifdef MFEM_USE_GINKGO
  /*!
   * \brief register the solvers of Ginkgo (the distributed ones in parallel,
   * if Ginkgo is built with MPI)
   * \param[in, out] f: factory
   */
  template <bool parallel>
  void declareGinkgoLinearSolvers(LinearSolverFactory<parallel>& f);
  /*!
   * \brief register the preconditioners of Ginkgo
   * \param[in, out] f: factory
   */
  template <bool parallel>
  void declareGinkgoPreconditioners(PreconditionerFactory<parallel>& f);
#endif /* MFEM_USE_GINKGO */
#ifdef MFEM_MGIS_HAVE_AMGX
  /*!
   * \brief register the solvers of AmgX
   * \param[in, out] f: factory
   */
  template <bool parallel>
  void declareAmgXLinearSolvers(LinearSolverFactory<parallel>& f);
  /*!
   * \brief register the preconditioner of AmgX
   * \param[in, out] f: factory
   */
  template <bool parallel>
  void declareAmgXPreconditioners(PreconditionerFactory<parallel>& f);
#endif /* MFEM_MGIS_HAVE_AMGX */

  /*!
   * \return a message explaining why a linear solver or a preconditioner of
   * an optional library is not registered (the library is not available),
   * or an empty string if the name is not the one of such a solver
   * \param[in] n: name of the linear solver or of the preconditioner
   */
  [[nodiscard]] std::string getUnavailableLibraryHint(std::string_view n);

  /*!
   * \brief register a linear solver or a preconditioner in a factory, a
   * failure (two entries with the same name) being a contract violation
   * \param[in, out] f: factory
   * \param[in] n: name
   * \param[in] d: description
   * \param[in] g: generator
   */
  template <typename Factory, typename Generator>
  void declare(Factory& f,
               const std::string_view n,
               LinearSolverDescription d,
               Generator&& g) {
    auto ctx = Context{};
    if (!f.add(ctx, n, std::move(d), std::forward<Generator>(g))) {
      auto eh = mgis::ContractViolationHandler{};
      eh.registerErrorMessage(ctx.getErrorMessage().c_str());
    }
  }  // end of declare

  /*!
   * \return the preconditioner described by the given parameters: its name
   * (`Name`) and its options (`Options`), or a null pointer if the name is
   * `None` (without options) or on error
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] p: parameters of the preconditioner
   */
  template <bool parallel>
  [[nodiscard]] std::unique_ptr<LinearSolverPreconditioner>
  generatePreconditioner(Context& ctx,
                         FiniteElementSpace<parallel>& fespace,
                         const Parameters& p);

  /*!
   * \return the preconditioner given by the `Preconditioner` parameter of a
   * linear solver: a null pointer if this parameter is not given or selects
   * no preconditioner (name `None`), or an empty value on error. The
   * preconditioner must be accepted by the linear solver (see
   * `LinearSolverDescription::preconditioners`), and compute where the
   * execution settings place the linear solvers.
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] solver: name of the linear solver
   * \param[in] d: description of the linear solver
   * \param[in] params: parameters of the linear solver
   */
  template <bool parallel>
  [[nodiscard]] std::optional<std::unique_ptr<LinearSolverPreconditioner>>
  makePreconditioner(Context& ctx,
                     FiniteElementSpace<parallel>& fespace,
                     std::string_view solver,
                     const LinearSolverDescription& d,
                     const Parameters& params);

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_LINEARSOLVERMODULES_HXX */
