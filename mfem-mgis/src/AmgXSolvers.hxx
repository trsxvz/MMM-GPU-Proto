/*!
 * \file   src/AmgXSolvers.hxx
 * \brief  This file declares the solvers of the AmgX library of NVIDIA
 * built by the factory of linear solvers (internal header)
 * \date   10/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_AMGXSOLVERS_HXX
#define LIB_MFEM_MGIS_SRC_AMGXSOLVERS_HXX

#ifdef MFEM_MGIS_HAVE_AMGX

#include <memory>
#include <string>
#include <optional>
// the headers of AmgX compute modes by arithmetic between enumerations,
// deprecated in C++20
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-enum-enum-conversion"
#endif
#include <amgx_c.h>
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#include "mfem/config/config.hpp"
#include "mfem/general/array.hpp"
#include "mfem/linalg/vector.hpp"
#include "mfem/linalg/solvers.hpp"
#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief precision of the computations of AmgX (the mixed precision of
   * AmgX, matrices in single precision and vectors in double precision, is
   * not supported by AmgX 2.5 with CUDA 10.1 or later)
   */
  enum struct AmgXPrecision {
    //! \brief vectors and matrices in double precision (mode `dDDI`)
    DOUBLE,
    /*!
     * \brief vectors and matrices in single precision (mode `dFFI`), the
     * vectors exchanged with MFEM being converted
     */
    SINGLE
  };  // end of AmgXPrecision

  //! \brief options of a solver of AmgX
  struct AmgXSolverOptions {
    //! \brief precision of the computations of AmgX
    AmgXPrecision precision = AmgXPrecision::DOUBLE;
    /*!
     * \brief if true, the solver is an iterative solver whose number of
     * iterations is given by `getNumberOfIterations`
     */
    bool iterative = false;
    /*!
     * \brief if true, the right hand side is multiplied by a power of two
     * giving to the maximum of the absolute values of its components a
     * value between one and two, and the solution by its inverse: exact
     * operations, for a configuration whose result does not depend on the
     * scale of the right hand side (no absolute tolerance). The floors of
     * the convergence criteria of AmgX and the range of the single
     * precision then do not depend on the scale of the right hand side;
     * the statistics printed by AmgX are the ones of the scaled system.
     */
    bool scale_right_hand_side = false;
  };  // end of AmgXSolverOptions

  /*!
   * \brief solver of the AmgX library, built on its C interface. AmgX
   * computes on the current GPU of the process, whatever the device used by
   * MFEM. The resources of AmgX are shared by all the solvers of the
   * process, created on `MPI_COMM_WORLD` for the first solver and destroyed
   * with the last one: the communicator of a parallel solver must contain
   * all the processes, in the same order, or only one process. Each process
   * must own rows of the matrix. The local rows of the matrix, with the
   * global indices of their columns, are uploaded from the host at each call
   * to `SetOperator`. The resolutions start from a zero initial guess.
   *
   * The solvers must be destroyed before `MPI_Finalize` (AmgX frees its
   * communicators). With a MFEM compiled with AmgX, which contains its own
   * copy of AmgX, the library of AmgX given to mfem-mgis must precede MFEM
   * in the link line (as the target exported by mfem-mgis does, see
   * `src/CMakeLists.txt`), and a program exiting with a solver alive may
   * abort or hang (static objects of the two copies).
   */
  struct AmgXSolver final : public mfem::Solver {
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] c: communicator (`MPI_COMM_SELF` for a sequential
     * resolution)
     * \param[in] config: configuration of AmgX
     * \param[in] o: options
     */
    AmgXSolver(MPI_Comm c,
               const std::string& config,
               const AmgXSolverOptions& o);
#else  /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] config: configuration of AmgX
     * \param[in] o: options
     */
    AmgXSolver(const std::string& config, const AmgXSolverOptions& o);
#endif /* MFEM_USE_MPI */
    /*!
     * \brief uploads the local rows of the matrix from the host and sets up
     * the solver
     * \param[in] op: matrix (`mfem::SparseMatrix` or
     * `mfem::HypreParMatrix`), an exception being thrown for the other
     * operators
     */
    void SetOperator(const mfem::Operator& op) override;
    /*!
     * \brief solves the linear system, from a zero initial guess
     * \param[in] b: right hand side
     * \param[out] x: solution
     */
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override;
    /*!
     * \return the number of iterations of the last resolution of an
     * iterative solver, if a resolution was made
     */
    [[nodiscard]] std::optional<int> getNumberOfIterations() const noexcept;
    /*!
     * \return if the last resolution of an iterative solver succeeded
     * (status of AmgX), true for the other solvers
     */
    [[nodiscard]] bool hasConverged() const noexcept;
    //! \brief destructor
    ~AmgXSolver() override;
    // the handles of AmgX are owned by the solver
    AmgXSolver(const AmgXSolver&) = delete;
    AmgXSolver& operator=(const AmgXSolver&) = delete;

   private:
    //! \brief creates the objects of AmgX
    void create(const std::string&);
    //! \brief releases the objects of AmgX
    void release() noexcept;
    //! \brief uploads a sequential matrix
    void upload(const mfem::SparseMatrix&);
#ifdef MFEM_USE_MPI
    //! \brief uploads a parallel matrix
    void upload(const mfem::HypreParMatrix&);
    //! \brief communicator (duplicated)
    MPI_Comm comm = MPI_COMM_NULL;
#endif /* MFEM_USE_MPI */
    //! \brief mode of AmgX
    const AMGX_Mode mode;
    //! \brief if true, the solver is an iterative solver
    const bool iterative;
    //! \brief if true, the right hand side is scaled
    const bool scale_right_hand_side;
    //! \brief configuration
    AMGX_config_handle configuration = nullptr;
    //! \brief resources, shared by the solvers of the process
    AMGX_resources_handle resources = nullptr;
    //! \brief matrix
    AMGX_matrix_handle matrix = nullptr;
    //! \brief right hand side
    AMGX_vector_handle rhs = nullptr;
    //! \brief solution
    AMGX_vector_handle solution = nullptr;
    //! \brief solver
    AMGX_solver_handle solver = nullptr;
    //! \brief if true, the solver is set up
    bool setup = false;
    //! \brief status of the last resolution
    mutable std::optional<AMGX_SOLVE_STATUS> status;
    //! \brief number of iterations of the last resolution
    mutable std::optional<int> iterations;
    //! \brief scaled right hand side in double precision
    mutable mfem::Vector scaled_rhs;
    //! \brief workspace of the reduction giving the scaling factor
    mutable mfem::Array<double> reduction_workspace{
        mfem::MemoryType::HOST_PINNED};
    //! \brief right hand side in single precision
    mutable mfem::Array<float> single_rhs;
    //! \brief solution in single precision
    mutable mfem::Array<float> single_solution;
    //! \brief zero initial guess
    mutable mfem::Vector zeros;
    //! \brief zero initial guess in single precision
    mutable mfem::Array<float> single_zeros;
  };  // end of AmgXSolver

  /*!
   * \return cycles of the algebraic multigrid of AmgX, as a preconditioner
   * or an approximate linear solver (`AmgX` in the factory of linear
   * solvers).
   *
   * The options are:
   * - `Configuration`: configuration of AmgX, as a string;
   * - `ConfigurationFile`: file containing the configuration of AmgX (the
   *   configurations distributed with AmgX are in `lib/configs`);
   * - `VerbosityLevel`: if strictly positive, AmgX prints the statistics of
   *   the default configuration;
   * - `AmgXParameters`: parameters of AmgX, with their names in AmgX
   *   (`strength_threshold`, `interpolator`, ...), which complete or
   *   override the default configuration. A parameter without scope is
   *   given to the multigrid (scope `amg`), a parameter `scope:name` as is
   *   (the smoother has the scope `smoother`);
   * - `Precision`: `Double` (default) or `Single` (computations of AmgX in
   *   single precision).
   * `Configuration` and `ConfigurationFile` are exclusive; `VerbosityLevel`
   * and `AmgXParameters` only apply to the default configuration and are
   * refused with a configuration given.
   * The default configuration is the one of MFEM: one V-cycle of the
   * classical algebraic multigrid with a L1-Jacobi smoother (one sweep
   * before and one after the coarse correction), the D2 interpolation and
   * a strength threshold of 0.25. In single precision, the right hand side
   * of the default configuration is scaled by a power of two (see
   * `AmgXSolverOptions`), so that its values stay in the range of the single
   * precision: a tolerance given by `AmgXParameters` applies to the scaled
   * right hand side in single precision, to the right hand side as given in
   * double precision. A configuration given is used without scaling (its
   * criteria apply to the right hand side as given).
   *
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] params: options
   */
  template <bool parallel>
  [[nodiscard]] std::unique_ptr<AmgXSolver> makeAmgXPreconditioner(
      Context& ctx,
      FiniteElementSpace<parallel>& fespace,
      const Parameters& params) noexcept;
  /*!
   * \return the conjugate gradient of AmgX preconditioned by one V-cycle of
   * its algebraic multigrid in the default configuration of MFEM
   * (`AmgXPCG` in the factory of linear solvers).
   *
   * The options are the ones of the iterative solvers (`VerbosityLevel`,
   * `RelativeTolerance`, `AbsoluteTolerance`, `MaximumNumberOfIterations`,
   * by default 0, 1e-12, 0 and 150): the resolution converges when the
   * euclidean norm of the residual is lower than the absolute tolerance or
   * than the relative tolerance times the norm of the initial residual (the
   * norm of the right hand side, the initial guess being zero). The
   * conjugate gradient of MFEM checks instead the norm of the
   * preconditioned residual. The options `AmgXParameters` and `Precision` are
   * the ones of the `AmgX` solver; the parameters of the conjugate gradient
   * (scope `main`) are refused.
   *
   * AmgX also reports a convergence when the norm of the residual is lower
   * than the maximum of `e` times the norm of the initial residual and of
   * 1e-20, `e` being 1e-12 in double precision and 1e-6 in single
   * precision. A strictly positive relative tolerance must be at least `e`.
   * Without absolute tolerance, the right hand side is scaled by a power of
   * two (see `AmgXSolverOptions`): only the floor `e` times the norm of the
   * initial residual remains. With an absolute tolerance, the right hand
   * side is not scaled: the absolute tolerance must be at least 1e-20, and
   * an absolute tolerance alone (relative tolerance of 0) is reached or the
   * norm of the residual is lower than `e` times the norm of the initial
   * residual.
   * In single precision, the right hand side is always scaled: a relative
   * tolerance is required, and an absolute tolerance is refused.
   *
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] params: options
   */
  template <bool parallel>
  [[nodiscard]] std::unique_ptr<AmgXSolver> makeAmgXConjugateGradientSolver(
      Context& ctx,
      FiniteElementSpace<parallel>& fespace,
      const Parameters& params) noexcept;

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_MGIS_HAVE_AMGX */

#endif /* LIB_MFEM_MGIS_SRC_AMGXSOLVERS_HXX */
