/*!
 * \file   src/TrilinosBackend.hxx
 * \brief  This file declares the interface of the library of mfem-mgis
 * built against Trilinos (`MFEMMGISTrilinos`, see src/CMakeLists.txt): the
 * solvers and the preconditioners of Trilinos are seen through standard
 * types only, so that the other sources of mfem-mgis, compiled by the
 * compiler of the project, never see Kokkos, Teuchos or Tpetra, whose
 * headers require the compiler of Kokkos (internal header)
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_TRILINOSBACKEND_HXX
#define LIB_MFEM_MGIS_SRC_TRILINOSBACKEND_HXX

#include <mpi.h>
#include <memory>
#include <string>
#include <optional>
#include <vector>
#include <utility>
#include <variant>
#include <functional>

namespace mfem_mgis::trilinos {

  // forward declaration
  struct Option;

  //! \brief list of options, converted to a `Teuchos::ParameterList`
  using Options = std::vector<std::pair<std::string, Option>>;

  //! \brief value of an option (a list of options is a sublist)
  struct Option {
    //! \brief value
    std::variant<bool, int, double, std::string, Options> value;
  };

  /*!
   * \brief view of a distributed square matrix in the format of hypre, on
   * the host: the diagonal block of the process, whose columns are the local
   * indices of the rows of the process, and the off-diagonal block, whose
   * columns are indices in `offd_global_columns`
   */
  struct DistributedMatrixView {
    //! \brief communicator
    MPI_Comm communicator = MPI_COMM_NULL;
    //! \brief number of rows of the matrix
    long long global_size = 0;
    //! \brief first row owned by the process
    long long first_row = 0;
    //! \brief number of rows owned by the process
    int local_size = 0;
    //! \brief offsets of the rows of the diagonal block
    const int* diag_offsets = nullptr;
    //! \brief columns of the diagonal block (local indices)
    const int* diag_columns = nullptr;
    //! \brief values of the diagonal block
    const double* diag_values = nullptr;
    //! \brief offsets of the rows of the off-diagonal block
    const int* offd_offsets = nullptr;
    //! \brief columns of the off-diagonal block (indices in
    //! `offd_global_columns`)
    const int* offd_columns = nullptr;
    //! \brief values of the off-diagonal block
    const double* offd_values = nullptr;
    //! \brief number of columns of the off-diagonal block
    int number_of_offd_columns = 0;
    //! \brief global indices of the columns of the off-diagonal block
    const long long* offd_global_columns = nullptr;
    /*!
     * \brief number of components of the unknowns. The unknowns ordered by
     * nodes (`ordered_by_nodes`), the components of all the nodes of a
     * process being stored one after the other, are seen by Trilinos ordered
     * by vector dimension, the components of a node being contiguous, as
     * the algebraic multigrid of MueLu expects: the rows and the columns of
     * each process are permuted, and the vectors accordingly.
     */
    int number_of_components = 1;
    //! \brief if true, the unknowns are ordered by nodes
    bool ordered_by_nodes = true;
    /*!
     * \brief if true, the arrays of the blocks are in the memory of the GPU
     * which the operator computes on (the global indices of the columns of
     * the off-diagonal block being in the memory of the host)
     */
    bool device = false;
  };

  /*!
   * \brief coordinates of the nodes of the unknowns of the process, used to
   * build the rigid body modes, the near null space of the algebraic
   * multigrid of the linear elasticity
   */
  struct NodesCoordinates {
    //! \brief space dimension
    int dimension = 0;
    /*!
     * \brief coordinates, component by component (`dimension` arrays of the
     * number of nodes of the process), only given (on all the processes)
     * when the number of components of the unknowns is the space dimension
     */
    std::vector<double> values;
  };

  //! \brief report of a resolution
  struct Report {
    //! \brief if true, the resolution converged
    bool converged = false;
    //! \brief number of iterations (zero for the direct solvers)
    int iterations = 0;
    //! \brief tolerance achieved, if reported by the solver
    double achieved_tolerance = 0;
    /*!
     * \brief message of the failure of the solver (loss of positive
     * definiteness, breakdown...), empty if it did not fail
     */
    std::string failure;
  };

  /*!
   * \brief solver or preconditioner of Trilinos applied to vectors of the
   * host. The methods report the errors of Trilinos by exceptions
   * (`std::runtime_error`).
   */
  struct Operator {
    /*!
     * \brief set the matrix, and set up the solver or the preconditioner
     * (the matrix is converted and the setup made again for each matrix)
     * \param[in] m: matrix
     */
    virtual void setMatrix(const DistributedMatrixView& m) = 0;
    /*!
     * \brief solve the system or apply the preconditioner
     * \param[in] b: right hand side (local part)
     * \param[in, out] x: solution (local part), initial guess of the
     * iterative solvers which use one
     */
    virtual void apply(const double* b, double* x) const = 0;
    /*!
     * \brief solve the system or apply the preconditioner, the vectors
     * being in the memory of the GPU which the operator computes on (see
     * `computesOnTheGPU`): they are copied on the GPU, the GPU being
     * synchronized before returning
     * \param[in] b: right hand side (local part)
     * \param[in, out] x: solution (local part), initial guess of the
     * iterative solvers which use one
     */
    virtual void applyOnDevice(const double* b, double* x) const = 0;
    //! \return true if the operator computes on the GPU
    [[nodiscard]] virtual bool computesOnTheGPU() const noexcept = 0;
    //! \return the report of the last resolution
    [[nodiscard]] virtual Report getReport() const = 0;
    /*!
     * \return the relative tolerance of an iterative solver, or a negative
     * value if the operator has none
     */
    [[nodiscard]] virtual double getRelativeTolerance() const;
    /*!
     * \brief set the relative tolerance of an iterative solver
     * \return false if the operator has no relative tolerance
     * \param[in] tol: relative tolerance
     */
    [[nodiscard]] virtual bool setRelativeTolerance(const double tol);
    //! \brief destructor
    virtual ~Operator();
  };

  /*!
   * \brief preconditioner of mfem-mgis, given to the Krylov solvers of
   * Belos: it owns its right hand side and its result, which Trilinos
   * writes and reads, in the memory of the GPU of Tpetra if `device` is
   * true (Tpetra and mfem-mgis computing on the GPU), and on the host
   * otherwise
   */
  struct ExternalPreconditioner {
    /*!
     * \brief function returning the right hand side (local part), to be
     * written before calling `apply`
     */
    std::function<double*()> getRightHandSide;
    //! \brief function applying the preconditioner, returning the result
    std::function<const double*()> apply;
    //! \brief if true, the vectors are in the memory of the GPU
    bool device = false;
    //! \return true if the preconditioner is given
    explicit operator bool() const noexcept {
      return static_cast<bool>(this->apply);
    }
  };

  /*!
   * \return a preconditioner of Ifpack2
   * \param[in] type: type of the preconditioner of Ifpack2 (`RELAXATION`,
   * `CHEBYSHEV`, `RILUK`, `ILUT`, `SCHWARZ`...)
   * \param[in] options: parameters of the preconditioner
   * \param[in] gpu: if true, computes on the GPU if Trilinos is built for
   * one (see `hasGPUNode`), on the host otherwise
   */
  [[nodiscard]] std::unique_ptr<Operator> makeIfpack2Preconditioner(
      const std::string& type, const Options& options, const bool gpu);
  /*!
   * \return the smoothed aggregation multigrid of MueLu (one V-cycle), the
   * number of equations per node being the number of components of the
   * unknowns, unless the parameter `number of equations` is given
   * \param[in] options: parameters of MueLu
   * \param[in] c: coordinates of the nodes: the rigid body modes are given
   * as the near null space when the number of components of the unknowns
   * is the space dimension, 2 or 3, and the number of equations per node
   * \param[in] gpu: if true, computes on the GPU if Trilinos is built for
   * one, on the host otherwise
   */
  [[nodiscard]] std::unique_ptr<Operator> makeMueLuPreconditioner(
      const Options& options, NodesCoordinates c, const bool gpu);
  /*!
   * \return the two level overlapping Schwarz preconditioner of FROSch
   * (`FROSch::TwoLevelPreconditioner`), the number of degrees of freedom
   * per node being the number of components of the unknowns: the null
   * space of the coarse space is the one of the linear elasticity when the
   * coordinates of the nodes are given, the constants otherwise, unless the
   * parameter `Null Space Type` is given
   * \param[in] options: parameters of FROSch
   * \param[in] c: coordinates of the nodes
   * \param[in] gpu: ignored: FROSch computes with the default node of
   * Tpetra, on the GPU if Trilinos is built for one (FROSch 17.2.1 runs some
   * loops in the default execution space of Kokkos, whatever the node)
   */
  [[nodiscard]] std::unique_ptr<Operator> makeFROSchPreconditioner(
      const Options& options, NodesCoordinates c, const bool gpu);
  /*!
   * \return a direct solver of Amesos2
   * \param[in] solver: name of the solver (`KLU2`, `Basker`, `Tacho`,
   * `SuperLU_DIST`...)
   * \param[in] options: parameters of the solver
   * \param[in] gpu: if true, computes on the GPU if Trilinos is built for
   * one, on the host otherwise, except Tacho, which computes in the
   * default execution space of Kokkos whatever the node: it uses the
   * default node of Tpetra
   */
  [[nodiscard]] std::unique_ptr<Operator> makeAmesos2Solver(
      const std::string& solver, const Options& options, const bool gpu);
  /*!
   * \return if a direct solver of Amesos2 is available
   * \param[in] solver: name of the solver
   */
  [[nodiscard]] bool isAmesos2SolverAvailable(const std::string& solver);
  /*!
   * \brief settings of mfem-mgis given to a Krylov solver of Belos: they
   * are added to its parameters when the solver accepts them and when the
   * parameters do not give them. The convergence settings of the hybrid
   * GMRES are given to its outer solver.
   */
  struct BelosSettings {
    //! \brief relative tolerance (`Convergence Tolerance`), if given
    std::optional<double> relative_tolerance;
    //! \brief maximum number of iterations (`Maximum Iterations`), if given
    std::optional<int> maximum_number_of_iterations;
    /*!
     * \brief verbosity level of mfem-mgis, if given: 0 (silent), 1
     * (warnings and final summary) or more (details of each iteration)
     */
    std::optional<int> verbosity_level;
  };

  /*!
   * \return a Krylov solver of Belos
   * \param[in] solver: name of the solver (`CG`, `GMRES`, `BICGSTAB`,
   * `GCRODR`, `MINRES`, `TFQMR`...)
   * \param[in] options: parameters of the solver
   * \param[in] settings: settings of mfem-mgis. The GMRES solvers keep 50
   * Krylov vectors and are restarted as many times as the maximum number of
   * iterations allows, unless their parameters give them.
   * \param[in] p: preconditioner of Trilinos, if any, applied to the
   * vectors of Tpetra without copy. It is not owned by the solver: its
   * matrix is set by the caller before the one of the solver, which then
   * shares it.
   * \param[in] e: preconditioner of mfem-mgis, if any, used if `p` is null
   * \param[in] gpu: if true, computes on the GPU if Trilinos is built for
   * one, on the host otherwise, unless a preconditioner of Trilinos is
   * given: the solver then computes at the same place as the preconditioner
   */
  [[nodiscard]] std::unique_ptr<Operator> makeBelosSolver(
      const std::string& solver,
      const Options& options,
      const BelosSettings& settings,
      const Operator* p,
      ExternalPreconditioner e,
      const bool gpu);
  /*!
   * \return true if Trilinos is built for a GPU (default node of Tpetra
   * whose execution space can't access the memory of the host). Tpetra is
   * then also used on the host when it is built for a node of the host
   * (serial or OpenMP).
   */
  [[nodiscard]] bool hasGPUNode() noexcept;
  /*!
   * \return true if Tpetra computes on the host when asked: Trilinos is
   * built without GPU, or Tpetra is also instantiated for its serial node
   */
  [[nodiscard]] bool hasHostNode() noexcept;

}  // end of namespace mfem_mgis::trilinos

#endif /* LIB_MFEM_MGIS_SRC_TRILINOSBACKEND_HXX */
