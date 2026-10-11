/*!
 * \file   src/AmgXSolvers.cxx
 * \brief  This file implements the solvers of the AmgX library of NVIDIA
 * built by the factory of linear solvers
 * \date   10/10/2026
 */

#include "AmgXSolvers.hxx"

#ifdef MFEM_MGIS_HAVE_AMGX

#include <array>
#include <cmath>
#include <limits>
#include <locale>
#include <vector>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <exception>
#include <type_traits>
#include <cuda_runtime_api.h>
#include "MGIS/Raise.hxx"
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "AmgXKernels.hxx"

namespace mfem_mgis::internals {

  namespace {

    /*!
     * \brief raises an exception if a call to AmgX failed
     * \param[in] rc: return code of the call
     * \param[in] what: description of the call
     */
    void checkAmgX(const AMGX_RC rc, const char* const what) {
      if (rc == AMGX_RC_OK) {
        return;
      }
      auto msg = std::array<char, 4096>{};
      AMGX_get_error_string(rc, msg.data(), static_cast<int>(msg.size()));
      mgis::raise(std::string{"AmgX: "} + what + " failed (" + msg.data() +
                  ")");
    }  // end of checkAmgX

    /*!
     * \brief prints the messages of AmgX on the first process
     * \param[in] msg: message
     * \param[in] length: length of the message
     */
    void printAmgXMessage(const char* msg, int length) {
      auto rank = 0;
#ifdef MFEM_USE_MPI
      auto initialized = 0;
      auto finalized = 0;
      MPI_Initialized(&initialized);
      MPI_Finalized(&finalized);
      if ((initialized != 0) && (finalized == 0)) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      }
#endif /* MFEM_USE_MPI */
      if (rank == 0) {
        getOutputStream().write(msg, length);
      }
    }  // end of printAmgXMessage

    /*!
     * \brief objects of AmgX shared by all the solvers of the process. The
     * resources of AmgX hold global objects (pool of memory of the device,
     * handles of cuSPARSE and cuBLAS) which are released when any resources
     * are destroyed, and keep the address of their communicator: the
     * resources are created once, for the first solver, on `MPI_COMM_WORLD`
     * (whose handle is kept here; AmgX duplicates it for the exchanges of
     * the distributed matrices), and destroyed with the last solver. Their
     * creation and destruction are not collective.
     */
    struct AmgXProcessResources {
      //! \brief number of solvers using the resources
      int users = 0;
      //! \brief if true, AmgX is initialised
      bool initialised = false;
      //! \brief configuration of the resources
      AMGX_config_handle configuration = nullptr;
      //! \brief resources
      AMGX_resources_handle resources = nullptr;
#ifdef MFEM_USE_MPI
      //! \brief communicator of the resources, whose address AmgX keeps
      MPI_Comm comm = MPI_COMM_NULL;
#endif /* MFEM_USE_MPI */
    };

    //! \return the objects of AmgX shared by the solvers of the process
    AmgXProcessResources& getAmgXProcessResources() noexcept {
      static auto r = AmgXProcessResources{};
      return r;
    }  // end of getAmgXProcessResources

    //! \brief destroys the objects of AmgX shared by the solvers
    void destroyAmgXProcessResources() noexcept {
      auto& r = getAmgXProcessResources();
      if (r.resources != nullptr) {
        static_cast<void>(AMGX_resources_destroy(r.resources));
        r.resources = nullptr;
      }
      if (r.configuration != nullptr) {
        static_cast<void>(AMGX_config_destroy(r.configuration));
        r.configuration = nullptr;
      }
#ifdef MFEM_USE_MPI
      r.comm = MPI_COMM_NULL;
#endif /* MFEM_USE_MPI */
      if (r.initialised) {
        static_cast<void>(AMGX_finalize());
        r.initialised = false;
      }
    }  // end of destroyAmgXProcessResources

    /*!
     * \return the resources of AmgX shared by the solvers of the process,
     * created for the first solver
     */
    AMGX_resources_handle acquireAmgXProcessResources() {
      auto& r = getAmgXProcessResources();
      if (r.users == 0) {
        try {
          checkAmgX(AMGX_register_print_callback(&printAmgXMessage),
                    "registering the print function");
          checkAmgX(AMGX_initialize(), "initialisation");
          r.initialised = true;
          // without the parameter `exception_handling`, the errors are
          // reported by the return codes (with it, AmgX aborts)
          checkAmgX(AMGX_config_create(&(r.configuration), "config_version=2"),
                    "creation of the configuration of the resources");
          // the current device of the process (see `mfem_mgis::initialize`)
          auto device = 0;
          if (cudaGetDevice(&device) != cudaSuccess) {
            mgis::raise("AmgX: no current GPU");
          }
#ifdef MFEM_USE_MPI
          r.comm = MPI_COMM_WORLD;
          checkAmgX(AMGX_resources_create(&(r.resources), r.configuration,
                                          &(r.comm), 1, &device),
                    "creation of the resources");
#else  /* MFEM_USE_MPI */
          checkAmgX(AMGX_resources_create(&(r.resources), r.configuration,
                                          nullptr, 1, &device),
                    "creation of the resources");
#endif /* MFEM_USE_MPI */
        } catch (...) {
          destroyAmgXProcessResources();
          throw;
        }
      }
      ++(r.users);
      return r.resources;
    }  // end of acquireAmgXProcessResources

    //! \brief releases the resources of AmgX for a destroyed solver
    void releaseAmgXProcessResources() noexcept {
      auto& r = getAmgXProcessResources();
      --(r.users);
      if (r.users == 0) {
        destroyAmgXProcessResources();
      }
    }  // end of releaseAmgXProcessResources

    //! \return the mode of AmgX associated with a precision
    AMGX_Mode getMode(const AmgXPrecision p) noexcept {
      return p == AmgXPrecision::SINGLE ? AMGX_mode_dFFI : AMGX_mode_dDDI;
    }  // end of getMode

    /*!
     * \return the values of a matrix in the precision of the matrices of
     * AmgX: the given values in double precision, or their copy in single
     * precision
     * \param[out] copy: values in single precision, if needed
     * \param[in] mode: mode of AmgX
     * \param[in] values: values in double precision
     * \param[in] n: number of values
     */
    const void* getMatrixValues(std::vector<float>& copy,
                                const AMGX_Mode mode,
                                const real* const values,
                                const size_type n) {
      if (mode == AMGX_mode_dFFI) {
        copy.assign(values, values + n);
        return copy.data();
      }
      return values;
    }  // end of getMatrixValues

    /*!
     * \return a power of two giving to the given positive value a value
     * between one and two, if the value is finite and if this power is
     * representable
     * \param[in] v: value, strictly positive
     */
    std::optional<real> getScalingFactor(const real v) noexcept {
      if (!std::isfinite(v)) {
        return {};
      }
      const auto f = std::ldexp(real{1}, -std::ilogb(v));
      if (!std::isfinite(f)) {
        return {};
      }
      return f;
    }  // end of getScalingFactor

  }  // end of namespace

#ifdef MFEM_USE_MPI

  AmgXSolver::AmgXSolver(MPI_Comm c,
                         const std::string& config,
                         const AmgXSolverOptions& o)
      : mode(getMode(o.precision)),
        iterative(o.iterative),
        scale_right_hand_side(o.scale_right_hand_side) {
    this->iterative_mode = false;
    auto size = 0;
    MPI_Comm_size(c, &size);
    if (size != 1) {
      // the resources of AmgX are created on `MPI_COMM_WORLD`
      int result = MPI_UNEQUAL;
      MPI_Comm_compare(c, MPI_COMM_WORLD, &result);
      if ((result != MPI_IDENT) && (result != MPI_CONGRUENT)) {
        mgis::raise(
            "AmgX: the parallel solvers of AmgX require a communicator "
            "containing all the processes, in the same order");
      }
    }
    MPI_Comm_dup(c, &(this->comm));
    try {
      this->create(config);
    } catch (...) {
      this->release();
      throw;
    }
  }  // end of AmgXSolver

#else /* MFEM_USE_MPI */

  AmgXSolver::AmgXSolver(const std::string& config, const AmgXSolverOptions& o)
      : mode(getMode(o.precision)),
        iterative(o.iterative),
        scale_right_hand_side(o.scale_right_hand_side) {
    this->iterative_mode = false;
    try {
      this->create(config);
    } catch (...) {
      this->release();
      throw;
    }
  }  // end of AmgXSolver

#endif /* MFEM_USE_MPI */

  void AmgXSolver::create(const std::string& config) {
    this->resources = acquireAmgXProcessResources();
    checkAmgX(AMGX_config_create(&(this->configuration), config.c_str()),
              "creation of the configuration");
    checkAmgX(AMGX_matrix_create(&(this->matrix), this->resources, this->mode),
              "creation of the matrix");
    checkAmgX(AMGX_vector_create(&(this->rhs), this->resources, this->mode),
              "creation of a vector");
    checkAmgX(
        AMGX_vector_create(&(this->solution), this->resources, this->mode),
        "creation of a vector");
    checkAmgX(AMGX_solver_create(&(this->solver), this->resources, this->mode,
                                 this->configuration),
              "creation of the solver");
  }  // end of create

  void AmgXSolver::release() noexcept {
    if (this->solver != nullptr) {
      static_cast<void>(AMGX_solver_destroy(this->solver));
      this->solver = nullptr;
    }
    if (this->matrix != nullptr) {
      static_cast<void>(AMGX_matrix_destroy(this->matrix));
      this->matrix = nullptr;
    }
    for (auto* const v : {&(this->rhs), &(this->solution)}) {
      if (*v != nullptr) {
        static_cast<void>(AMGX_vector_destroy(*v));
        *v = nullptr;
      }
    }
    if (this->configuration != nullptr) {
      static_cast<void>(AMGX_config_destroy(this->configuration));
      this->configuration = nullptr;
    }
    if (this->resources != nullptr) {
      releaseAmgXProcessResources();
      this->resources = nullptr;
    }
#ifdef MFEM_USE_MPI
    if (this->comm != MPI_COMM_NULL) {
      auto finalized = 0;
      MPI_Finalized(&finalized);
      if (finalized == 0) {
        MPI_Comm_free(&(this->comm));
      }
      this->comm = MPI_COMM_NULL;
    }
#endif /* MFEM_USE_MPI */
  }  // end of release

  AmgXSolver::~AmgXSolver() { this->release(); }

  void AmgXSolver::upload(const mfem::SparseMatrix& m) {
    auto values = std::vector<float>{};
    const auto nnz = m.NumNonZeroElems();
    checkAmgX(
        AMGX_matrix_upload_all(
            this->matrix, m.Height(), nnz, 1, 1, m.HostReadI(), m.HostReadJ(),
            getMatrixValues(values, this->mode, m.HostReadData(), nnz),
            nullptr),
        "upload of the matrix");
  }  // end of upload

#ifdef MFEM_USE_MPI

  void AmgXSolver::upload(const mfem::HypreParMatrix& m) {
    // limits of AmgX: global indices of rows of type int, rows on each
    // process (checked by all the processes, which fail together)
    if (m.GetGlobalNumRows() >
        static_cast<HYPRE_BigInt>(std::numeric_limits<int>::max())) {
      throw CollectiveFailure(
          "AmgX: too many rows (more than the range of int)");
    }
    auto min_rows = static_cast<int>(m.Height());
    MPI_Allreduce(MPI_IN_PLACE, &min_rows, 1, MPI_INT, MPI_MIN, this->comm);
    if (min_rows == 0) {
      throw CollectiveFailure("AmgX: a process owns no row of the matrix");
    }
    // local rows, with the global indices of the columns, merged on the
    // host
    m.HostRead();
    auto merged = std::unique_ptr<hypre_CSRMatrix, void (*)(hypre_CSRMatrix*)>(
        hypre_MergeDiagAndOffd(static_cast<hypre_ParCSRMatrix*>(
            const_cast<mfem::HypreParMatrix&>(m))),
        [](hypre_CSRMatrix* ptr) {
          static_cast<void>(hypre_CSRMatrixDestroy(ptr));
        });
    m.HypreRead();
    const auto n = hypre_CSRMatrixNumRows(merged.get());
    const auto nnz = hypre_CSRMatrixNumNonzeros(merged.get());
    const auto* const columns = hypre_CSRMatrixBigJ(merged.get());
    auto values = std::vector<float>{};
    const auto* const data = getMatrixValues(
        values, this->mode, hypre_CSRMatrixData(merged.get()), nnz);
    auto size = 0;
    MPI_Comm_size(this->comm, &size);
    if (size == 1) {
      // matrix of one process, whose communicator may differ from the one
      // of the resources: uploaded as a sequential matrix
      const auto local_columns = std::vector<int>(columns, columns + nnz);
      checkAmgX(AMGX_matrix_upload_all(this->matrix, n, nnz, 1, 1,
                                       hypre_CSRMatrixI(merged.get()),
                                       local_columns.data(), data, nullptr),
                "upload of the matrix");
      return;
    }
    // offsets of the rows of the processes, with the type of the global
    // indices of the columns
    auto offsets =
        std::vector<HYPRE_BigInt>(static_cast<std::size_t>(size) + 1);
    const auto first = m.GetRowStarts()[0];
    MPI_Allgather(&first, 1, HYPRE_MPI_BIG_INT, offsets.data(), 1,
                  HYPRE_MPI_BIG_INT, this->comm);
    offsets.back() = m.GetGlobalNumRows();
    auto distribution = AMGX_distribution_handle{};
    checkAmgX(AMGX_distribution_create(&distribution, this->configuration),
              "creation of the distribution");
    const auto d =
        std::unique_ptr<std::remove_pointer_t<AMGX_distribution_handle>,
                        void (*)(AMGX_distribution_handle)>(
            distribution, [](AMGX_distribution_handle ptr) {
              static_cast<void>(AMGX_distribution_destroy(ptr));
            });
    checkAmgX(AMGX_distribution_set_partition_data(
                  distribution, AMGX_DIST_PARTITION_OFFSETS, offsets.data()),
              "definition of the partition");
    checkAmgX(AMGX_distribution_set_32bit_colindices(
                  distribution, sizeof(HYPRE_BigInt) == 4 ? 1 : 0),
              "definition of the type of the indices");
    checkAmgX(AMGX_matrix_upload_distributed(
                  this->matrix, static_cast<int>(m.GetGlobalNumRows()), n, nnz,
                  1, 1, hypre_CSRMatrixI(merged.get()), columns, data, nullptr,
                  distribution),
              "upload of the matrix");
  }  // end of upload

#endif /* MFEM_USE_MPI */

  void AmgXSolver::SetOperator(const mfem::Operator& op) {
    // the solver is set up again only if the setup succeeds
    this->setup = false;
    this->status.reset();
    this->iterations.reset();
    if (const auto* const m = dynamic_cast<const mfem::SparseMatrix*>(&op);
        m != nullptr) {
      this->upload(*m);
#ifdef MFEM_USE_MPI
    } else if (const auto* const pm =
                   dynamic_cast<const mfem::HypreParMatrix*>(&op);
               pm != nullptr) {
      this->upload(*pm);
#endif /* MFEM_USE_MPI */
    } else {
      mgis::raise(
          "AmgXSolver::SetOperator: unsupported operator (a sparse matrix "
          "or a parallel matrix of hypre is expected)");
    }
    checkAmgX(AMGX_solver_setup(this->solver, this->matrix),
              "setup of the solver");
    checkAmgX(AMGX_vector_bind(this->rhs, this->matrix), "binding of a vector");
    checkAmgX(AMGX_vector_bind(this->solution, this->matrix),
              "binding of a vector");
    this->height = op.Height();
    this->width = op.Width();
    this->setup = true;
  }  // end of SetOperator

  void AmgXSolver::Mult(const mfem::Vector& b, mfem::Vector& x) const {
    // the resolution fails until AmgX gives its status
    this->status = AMGX_SOLVE_FAILED;
    this->iterations.reset();
    if (!this->setup) {
      mgis::raise("AmgXSolver::Mult: no operator set up");
    }
    if ((b.Size() != this->Height()) || (x.Size() != this->Width())) {
      mgis::raise("AmgXSolver::Mult: invalid sizes of the vectors");
    }
    const auto n = b.Size();
    const auto single = this->mode == AMGX_mode_dFFI;
    x.UseDevice(true);
    // factor of the right hand side, exact, given by the maximum of the
    // absolute values of its components (its norm may underflow or
    // overflow)
    auto factor = std::optional<real>{};
    if (this->scale_right_hand_side) {
      auto m = getMaximumAbsoluteValue(b, this->reduction_workspace);
#ifdef MFEM_USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &m, 1, MPI_DOUBLE, MPI_MAX, this->comm);
#endif /* MFEM_USE_MPI */
      if (std::fpclassify(m) == FP_ZERO) {
        // zero solution, as given by AmgX
        x = real{};
        this->status = AMGX_SOLVE_SUCCESS;
        if (this->iterative) {
          this->iterations = 0;
        }
        return;
      }
      factor = getScalingFactor(m);
    }
    if (single) {
      convertToSinglePrecision(this->single_rhs, b, factor.value_or(1));
      checkAmgX(AMGX_vector_upload(this->rhs, n, 1, this->single_rhs.Read()),
                "upload of the right hand side");
    } else if (factor.has_value()) {
      scale(this->scaled_rhs, b, *factor);
      checkAmgX(AMGX_vector_upload(this->rhs, n, 1, this->scaled_rhs.Read()),
                "upload of the right hand side");
    } else {
      checkAmgX(AMGX_vector_upload(this->rhs, n, 1, b.Read()),
                "upload of the right hand side");
    }
    // zero initial guess, uploaded: `AMGX_vector_set_zero` resizes the
    // vector without the halo of a distributed matrix
    if (this->zeros.Size() != n) {
      this->zeros.SetSize(n);
      this->zeros.UseDevice(true);
      this->zeros = real{};
      if (single) {
        convertToSinglePrecision(this->single_zeros, this->zeros);
      }
    }
    checkAmgX(AMGX_vector_upload(
                  this->solution, n, 1,
                  single ? static_cast<const void*>(this->single_zeros.Read())
                         : static_cast<const void*>(this->zeros.Read())),
              "initialisation of the solution");
    checkAmgX(AMGX_solver_solve_with_0_initial_guess(this->solver, this->rhs,
                                                     this->solution),
              "resolution");
    auto s = AMGX_SOLVE_STATUS{};
    checkAmgX(AMGX_solver_get_status(this->solver, &s), "access to the status");
    if (this->iterative) {
      auto it = 0;
      checkAmgX(AMGX_solver_get_iterations_number(this->solver, &it),
                "access to the number of iterations");
      this->iterations = it;
    }
    if (single) {
      this->single_solution.SetSize(n);
      checkAmgX(
          AMGX_vector_download(this->solution, this->single_solution.Write()),
          "download of the solution");
      convertToDoublePrecision(x, this->single_solution,
                               1 / factor.value_or(1));
    } else {
      checkAmgX(AMGX_vector_download(this->solution, x.Write()),
                "download of the solution");
      if (factor.has_value()) {
        x *= 1 / *factor;
      }
    }
    this->status = s;
  }  // end of Mult

  std::optional<int> AmgXSolver::getNumberOfIterations() const noexcept {
    return this->iterations;
  }  // end of getNumberOfIterations

  bool AmgXSolver::hasConverged() const noexcept {
    if ((!this->iterative) || (!this->status.has_value())) {
      return true;
    }
    return *(this->status) == AMGX_SOLVE_SUCCESS;
  }  // end of hasConverged

  //! \brief name of the option giving parameters of AmgX
  static const char* const AmgXParameters = "AmgXParameters";
  //! \brief name of the option giving the precision
  static const char* const Precision = "Precision";

  /*!
   * \return true if a GPU is visible, an error otherwise
   * \param[in, out] ctx: execution context
   */
  static bool checkGPU(Context& ctx) noexcept {
    auto n = 0;
    if ((cudaGetDeviceCount(&n) != cudaSuccess) || (n == 0)) {
      static_cast<void>(cudaGetLastError());
      return ctx.registerErrorMessage("AmgX: no GPU available");
    }
    return true;
  }  // end of checkGPU

  /*!
   * \return the precision given by the options
   * \param[in, out] ctx: execution context
   * \param[in] params: options
   */
  static std::optional<AmgXPrecision> getPrecision(
      Context& ctx, const Parameters& params) noexcept {
    const auto oprecision =
        get_if<std::string>(ctx, params, Precision, std::string{"Double"});
    if (isInvalid(oprecision)) {
      return {};
    }
    if (*oprecision == "Double") {
      return AmgXPrecision::DOUBLE;
    } else if (*oprecision == "Single") {
      return AmgXPrecision::SINGLE;
    }
    return ctx.registerErrorMessage("AmgX: invalid precision '" + *oprecision +
                                    "' (expected 'Double' or 'Single')");
  }  // end of getPrecision

  /*!
   * \return a solver, or an error
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] c: configuration
   * \param[in] o: options of the solver
   */
  template <bool parallel>
  static std::unique_ptr<AmgXSolver> makeSolver(
      Context& ctx,
      FiniteElementSpace<parallel>& fespace,
      const std::string& c,
      const AmgXSolverOptions& o) noexcept {
    try {
#ifdef MFEM_USE_MPI
      auto comm = MPI_Comm{MPI_COMM_SELF};
      if constexpr (parallel) {
        comm = fespace.GetComm();
      } else {
        static_cast<void>(fespace);
      }
      return std::make_unique<AmgXSolver>(comm, c, o);
#else  /* MFEM_USE_MPI */
      static_cast<void>(fespace);
      return std::make_unique<AmgXSolver>(c, o);
#endif /* MFEM_USE_MPI */
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(e.what());
    } catch (...) {
      return ctx.registerErrorMessage("AmgX: unknown error");
    }
  }  // end of makeSolver

  /*!
   * \return a stream writing a configuration of AmgX in its format
   * `name=value`, the reals being written exactly
   */
  static std::ostringstream makeConfigurationStream() {
    auto c = std::ostringstream{};
    c.imbue(std::locale::classic());
    c << std::setprecision(std::numeric_limits<real>::max_digits10);
    return c;
  }  // end of makeConfigurationStream

  /*!
   * \brief writes the parameters of the default multigrid of MFEM (scope
   * `amg`): one V-cycle of the classical algebraic multigrid with a
   * L1-Jacobi smoother
   * \param[in, out] c: configuration
   * \param[in] verbose: print the statistics of the grids
   */
  static void writeDefaultMultigridParameters(std::ostringstream& c,
                                              const bool verbose) {
    c << "amg:max_iters=1, amg:cycle=V, amg:max_levels=100, "
      << "amg:smoother(smoother)=JACOBI_L1, amg:presweeps=1, "
      << "amg:postsweeps=1, amg:interpolator=D2, amg:max_row_sum=0.9, "
      << "amg:strength_threshold=0.25";
    if (verbose) {
      c << ", amg:print_grid_stats=1";
    }
  }  // end of writeDefaultMultigridParameters

  /*!
   * \return if the given string can be written in a configuration of AmgX
   * as a name or a value (no separator of AmgX, no space)
   * \param[in] s: string
   */
  static bool isValidConfigurationString(const std::string& s) noexcept {
    return (!s.empty()) && (s.find_first_of(",;= \t\n\r") == std::string::npos);
  }  // end of isValidConfigurationString

  /*!
   * \brief writes parameters of AmgX given by the user (option
   * `AmgXParameters`), which override the previous values: a parameter
   * without scope is given to the multigrid (scope `amg`), a parameter
   * `scope:name` as is
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] c: configuration
   * \param[in] params: options
   * \param[in] forbidden_scope: scope of the parameters refused, if not
   * empty
   */
  static bool writeAmgXParameters(
      Context& ctx,
      std::ostringstream& c,
      const Parameters& params,
      const std::string& forbidden_scope = "") noexcept {
    if (!contains(params, AmgXParameters)) {
      return true;
    }
    const auto oparameters = get<Parameters>(ctx, params, AmgXParameters);
    if (isInvalid(oparameters)) {
      return false;
    }
    for (const auto& [name, value] : *oparameters) {
      if (!isValidConfigurationString(name)) {
        return ctx.registerErrorMessage(
            "AmgX: invalid name of a parameter of AmgX '" + name + "'");
      }
      if ((!forbidden_scope.empty()) &&
          (name.starts_with(forbidden_scope + ":"))) {
        return ctx.registerErrorMessage(
            "AmgX: the parameter of AmgX '" + name +
            "' is given by the options of the solver");
      }
      c << ", " << (name.find(':') == std::string::npos ? "amg:" : "") << name
        << '=';
      if (is<bool>(value)) {
        c << (*get<bool>(ctx, value) ? 1 : 0);
      } else if (is<int>(value)) {
        c << *get<int>(ctx, value);
      } else if (is<real>(value)) {
        c << *get<real>(ctx, value);
      } else if (is<std::string>(value)) {
        const auto& v = *get<std::string>(ctx, value);
        if (!isValidConfigurationString(v)) {
          return ctx.registerErrorMessage(
              "AmgX: invalid value of the parameter of AmgX '" + name + "'");
        }
        c << v;
      } else {
        return ctx.registerErrorMessage(
            "AmgX: invalid type of the parameter of AmgX '" + name +
            "' (boolean, integer, real or string expected)");
      }
    }
    return true;
  }  // end of writeAmgXParameters

  template <bool parallel>
  std::unique_ptr<AmgXSolver> makeAmgXPreconditioner(
      Context& ctx,
      FiniteElementSpace<parallel>& fespace,
      const Parameters& params) noexcept {
    const char* const Configuration = "Configuration";
    const char* const ConfigurationFile = "ConfigurationFile";
    const char* const VerbosityLevel = "VerbosityLevel";
    if (!checkParameters(ctx, params,
                         std::vector<std::string>{
                             Configuration, ConfigurationFile, VerbosityLevel,
                             AmgXParameters, Precision})) {
      return {};
    }
    const auto oconfiguration =
        get_if<std::string>(ctx, params, Configuration, std::string{});
    const auto ofile =
        get_if<std::string>(ctx, params, ConfigurationFile, std::string{});
    const auto overbosity = get_if<int>(ctx, params, VerbosityLevel, 0);
    const auto oprecision = getPrecision(ctx, params);
    if (isInvalid(oconfiguration) || isInvalid(ofile) ||
        isInvalid(overbosity) || isInvalid(oprecision)) {
      return {};
    }
    const auto has_configuration = contains(params, Configuration);
    const auto has_file = contains(params, ConfigurationFile);
    if (has_configuration && has_file) {
      return ctx.registerErrorMessage(
          "AmgX: the options 'Configuration' and 'ConfigurationFile' are "
          "exclusive");
    }
    if ((has_configuration || has_file) && (contains(params, VerbosityLevel) ||
                                            contains(params, AmgXParameters))) {
      return ctx.registerErrorMessage(
          "AmgX: the verbosity level and the parameters of AmgX only apply "
          "to the default configuration (a configuration given defines all "
          "the parameters of AmgX)");
    }
    if (has_configuration && oconfiguration->empty()) {
      return ctx.registerErrorMessage("AmgX: empty configuration");
    }
    if (!checkGPU(ctx)) {
      return {};
    }
    if (has_file) {
      auto f = std::ifstream{*ofile};
      if (!f) {
        return ctx.registerErrorMessage(
            "AmgX: can't open the configuration file '" + *ofile + "'");
      }
      auto c = std::ostringstream{};
      c << f.rdbuf();
      return makeSolver<parallel>(ctx, fespace, c.str(),
                                  AmgXSolverOptions{.precision = *oprecision});
    }
    if (has_configuration) {
      return makeSolver<parallel>(ctx, fespace, *oconfiguration,
                                  AmgXSolverOptions{.precision = *oprecision});
    }
    // default configuration, completed by the parameters given
    auto c = makeConfigurationStream();
    c << "config_version=2, solver(amg)=AMG, ";
    writeDefaultMultigridParameters(c, *overbosity > 0);
    if (*overbosity > 0) {
      // the residual is computed to print the statistics: a negative
      // tolerance prevents a convergence before the cycles (a small right
      // hand side would give a zero solution)
      c << ", amg:monitor_residual=1, amg:tolerance=-1, "
        << "amg:print_solve_stats=1, amg:obtain_timings=1";
    }
    if (!writeAmgXParameters(ctx, c, params)) {
      return {};
    }
    // the cycles of the default configuration do not depend on the scale of
    // the right hand side: in single precision, it is scaled to stay in the
    // range of the single precision
    return makeSolver<parallel>(
        ctx, fespace, c.str(),
        AmgXSolverOptions{
            .precision = *oprecision,
            .scale_right_hand_side = *oprecision == AmgXPrecision::SINGLE});
  }  // end of makeAmgXPreconditioner

  template <bool parallel>
  std::unique_ptr<AmgXSolver> makeAmgXConjugateGradientSolver(
      Context& ctx,
      FiniteElementSpace<parallel>& fespace,
      const Parameters& params) noexcept {
    using Problem = AbstractNonLinearEvolutionProblem;
    auto allowed_parameters = getIterativeSolverParametersList();
    allowed_parameters.push_back(AmgXParameters);
    allowed_parameters.push_back(Precision);
    if (!checkParameters(ctx, params, allowed_parameters)) {
      return {};
    }
    // default values of the solver configuration of MFEM
    const auto overbosity =
        get_if<int>(ctx, params, Problem::SolverVerbosityLevel, 0);
    const auto ortol =
        get_if<real>(ctx, params, Problem::SolverRelativeTolerance, 1e-12);
    const auto oatol =
        get_if<real>(ctx, params, Problem::SolverAbsoluteTolerance, 0);
    const auto on =
        get_if<int>(ctx, params, Problem::SolverMaximumNumberOfIterations, 150);
    const auto oprecision = getPrecision(ctx, params);
    if (isInvalid(overbosity) || isInvalid(ortol) || isInvalid(oatol) ||
        isInvalid(on) || isInvalid(oprecision)) {
      return {};
    }
    if ((*ortol < 0) || (*oatol < 0) || ((*ortol <= 0) && (*oatol <= 0))) {
      return ctx.registerErrorMessage(
          "AmgXPCG: the tolerances must be positive, one of them strictly");
    }
    if (*on < 1) {
      return ctx.registerErrorMessage(
          "AmgXPCG: the maximum number of iterations must be strictly "
          "positive");
    }
    // AmgX reports a convergence when the norm of the residual is lower
    // than the maximum of 1e-20 and of 1e-12 times the norm of the initial
    // residual in double precision, 1e-6 in single precision: lower
    // relative tolerances would be ignored. Without absolute tolerance, the
    // right hand side is scaled (its norm is at least one): the floor of
    // 1e-20 is negligible. With an absolute tolerance (double precision
    // only), the absolute tolerance must be greater than this floor; an
    // absolute tolerance alone is reached or the norm of the residual is
    // reduced by 1e-12.
    const auto floor =
        (*oprecision == AmgXPrecision::SINGLE) ? real{1e-6} : real{1e-12};
    if ((*ortol > 0) && (*ortol < floor)) {
      return ctx.registerErrorMessage(
          "AmgXPCG: the relative tolerance must be at least " +
          std::string{*oprecision == AmgXPrecision::SINGLE ? "1e-6" : "1e-12"} +
          " (floor of the convergence criteria of AmgX)");
    }
    if ((*oprecision == AmgXPrecision::SINGLE) && (*ortol <= 0)) {
      return ctx.registerErrorMessage(
          "AmgXPCG: in single precision, a relative tolerance is required");
    }
    if ((*oprecision == AmgXPrecision::SINGLE) && (*oatol > 0)) {
      return ctx.registerErrorMessage(
          "AmgXPCG: in single precision, the right hand side is scaled to "
          "stay in the range of the single precision: an absolute tolerance "
          "is not supported");
    }
    if ((*oatol > 0) && (*oatol < 1e-20)) {
      return ctx.registerErrorMessage(
          "AmgXPCG: the absolute tolerance must be at least 1e-20 (floor of "
          "the convergence criteria of AmgX)");
    }
    // conjugate gradient preconditioned by the default multigrid, completed
    // by the parameters given
    auto c = makeConfigurationStream();
    c << "config_version=2, solver(main)=PCG, main:max_iters=" << *on
      << ", main:monitor_residual=1, main:norm=L2, ";
    if (*oatol > 0) {
      // converged if the norm of the residual is lower than `tolerance` or
      // than `alt_rel_tolerance` times the norm of the initial residual
      c << "main:convergence=COMBINED_REL_INI_ABS, main:tolerance=" << *oatol
        << ", main:alt_rel_tolerance=" << *ortol << ", ";
    } else {
      c << "main:convergence=RELATIVE_INI_CORE, main:tolerance=" << *ortol
        << ", ";
    }
    if (*overbosity > 0) {
      c << "main:obtain_timings=1, main:print_solve_stats=1, ";
    }
    c << "main:preconditioner(amg)=AMG, ";
    writeDefaultMultigridParameters(c, *overbosity > 0);
    // the parameters of the conjugate gradient (scope `main`) are given by
    // the options of the solver
    if (!writeAmgXParameters(ctx, c, params, "main")) {
      return {};
    }
    if (!checkGPU(ctx)) {
      return {};
    }
    // the right hand side is scaled for a relative criterion only
    return makeSolver<parallel>(
        ctx, fespace, c.str(),
        AmgXSolverOptions{.precision = *oprecision,
                          .iterative = true,
                          .scale_right_hand_side = !(*oatol > 0)});
  }  // end of makeAmgXConjugateGradientSolver

  template std::unique_ptr<AmgXSolver> makeAmgXPreconditioner<false>(
      Context&, FiniteElementSpace<false>&, const Parameters&) noexcept;
  template std::unique_ptr<AmgXSolver> makeAmgXConjugateGradientSolver<false>(
      Context&, FiniteElementSpace<false>&, const Parameters&) noexcept;
#ifdef MFEM_USE_MPI
  template std::unique_ptr<AmgXSolver> makeAmgXPreconditioner<true>(
      Context&, FiniteElementSpace<true>&, const Parameters&) noexcept;
  template std::unique_ptr<AmgXSolver> makeAmgXConjugateGradientSolver<true>(
      Context&, FiniteElementSpace<true>&, const Parameters&) noexcept;
#endif /* MFEM_USE_MPI */

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_MGIS_HAVE_AMGX */
