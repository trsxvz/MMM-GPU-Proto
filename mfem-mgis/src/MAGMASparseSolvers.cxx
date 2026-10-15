/*!
 * \file   src/MAGMASparseSolvers.cxx
 * \brief  This file implements the module of the sparse part of MAGMA
 * (MAGMA sparse), which computes on one GPU: Krylov solvers in sequential,
 * and preconditioners (Jacobi, incomplete factorizations of cuSPARSE or
 * hipSPARSE) applied to the matrix in sequential and to the diagonal block
 * of each process in parallel (see `LinearSolverModules.hxx`)
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#ifdef MFEM_MGIS_HAVE_MAGMA_SPARSE

#include <cmath>
#include <array>
#include <tuple>
#include <memory>
#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <algorithm>
#include <exception>
#include <type_traits>
#include "magma_v2.h"
#include "magmasparse.h"
#include "mfem/general/backends.hpp"
#include "mfem/general/device.hpp"
#include "mfem/linalg/sparsemat.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverConvergenceReport.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis::internals {

  static_assert(std::is_same_v<real, double>,
                "the module of MAGMA is only built in double precision");

  /*!
   * \brief raise an exception if a call to MAGMA failed
   * \param[in] info: return code of MAGMA
   * \param[in] what: description of the call
   */
  static void checkMAGMA(const magma_int_t info, const char* const what) {
    if (info != MAGMA_SUCCESS) {
      mgis::raise(std::string{"MAGMA: "} + what + " failed (" +
                  magma_strerror(info) + ")");
    }
  }  // end of checkMAGMA

  /*!
   * \brief MAGMA, initialized for each queue (MAGMA counts its
   * initializations), and a queue of MAGMA on the current device. The queue
   * has its own stream: the device is synchronized before the vectors
   * written by MFEM are given to MAGMA, and the queue before they are given
   * back to MFEM.
   */
  struct MAGMAQueue {
    MAGMAQueue() {
      const auto info = magma_init();
      if (info != MAGMA_SUCCESS) {
        // MAGMA counts the failed initializations too
        magma_finalize();
        checkMAGMA(info, "initialization");
      }
      // MAGMA is initialized without GPU: the failures would only appear
      // when using it
      auto device = magma_device_t{};
      auto number_of_devices = magma_int_t{};
      magma_getdevices(&device, 1, &number_of_devices);
      if (number_of_devices == 0) {
        magma_finalize();
        mgis::raise("MAGMA: no GPU available");
      }
      magma_getdevice(&device);
      magma_queue_create(device, &(this->queue));
    }
    MAGMAQueue(MAGMAQueue&&) = delete;
    MAGMAQueue(const MAGMAQueue&) = delete;
    MAGMAQueue& operator=(MAGMAQueue&&) = delete;
    MAGMAQueue& operator=(const MAGMAQueue&) = delete;
    //! \brief destructor
    ~MAGMAQueue() {
      magma_queue_destroy(this->queue);
      magma_finalize();
    }
    //! \brief queue
    magma_queue_t queue = nullptr;
  };  // end of struct MAGMAQueue

  /*!
   * \brief matrix in the compressed row format of MAGMA, on the host, the
   * columns of each row being sorted (as the incomplete factorizations of
   * cuSPARSE and hipSPARSE require)
   */
  struct MAGMAHostMatrix {
    //! \brief offsets of the rows
    std::vector<magma_index_t> rows;
    //! \brief columns
    std::vector<magma_index_t> columns;
    //! \brief values
    std::vector<double> values;
    //! \brief number of rows
    magma_int_t size = 0;
    /*!
     * \return the matrix of MAGMA referring to this one (not to be freed)
     * \param[in] queue: queue
     */
    [[nodiscard]] magma_d_matrix wrap(magma_queue_t queue) {
      auto m = magma_d_matrix{};
      checkMAGMA(
          magma_dcsrset(this->size, this->size, this->rows.data(),
                        this->columns.data(), this->values.data(), &m, queue),
          "creation of a matrix");
      return m;
    }
  };  // end of struct MAGMAHostMatrix

  /*!
   * \return the matrix of MAGMA of a sparse matrix of MFEM (read on the
   * host)
   * \param[in] a: matrix
   */
  [[nodiscard]] static MAGMAHostMatrix makeMAGMAHostMatrix(
      const mfem::SparseMatrix& a) {
    if (a.Height() != a.Width()) {
      mgis::raise("MAGMA: the matrix is not square");
    }
    const auto n = a.Height();
    const auto* const I = a.HostReadI();
    const auto* const J = a.HostReadJ();
    const auto* const A = a.HostReadData();
    auto m = MAGMAHostMatrix{};
    m.size = static_cast<magma_int_t>(n);
    m.rows.reserve(n + 1);
    m.columns.reserve(I[n]);
    m.values.reserve(I[n]);
    m.rows.push_back(0);
    auto row = std::vector<std::pair<magma_index_t, double>>{};
    for (int i = 0; i != n; ++i) {
      row.clear();
      for (int k = I[i]; k != I[i + 1]; ++k) {
        row.emplace_back(static_cast<magma_index_t>(J[k]), A[k]);
      }
      std::sort(row.begin(), row.end(),
                [](const auto& u, const auto& v) { return u.first < v.first; });
      for (const auto& [c, v] : row) {
        m.columns.push_back(c);
        m.values.push_back(v);
      }
      m.rows.push_back(static_cast<magma_index_t>(m.columns.size()));
    }
    return m;
  }  // end of makeMAGMAHostMatrix

  /*!
   * \brief check that the diagonal entries of a matrix are not null, as the
   * preconditioners of MAGMA require: MAGMA only detects the missing
   * diagonal entries of the Jacobi preconditioner, and none of the
   * incomplete factorizations
   * \param[in] m: matrix
   */
  static void checkDiagonal(const MAGMAHostMatrix& m) {
    for (magma_int_t i = 0; i != m.size; ++i) {
      const auto b = m.columns.begin() + m.rows[i];
      const auto e = m.columns.begin() + m.rows[i + 1];
      const auto p = std::lower_bound(b, e, static_cast<magma_index_t>(i));
      const auto v =
          (p != e) && (*p == i)
              ? m.values[static_cast<std::size_t>(p - m.columns.begin())]
              : 0.;
      if (!(std::abs(v) > 0)) {
        mgis::raise("MAGMA: null or missing diagonal entry in the row " +
                    std::to_string(i) +
                    " of the matrix given to a preconditioner of MAGMA");
      }
    }
  }  // end of checkDiagonal

  /*!
   * \brief vector of MAGMA allocated on the device
   */
  struct MAGMADeviceVector {
    //! \param[in] q: queue
    explicit MAGMADeviceVector(magma_queue_t q) : queue(q) {}
    MAGMADeviceVector(MAGMADeviceVector&&) = delete;
    MAGMADeviceVector(const MAGMADeviceVector&) = delete;
    MAGMADeviceVector& operator=(MAGMADeviceVector&&) = delete;
    MAGMADeviceVector& operator=(const MAGMADeviceVector&) = delete;
    /*!
     * \brief allocate the vector, if its size changes
     * \param[in] n: size
     */
    void resize(const magma_int_t n) {
      if (this->allocated && (this->vector.num_rows == n)) {
        return;
      }
      this->release();
      checkMAGMA(
          magma_dvinit(&(this->vector), Magma_DEV, n, 1, 0., this->queue),
          "allocation of a vector");
      this->allocated = true;
    }
    //! \brief release the vector
    void release() noexcept {
      if (this->allocated) {
        magma_dmfree(&(this->vector), this->queue);
        this->allocated = false;
      }
    }
    //! \brief destructor
    ~MAGMADeviceVector() { this->release(); }
    //! \brief vector
    magma_d_matrix vector = {};

   private:
    //! \brief queue
    magma_queue_t queue;
    //! \brief if true, the vector is allocated
    bool allocated = false;
  };  // end of struct MAGMADeviceVector

  /*!
   * \brief vectors of MAGMA on the device, kept between the resolutions,
   * holding the copies of the vectors of MFEM when MFEM does not compute on
   * the GPU
   */
  struct MAGMAWorkspace {
    //! \param[in] q: queue
    explicit MAGMAWorkspace(magma_queue_t q) : rhs(q), unknowns(q) {}
    //! \brief right hand side
    MAGMADeviceVector rhs;
    //! \brief unknowns
    MAGMADeviceVector unknowns;
  };  // end of struct MAGMAWorkspace

  //! \brief initialization of the unknowns given to MAGMA
  enum struct MAGMAUnknownsInitialization {
    //! \brief the unknowns of MFEM are the initial guess
    INITIAL_GUESS,
    //! \brief the initial guess is null
    NULL_INITIAL_GUESS,
    //! \brief the unknowns are only written by MAGMA
    NONE
  };

  /*!
   * \brief vectors given to MAGMA on the device: the memory of the vectors
   * of MFEM when MFEM computes on the GPU and the vectors use its memory
   * (`UseDevice`), copies in a workspace otherwise
   */
  struct MAGMADeviceVectors {
    /*!
     * \param[in] b: right hand side
     * \param[in, out] x: unknowns
     * \param[in] init: initialization of the unknowns
     * \param[in, out] w: workspace
     * \param[in] q: queue
     */
    MAGMADeviceVectors(const mfem::Vector& b,
                       mfem::Vector& x,
                       const MAGMAUnknownsInitialization init,
                       MAGMAWorkspace& w,
                       magma_queue_t q)
        : solution(x),
          queue(q),
          on_device(mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                                         mfem::Backend::HIP_MASK) &&
                    b.UseDevice() && x.UseDevice()) {
      using Initialization = MAGMAUnknownsInitialization;
      const auto n = static_cast<magma_int_t>(b.Size());
      if (this->on_device) {
        if (init == Initialization::NULL_INITIAL_GUESS) {
          x = 0.;
        }
        auto* const pb = const_cast<double*>(b.Read());
        auto* const px =
            init == Initialization::NONE ? x.Write() : x.ReadWrite();
        // the kernels of MFEM writing the vectors are completed (the queue
        // of MAGMA does not wait for them)
        MFEM_DEVICE_SYNC;
        checkMAGMA(magma_dvset_dev(n, 1, pb, &(this->rhs), this->queue),
                   "creation of a vector");
        checkMAGMA(magma_dvset_dev(n, 1, px, &(this->unknowns), this->queue),
                   "creation of a vector");
        return;
      }
      w.rhs.resize(n);
      w.unknowns.resize(n);
      this->rhs = w.rhs.vector;
      this->unknowns = w.unknowns.vector;
      magma_dsetvector(n, b.HostRead(), 1, this->rhs.dval, 1, this->queue);
      if (init == Initialization::INITIAL_GUESS) {
        magma_dsetvector(n, x.HostRead(), 1, this->unknowns.dval, 1,
                         this->queue);
      } else if (init == Initialization::NULL_INITIAL_GUESS) {
        this->setNullUnknowns();
      }
    }
    MAGMADeviceVectors(MAGMADeviceVectors&&) = delete;
    MAGMADeviceVectors(const MAGMADeviceVectors&) = delete;
    MAGMADeviceVectors& operator=(MAGMADeviceVectors&&) = delete;
    MAGMADeviceVectors& operator=(const MAGMADeviceVectors&) = delete;
    //! \brief set the unknowns to zero
    void setNullUnknowns() {
      magmablas_dlaset(MagmaFull, this->unknowns.num_rows, 1, 0., 0.,
                       this->unknowns.dval, this->unknowns.num_rows,
                       this->queue);
    }
    //! \brief give the unknowns back to MFEM
    void finish() {
      if (this->on_device) {
        magma_queue_sync(this->queue);
        return;
      }
      magma_dgetvector(this->unknowns.num_rows, this->unknowns.dval, 1,
                       this->solution.HostWrite(), 1, this->queue);
    }
    //! \brief right hand side
    magma_d_matrix rhs = {};
    //! \brief unknowns
    magma_d_matrix unknowns = {};

   private:
    //! \brief unknowns of MFEM
    mfem::Vector& solution;
    //! \brief queue
    magma_queue_t queue;
    //! \brief if true, the vectors of MFEM are given to MAGMA
    const bool on_device;
  };  // end of struct MAGMADeviceVectors

  /*!
   * \return true if a return code of a Krylov solver of MAGMA only reports
   * a non convergence
   * \param[in] info: return code
   */
  [[nodiscard]] static bool isNonConvergence(const magma_int_t info) noexcept {
    return (info == MAGMA_SLOW_CONVERGENCE) || (info == MAGMA_DIVERGENCE) ||
           (info == MAGMA_NONSPD) || (info == MAGMA_NOTCONVERGED);
  }  // end of isNonConvergence

  /*!
   * \brief preconditioner of MAGMA: Jacobi, incomplete LU or Cholesky
   * factorizations (of cuSPARSE or hipSPARSE, without fill-in), applied by
   * the Krylov solvers of MAGMA or, as an `mfem::Solver`, by any solver. In
   * parallel, the preconditioner is applied to the diagonal block of each
   * process (block Jacobi preconditioner whose blocks are the processes).
   */
  struct MAGMAPreconditioner final : public mfem::Solver {
    //! \param[in] t: type of the preconditioner
    explicit MAGMAPreconditioner(const magma_solver_type t)
        : shared_queue(std::make_shared<MAGMAQueue>()),
          workspace(shared_queue->queue),
          work(shared_queue->queue),
          type(t) {}
    /*!
     * \brief prepare the preconditioner for a Krylov solver of MAGMA
     * \param[in] s: preconditioned Krylov solver (`Magma_PQMR` requires the
     * transposed factors)
     */
    void setKrylovSolver(const magma_solver_type s) noexcept {
      this->krylov_solver = s;
    }
    //! \return the type of the preconditioner
    [[nodiscard]] magma_solver_type getType() const noexcept {
      return this->type;
    }
    //! \return the queue of MAGMA
    [[nodiscard]] std::shared_ptr<MAGMAQueue> getQueue() const noexcept {
      return this->shared_queue;
    }
    //! \return true if the preconditioner is set up
    [[nodiscard]] bool isReady() const noexcept { return this->ready; }
    //! \return the preconditioner of MAGMA, updated when it is applied
    [[nodiscard]] magma_d_preconditioner* getMAGMAPreconditioner() noexcept {
      return &(this->preconditioner);
    }
    /*!
     * \brief set up the preconditioner of a matrix
     * \param[in] m: matrix
     */
    void setup(MAGMAHostMatrix& m) {
      auto* const q = this->shared_queue->queue;
      this->release();
      if (m.size == 0) {
        // empty diagonal block of a process: nothing to apply
        this->height = this->width = 0;
        this->ready = true;
        return;
      }
      checkDiagonal(m);
      auto solver = magma_d_solver_par{};
      solver.solver = this->krylov_solver;
      this->preconditioner = magma_d_preconditioner{};
      checkMAGMA(magma_dsolverinfo_init(&solver, &(this->preconditioner), q),
                 "initialization of a preconditioner");
      this->preconditioner.solver = this->type;
      this->preconditioner.trisolver = Magma_CUSOLVE;
      this->preconditioner.maxiter = 1;
      this->allocated = true;
      auto a = m.wrap(q);
      auto unused = magma_d_matrix{};
      checkMAGMA(
          magma_d_precondsetup(a, unused, &solver, &(this->preconditioner), q),
          "setup of the preconditioner");
      this->work.resize(m.size);
      this->height = this->width = static_cast<int>(m.size);
      this->ready = true;
    }
    void SetOperator(const mfem::Operator& op) override {
      this->release();
      if (const auto* const a = dynamic_cast<const mfem::SparseMatrix*>(&op);
          a != nullptr) {
        auto m = makeMAGMAHostMatrix(*a);
        this->setup(m);
        return;
      }
#ifdef MFEM_USE_MPI
      if (const auto* const a = dynamic_cast<const mfem::HypreParMatrix*>(&op);
          a != nullptr) {
        // collective: the setup may fail on some processes only (a null
        // diagonal entry, for instance)
        auto diagonal = mfem::SparseMatrix{};
        a->GetDiag(diagonal);
        auto msg = std::string{};
        try {
          auto m = makeMAGMAHostMatrix(diagonal);
          this->setup(m);
        } catch (std::exception& e) {
          msg = e.what();
        }
        auto failed = int{msg.empty() ? 0 : 1};
        MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_MAX, a->GetComm());
        if (failed != 0) {
          this->release();
          throw CollectiveFailure(
              "MAGMAPreconditioner::SetOperator: the setup of the "
              "preconditioner of MAGMA failed" +
              (msg.empty() ? std::string{" on another process"}
                           : " (" + msg + ")"));
        }
        return;
      }
#endif /* MFEM_USE_MPI */
      mgis::raise(
          "MAGMAPreconditioner::SetOperator: the operator is not an "
          "assembled matrix");
    }
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override {
      if (!this->ready) {
        mgis::raise("MAGMAPreconditioner::Mult: the operator is not set");
      }
      if ((b.Size() != this->Width()) || (x.Size() != this->Height())) {
        mgis::raise("MAGMAPreconditioner::Mult: invalid sizes");
      }
      if (this->Height() == 0) {
        return;
      }
      auto* const q = this->shared_queue->queue;
      auto v = MAGMADeviceVectors(b, x, MAGMAUnknownsInitialization::NONE,
                                  this->workspace, q);
      auto unused = magma_d_matrix{};
      checkMAGMA(magma_d_applyprecond_left(MagmaNoTrans, unused, v.rhs,
                                           &(this->work.vector),
                                           &(this->preconditioner), q),
                 "application of the preconditioner");
      checkMAGMA(
          magma_d_applyprecond_right(MagmaNoTrans, unused, this->work.vector,
                                     &(v.unknowns), &(this->preconditioner), q),
          "application of the preconditioner");
      v.finish();
    }
    //! \brief destructor
    ~MAGMAPreconditioner() override { this->release(); }

   private:
    //! \brief release the memory of MAGMA
    void release() noexcept {
      auto* const q = this->shared_queue->queue;
      if (this->allocated) {
        magma_dprecondfree(&(this->preconditioner), q);
        this->allocated = false;
      }
      this->ready = false;
    }
    //! \brief queue of MAGMA, shared with the Krylov solver of MAGMA
    const std::shared_ptr<MAGMAQueue> shared_queue;
    //! \brief copies of the vectors of MFEM
    mutable MAGMAWorkspace workspace;
    //! \brief work vector
    mutable MAGMADeviceVector work;
    //! \brief type of the preconditioner
    const magma_solver_type type;
    //! \brief Krylov solver of MAGMA using the preconditioner, if any
    magma_solver_type krylov_solver = Magma_NONE;
    //! \brief preconditioner of MAGMA (its timers are updated when applied)
    mutable magma_d_preconditioner preconditioner = {};
    //! \brief if true, the preconditioner holds memory of MAGMA
    bool allocated = false;
    //! \brief if true, the preconditioner is set up
    bool ready = false;
  };  // end of struct MAGMAPreconditioner

  //! \brief settings of a Krylov solver of MAGMA
  struct MAGMAKrylovSolverSettings {
    //! \brief name of the solver
    std::string name;
    //! \brief Krylov solver
    magma_solver_type solver = Magma_CG;
    //! \brief relative tolerance (default of MAGMA)
    real relative_tolerance = 1e-10;
    //! \brief absolute tolerance (default of MAGMA)
    real absolute_tolerance = 1e-16;
    //! \brief maximum number of iterations (default of MAGMA)
    magma_int_t maximum_number_of_iterations = 1000;
    /*!
     * \brief dimension of the Krylov space of GMRES (50 by default), or of
     * the shadow space of IDR (1 by default, given by any value but 50)
     */
    magma_int_t restart = 50;
    //! \brief verbosity level
    int verbosity = 0;
  };  // end of struct MAGMAKrylovSolverSettings

  /*!
   * \brief Krylov solver of MAGMA, in sequential: the sparse matrix of MFEM
   * is copied on the GPU for each operator, the preconditioner of MAGMA
   * being set up at the same time.
   *
   * A resolution is reported as converged if MAGMA reports a success and if
   * the residual recomputed by MAGMA at the end is finite and lower than
   * ten times the tolerance, `max(rtol ||b||, atol)`: the preconditioned
   * TFQMR of MAGMA 2.10 reports its breakdowns as successes.
   */
  struct MAGMAKrylovSolver final : public mfem::Solver,
                                   public LinearSolverConvergenceReport {
    /*!
     * \param[in] s: settings
     * \param[in] p: preconditioner of MAGMA (not owned), if any
     */
    MAGMAKrylovSolver(MAGMAKrylovSolverSettings s, MAGMAPreconditioner* p)
        : settings(std::move(s)),
          preconditioner(p),
          shared_queue(p != nullptr ? p->getQueue()
                                    : std::make_shared<MAGMAQueue>()),
          workspace(shared_queue->queue) {
      if (this->preconditioner != nullptr) {
        this->preconditioner->setKrylovSolver(this->getPreconditionedSolver());
      }
    }
    void SetOperator(const mfem::Operator& op) override {
      const auto* const a = dynamic_cast<const mfem::SparseMatrix*>(&op);
      this->release();
      if (a == nullptr) {
        mgis::raise(
            "MAGMAKrylovSolver::SetOperator: the operator is not a sparse "
            "matrix of MFEM");
      }
      auto* const q = this->shared_queue->queue;
      auto m = makeMAGMAHostMatrix(*a);
      if (this->preconditioner != nullptr) {
        this->preconditioner->setup(m);
      }
      checkMAGMA(
          magma_dmtransfer(m.wrap(q), &(this->matrix), Magma_CPU, Magma_DEV, q),
          "transfer of the matrix");
      this->has_matrix = true;
      this->height = this->width = static_cast<int>(m.size);
    }
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override {
      using Initialization = MAGMAUnknownsInitialization;
      this->converged = false;
      this->iterations = 0;
      if ((!this->has_matrix) || ((this->preconditioner != nullptr) &&
                                  (!this->preconditioner->isReady()))) {
        mgis::raise("MAGMAKrylovSolver::Mult: the operator is not set");
      }
      if ((b.Size() != this->Width()) || (x.Size() != this->Height())) {
        mgis::raise("MAGMAKrylovSolver::Mult: invalid sizes");
      }
      auto* const q = this->shared_queue->queue;
      auto v = MAGMADeviceVectors(b, x,
                                  this->iterative_mode
                                      ? Initialization::INITIAL_GUESS
                                      : Initialization::NULL_INITIAL_GUESS,
                                  this->workspace, q);
      const auto nb = magma_dnrm2(v.rhs.num_rows, v.rhs.dval, 1, q);
      if (!(nb > 0)) {
        // null right hand side: the solution is null (MAGMA computes
        // residuals relative to the norm of the right hand side)
        v.setNullUnknowns();
        v.finish();
        this->converged = std::isfinite(nb);
        return;
      }
      auto parameters = magma_d_solver_par{};
      parameters.solver = this->settings.solver;
      parameters.rtol = this->settings.relative_tolerance;
      parameters.atol = this->settings.absolute_tolerance;
      parameters.maxiter = this->settings.maximum_number_of_iterations;
      parameters.restart = this->settings.restart;
      const auto info = this->solve(v, parameters);
      v.finish();
      if ((info != MAGMA_SUCCESS) && (!isNonConvergence(info))) {
        checkMAGMA(info, (this->settings.name + " resolution").c_str());
      }
      const auto bound = 10 * std::max(this->settings.relative_tolerance * nb,
                                       this->settings.absolute_tolerance);
      this->converged = (info == MAGMA_SUCCESS) &&
                        std::isfinite(parameters.final_res) &&
                        (parameters.final_res <= bound);
      this->iterations = static_cast<int>(parameters.numiter);
      if (this->settings.verbosity > 0) {
        getOutputStream() << this->settings.name << ": "
                          << (this->converged ? "converged" : "not converged")
                          << " in " << this->iterations
                          << " iterations, residual " << parameters.init_res
                          << " -> " << parameters.final_res << '\n';
      }
    }
    [[nodiscard]] bool hasConverged() const noexcept override {
      return this->converged;
    }
    [[nodiscard]] std::optional<int> getNumberOfIterations()
        const noexcept override {
      return this->iterations;
    }
    [[nodiscard]] std::optional<real> getRelativeTolerance()
        const noexcept override {
      return this->settings.relative_tolerance;
    }
    [[nodiscard]] bool setRelativeTolerance(const real tol) noexcept override {
      this->settings.relative_tolerance = tol;
      return true;
    }
    //! \brief destructor
    ~MAGMAKrylovSolver() override { this->release(); }

   private:
    //! \return the preconditioned variant of the Krylov solver
    [[nodiscard]] magma_solver_type getPreconditionedSolver() const noexcept {
      switch (this->settings.solver) {
        case Magma_CG:
          return Magma_PCG;
        case Magma_BICGSTAB:
          return Magma_PBICGSTAB;
        case Magma_GMRES:
          return Magma_PGMRES;
        case Magma_IDR:
          return Magma_PIDR;
        case Magma_CGS:
          return Magma_PCGS;
        case Magma_QMR:
          return Magma_PQMR;
        case Magma_TFQMR:
          return Magma_PTFQMR;
        default:
          break;
      }
      return this->settings.solver;
    }
    /*!
     * \return the return code of the Krylov solver of MAGMA
     * \param[in, out] v: vectors
     * \param[in, out] s: parameters of the solver, and its results
     */
    [[nodiscard]] magma_int_t solve(MAGMADeviceVectors& v,
                                    magma_d_solver_par& s) const {
      auto* const q = this->shared_queue->queue;
      auto* const p = this->preconditioner != nullptr
                          ? this->preconditioner->getMAGMAPreconditioner()
                          : nullptr;
      const auto& a = this->matrix;
      auto* const x = &(v.unknowns);
      switch (this->settings.solver) {
        case Magma_CG:
          return p != nullptr ? magma_dpcg(a, v.rhs, x, &s, p, q)
                              : magma_dcg_res(a, v.rhs, x, &s, q);
        case Magma_BICGSTAB:
          return p != nullptr ? magma_dpbicgstab(a, v.rhs, x, &s, p, q)
                              : magma_dbicgstab(a, v.rhs, x, &s, q);
        case Magma_GMRES: {
          // the flexible GMRES of MAGMA always applies a preconditioner
          auto none = magma_d_preconditioner{};
          none.solver = Magma_NONE;
          return magma_dfgmres(a, v.rhs, x, &s, p != nullptr ? p : &none, q);
        }
        case Magma_IDR:
          return p != nullptr ? magma_dpidr(a, v.rhs, x, &s, p, q)
                              : magma_didr(a, v.rhs, x, &s, q);
        case Magma_CGS:
          return p != nullptr ? magma_dpcgs(a, v.rhs, x, &s, p, q)
                              : magma_dcgs(a, v.rhs, x, &s, q);
        case Magma_QMR:
          return p != nullptr ? magma_dpqmr(a, v.rhs, x, &s, p, q)
                              : magma_dqmr(a, v.rhs, x, &s, q);
        case Magma_TFQMR:
          return p != nullptr ? magma_dptfqmr(a, v.rhs, x, &s, p, q)
                              : magma_dtfqmr(a, v.rhs, x, &s, q);
        default:
          break;
      }
      mgis::raise("MAGMAKrylovSolver::solve: unsupported solver");
    }
    //! \brief release the matrix of MAGMA
    void release() noexcept {
      if (this->has_matrix) {
        magma_dmfree(&(this->matrix), this->shared_queue->queue);
        this->has_matrix = false;
      }
    }
    //! \brief settings
    MAGMAKrylovSolverSettings settings;
    //! \brief preconditioner of MAGMA (not owned), if any
    MAGMAPreconditioner* const preconditioner;
    //! \brief queue of MAGMA, shared with the preconditioner
    const std::shared_ptr<MAGMAQueue> shared_queue;
    //! \brief copies of the vectors of MFEM
    mutable MAGMAWorkspace workspace;
    //! \brief matrix on the GPU
    magma_d_matrix matrix = {};
    //! \brief if true, a matrix is set
    bool has_matrix = false;
    //! \brief if true, the last resolution converged
    mutable bool converged = false;
    //! \brief number of iterations of the last resolution
    mutable int iterations = 0;
  };  // end of struct MAGMAKrylovSolver

  //! \brief entry of a Krylov solver of MAGMA in the factory
  struct MAGMAKrylovSolverEntry {
    //! \brief name of the solver
    const char* name;
    //! \brief Krylov solver of MAGMA
    magma_solver_type solver;
    //! \brief option of the dimension of the Krylov or shadow space, if any
    const char* restart_option;
    //! \brief description
    const char* description;
  };  // end of struct MAGMAKrylovSolverEntry

  //! \return the Krylov solvers of MAGMA
  [[nodiscard]] static std::array<MAGMAKrylovSolverEntry, 7>
  getMAGMAKrylovSolverEntries() {
    return {{{"MAGMACG", Magma_CG, nullptr,
              "conjugate gradient of MAGMA (symmetric positive definite "
              "matrices)"},
             {"MAGMABiCGSTAB", Magma_BICGSTAB, nullptr, "BiCGSTAB of MAGMA"},
             {"MAGMAGMRES", Magma_GMRES, "KDim",
              "restarted flexible GMRES of MAGMA"},
             {"MAGMAIDR", Magma_IDR, "ShadowSpaceDimension",
              "induced dimension reduction method IDR(s) of MAGMA"},
             {"MAGMACGS", Magma_CGS, nullptr,
              "conjugate gradient squared of MAGMA"},
             {"MAGMAQMR", Magma_QMR, nullptr,
              "quasi minimal residual method of MAGMA"},
             {"MAGMATFQMR", Magma_TFQMR, nullptr,
              "transpose free quasi minimal residual method of MAGMA"}}};
  }  // end of getMAGMAKrylovSolverEntries

  /*!
   * \return the generator of a Krylov solver of MAGMA, whose options are the
   * ones of the iterative solvers (tolerances relative to the norm of the
   * right hand side, as MAGMA), `Preconditioner` (a preconditioner of
   * MAGMA), and `KDim` for GMRES (dimension of the Krylov space, 50 by
   * default) or `ShadowSpaceDimension` for IDR (dimension s of the shadow
   * space, 1 by default)
   * \param[in] e: entry of the solver
   * \param[in] d: description of the solver
   */
  [[nodiscard]] static LinearSolverFactory<false>::Generator
  buildMAGMAKrylovSolverGenerator(const MAGMAKrylovSolverEntry& e,
                                  const LinearSolverDescription& d) {
    return [e, d](Context& ctx, FiniteElementSpace<false>& fespace,
                  const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      if (!checkParameters(ctx, params, d.options)) {
        return {};
      }
      auto settings = MAGMAKrylovSolverSettings{};
      settings.name = e.name;
      settings.solver = e.solver;
      const auto olvl =
          get_if<int>(ctx, params, Problem::SolverVerbosityLevel, 0);
      const auto ortol =
          get_if<double>(ctx, params, Problem::SolverRelativeTolerance,
                         settings.relative_tolerance);
      const auto oatol =
          get_if<double>(ctx, params, Problem::SolverAbsoluteTolerance,
                         settings.absolute_tolerance);
      const auto omaxit =
          get_if<int>(ctx, params, Problem::SolverMaximumNumberOfIterations,
                      static_cast<int>(settings.maximum_number_of_iterations));
      // the defaults of MAGMA: 50 for GMRES, and IDR(1) (MAGMA reads the
      // dimension of the shadow space of IDR from its restart parameter, the
      // value 50 standing for its default dimension, 1)
      const auto default_restart = e.solver == Magma_IDR ? 1 : 50;
      const auto orestart =
          e.restart_option != nullptr
              ? get_if<int>(ctx, params, e.restart_option, default_restart)
              : std::optional<int>{default_restart};
      if (isInvalid(olvl) || isInvalid(ortol) || isInvalid(oatol) ||
          isInvalid(omaxit) || isInvalid(orestart)) {
        return {};
      }
      if ((*ortol < 0) || (*oatol < 0) || (*omaxit < 1) || (*orestart < 1)) {
        return ctx.registerErrorMessage(
            "invalid options for the solver '" + std::string{e.name} +
            "' (negative tolerance, or maximum number of iterations or "
            "dimension of space smaller than one)");
      }
      if ((e.solver == Magma_IDR) && (*orestart == 50)) {
        return ctx.registerErrorMessage(
            "the solver 'MAGMAIDR' does not accept a shadow space of "
            "dimension 50 (MAGMA reads it as the default dimension, 1)");
      }
      settings.verbosity = *olvl;
      settings.relative_tolerance = boundRelativeTolerance(*ortol);
      settings.absolute_tolerance = *oatol;
      settings.maximum_number_of_iterations = *omaxit;
      settings.restart = *orestart;
      auto oprec = makePreconditioner<false>(ctx, fespace, e.name, d, params);
      if (isInvalid(oprec)) {
        return {};
      }
      auto* const p = dynamic_cast<MAGMAPreconditioner*>(oprec->get());
      if ((*oprec != nullptr) && (p == nullptr)) {
        return ctx.registerErrorMessage("the solver '" + std::string{e.name} +
                                        "' requires a preconditioner of "
                                        "MAGMA");
      }
      // QMR applies the transposed preconditioner, for which MAGMA 2.10
      // applies the factors of the incomplete Cholesky factorization in the
      // wrong order (the ones of the preconditioner itself)
      if ((e.solver == Magma_QMR) && (p != nullptr) &&
          (p->getType() == Magma_ICC)) {
        return ctx.registerErrorMessage(
            "the solver 'MAGMAQMR' does not accept the preconditioner "
            "'MAGMAIC', whose transpose MAGMA applies wrongly: use "
            "'MAGMAILU', the same factorization for a symmetric matrix");
      }
      try {
        return LinearSolverHandler{
            std::make_unique<MAGMAKrylovSolver>(std::move(settings), p),
            std::move(*oprec)};
      } catch (std::exception& ex) {
        return ctx.registerErrorMessage(ex.what());
      }
    };
  }  // end of buildMAGMAKrylovSolverGenerator

  void declareMAGMALinearSolvers(LinearSolverFactory<false>& f) {
    for (const auto& e : getMAGMAKrylovSolverEntries()) {
      auto options = getIterativeSolverParametersList();
      options.push_back("Preconditioner");
      if (e.restart_option != nullptr) {
        options.push_back(e.restart_option);
      }
      const auto d = LinearSolverDescription{
          .library = "MAGMA",
          .kind = LinearSolverKind::ITERATIVE,
          .placement = LinearSolverPlacement::DEVICE,
          .preconditioners = AcceptedPreconditioners::SAME_LIBRARY,
          .options = std::move(options),
          .description = e.description};
      declare(f, e.name, d, buildMAGMAKrylovSolverGenerator(e, d));
    }
  }  // end of declareMAGMALinearSolvers

  template <bool parallel>
  void declareMAGMAPreconditioners(PreconditionerFactory<parallel>& f) {
    const auto entries =
        std::array<std::tuple<const char*, magma_solver_type, const char*>, 3>{
            {{"MAGMAJacobi", Magma_JACOBI, "Jacobi preconditioner of MAGMA"},
             {"MAGMAILU", Magma_ILU,
              "incomplete LU factorization without fill-in of MAGMA "
              "(cuSPARSE or hipSPARSE)"},
             {"MAGMAIC", Magma_ICC,
              "incomplete Cholesky factorization without fill-in of MAGMA "
              "(cuSPARSE or hipSPARSE), built from the lower triangle of the "
              "matrix"}}};
    for (const auto& [name, type, description] : entries) {
      const auto d = LinearSolverDescription{
          .library = "MAGMA",
          .kind = LinearSolverKind::PRECONDITIONER,
          .placement = LinearSolverPlacement::DEVICE,
          .preconditioners = AcceptedPreconditioners::NONE,
          .options = {},
          .description =
              parallel ? std::string{description} +
                             ", applied to the diagonal block of each process"
                       : std::string{description}};
      declare(f, name, d,
              [t = type](Context& ctx, FiniteElementSpace<parallel>&,
                         const Parameters& opts)
                  -> std::unique_ptr<LinearSolverPreconditioner> {
                if (!checkParameters(ctx, opts, std::vector<std::string>{})) {
                  return {};
                }
                try {
                  return std::make_unique<MAGMAPreconditioner>(t);
                } catch (std::exception& e) {
                  return ctx.registerErrorMessage(e.what());
                }
              });
    }
  }  // end of declareMAGMAPreconditioners

#ifdef MFEM_USE_MPI
  template void declareMAGMAPreconditioners<true>(PreconditionerFactory<true>&);
#endif /* MFEM_USE_MPI */
  template void declareMAGMAPreconditioners<false>(
      PreconditionerFactory<false>&);

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_MGIS_HAVE_MAGMA_SPARSE */
