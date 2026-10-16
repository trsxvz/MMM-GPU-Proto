/*!
 * \file   src/TrilinosLinearSolvers.cxx
 * \brief  This file implements the module of the linear solvers and of the
 * preconditioners of Trilinos (see `LinearSolverModules.hxx`): Krylov
 * solvers of Belos, direct solvers of Amesos2, preconditioners of Ifpack2
 * and the algebraic multigrid of MueLu, computed by the library
 * `MFEMMGISTrilinos` (see `TrilinosBackend.hxx`)
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_TRILINOS)

#include <array>
#include <tuple>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <optional>
#include <algorithm>
#include <exception>
#include <variant>
#include <type_traits>
#include "mfem/general/forall.hpp"
#include "mfem/general/device.hpp"
#include "mfem/general/communication.hpp"
#include "mfem/linalg/hypre.hpp"
#include "mfem/fem/pfespace.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverConvergenceReport.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "KokkosRuntime.hxx"
#include "NodesCoordinates.hxx"
#include "TrilinosBackend.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis::internals {

  static_assert(std::is_same_v<real, double>,
                "the module of Trilinos requires the double precision");
  static_assert(std::is_same_v<HYPRE_Int, int>,
                "the module of Trilinos does not support hypre built with "
                "64 bits integers (--enable-bigint)");

  /*!
   * \return the options of Trilinos given by parameters, or an empty value
   * on error: booleans, integers, reals, strings, and parameters as
   * sublists
   * \param[in, out] ctx: execution context
   * \param[in] params: parameters
   */
  [[nodiscard]] static std::optional<trilinos::Options>
  convertToTrilinosOptions(Context& ctx, const Parameters& params) {
    auto options = trilinos::Options{};
    for (const auto& [n, p] : params) {
      auto o = trilinos::Option{};
      if (is<bool>(p)) {
        o.value = *get<bool>(ctx, p);
      } else if (is<int>(p)) {
        o.value = *get<int>(ctx, p);
      } else if (is<double>(p)) {
        o.value = *get<double>(ctx, p);
      } else if (is<std::string>(p)) {
        o.value = *get<std::string>(ctx, p);
      } else if (is<Parameters>(p)) {
        auto osub = convertToTrilinosOptions(ctx, *get<Parameters>(ctx, p));
        if (isInvalid(osub)) {
          return {};
        }
        o.value = std::move(*osub);
      } else {
        std::ignore = ctx.registerErrorMessage(
            "invalid value of the parameter '" + n +
            "' of Trilinos (a boolean, an integer, a real, a string or a list "
            "of parameters is expected)");
        return {};
      }
      options.emplace_back(n, std::move(o));
    }
    return options;
  }  // end of convertToTrilinosOptions

  /*!
   * \brief view of a matrix of hypre on the host, given to Trilinos
   */
  struct HypreMatrixView {
    /*!
     * \param[in] a: matrix of hypre
     * \param[in] vdim: number of components of the unknowns
     * \param[in] by_nodes: if true, the unknowns are ordered by nodes
     * \param[in] gpu: if true, the solver of Trilinos computes on the GPU:
     * the arrays of hypre are read there when hypre computes on the GPU
     */
    HypreMatrixView(const mfem::HypreParMatrix& a,
                    const int vdim,
                    const bool by_nodes,
                    const bool gpu) {
      // the column of index `j` of the diagonal block is the row of index
      // `first_row + j`: the matrix must be square, its columns being
      // distributed as its rows (checked by all the processes, which fail
      // together)
      auto valid =
          static_cast<int>((a.Height() == a.Width()) &&
                           (a.GetGlobalNumRows() == a.GetGlobalNumCols()) &&
                           (a.GetRowStarts()[0] == a.GetColStarts()[0]) &&
                           (a.Height() % vdim == 0));
      MPI_Allreduce(MPI_IN_PLACE, &valid, 1, MPI_INT, MPI_LAND, a.GetComm());
      if (valid == 0) {
        mgis::raise(
            "Trilinos: the matrix must be square, its columns being "
            "distributed as its rows, the number of rows of each process "
            "being a multiple of the number of components");
      }
      // the arrays of hypre are read where Trilinos computes, if hypre
      // computes there, and on the host otherwise (the indices of the
      // columns of the off-diagonal block are on the host)
      const auto device = gpu && isHypreUsingGPU();
      if (device) {
        a.HypreRead();
      } else {
        a.HostRead();
      }
      auto* const h = static_cast<hypre_ParCSRMatrix*>(
          const_cast<mfem::HypreParMatrix&>(a));
      auto* const diag = hypre_ParCSRMatrixDiag(h);
      auto* const offd = hypre_ParCSRMatrixOffd(h);
      const auto* const cmap = hypre_ParCSRMatrixColMapOffd(h);
      const auto no = hypre_CSRMatrixNumCols(offd);
      this->offd_columns.assign(cmap, cmap + no);
      this->view.communicator = a.GetComm();
      this->view.global_size =
          static_cast<long long>(hypre_ParCSRMatrixGlobalNumRows(h));
      this->view.first_row =
          static_cast<long long>(hypre_ParCSRMatrixFirstRowIndex(h));
      this->view.local_size = hypre_CSRMatrixNumRows(diag);
      this->view.diag_offsets = hypre_CSRMatrixI(diag);
      this->view.diag_columns = hypre_CSRMatrixJ(diag);
      this->view.diag_values = hypre_CSRMatrixData(diag);
      this->view.offd_offsets = hypre_CSRMatrixI(offd);
      this->view.offd_columns = hypre_CSRMatrixJ(offd);
      this->view.offd_values = hypre_CSRMatrixData(offd);
      this->view.number_of_offd_columns = no;
      this->view.offd_global_columns = this->offd_columns.data();
      this->view.number_of_components = vdim;
      this->view.ordered_by_nodes = by_nodes;
      this->view.device = device;
    }
    //! \brief view
    trilinos::DistributedMatrixView view;

   private:
    //! \brief global indices of the columns of the off-diagonal block
    std::vector<long long> offd_columns;
  };  // end of struct HypreMatrixView

  /*!
   * \brief linear solver or preconditioner of Trilinos, seen as a solver of
   * MFEM: the matrix of hypre is given to Trilinos on the host (Tpetra
   * copies it on the GPU when it computes there), and the vectors of MFEM on
   * the GPU when MFEM and Tpetra compute there, on the host otherwise. The
   * errors of Trilinos are reported by exceptions.
   */
  struct TrilinosSolver final : public mfem::Solver,
                                public LinearSolverConvergenceReport {
    /*!
     * \param[in] o: solver or preconditioner of Trilinos
     * \param[in] fespace: finite element space
     * \param[in] is_iterative: if true, the solver is iterative (and its
     * convergence is reported)
     */
    TrilinosSolver(std::unique_ptr<trilinos::Operator> o,
                   const FiniteElementSpace<true>& fespace,
                   const bool is_iterative)
        : backend(std::move(o)),
          number_of_components(fespace.GetVDim()),
          ordered_by_nodes(fespace.GetOrdering() == mfem::Ordering::byNODES),
          iterative(is_iterative) {}
    //! \return the solver or the preconditioner of Trilinos
    [[nodiscard]] const trilinos::Operator& getBackend() const noexcept {
      return *(this->backend);
    }
    /*!
     * \brief set the preconditioner given to a Krylov solver of Belos,
     * whose operator is set with the one of the solver
     * \param[in] p: preconditioner
     */
    void setPreconditioner(mfem::Solver& p) noexcept {
      this->preconditioner = &p;
    }
    void SetOperator(const mfem::Operator& op) override {
      const auto* const a = dynamic_cast<const mfem::HypreParMatrix*>(&op);
      if (a == nullptr) {
        mgis::raise(
            "TrilinosSolver::SetOperator: the solvers and the "
            "preconditioners of Trilinos require a matrix of hypre");
      }
      if (this->preconditioner != nullptr) {
        this->preconditioner->SetOperator(op);
      }
      this->height = a->Height();
      this->width = a->Width();
      // the matrix, read on the host, is given back to hypre, whose
      // computations may be on the GPU, even on error (its values are
      // copied by Trilinos)
      struct HypreLocationGuard {
        ~HypreLocationGuard() { this->matrix.HypreRead(); }
        const mfem::HypreParMatrix& matrix;
      };
      const auto guard = HypreLocationGuard{*a};
      const auto v = HypreMatrixView(*a, this->number_of_components,
                                     this->ordered_by_nodes,
                                     this->backend->computesOnTheGPU());
      this->backend->setMatrix(v.view);
    }
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override {
      // the vectors are given on the GPU when MFEM and Tpetra compute there
      const auto on_device = this->backend->computesOnTheGPU() &&
                             b.UseDevice() && x.UseDevice() &&
                             mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                                                  mfem::Backend::HIP_MASK);
      if (!this->iterative_mode) {
        x = 0.;
      }
      const auto* const pb = b.Read(on_device);
      auto* const px = x.ReadWrite(on_device);
      if (on_device) {
        // Tpetra may compute on another stream than MFEM
        MFEM_DEVICE_SYNC;
        this->backend->applyOnDevice(pb, px);
      } else {
        this->backend->apply(pb, px);
      }
      // a failure of the solver is a resolution which did not converge
      const auto& failure = this->backend->getReport().failure;
      if ((!failure.empty()) && (mfem::Mpi::Root())) {
        mfem_mgis::getErrorStream()
            << "Trilinos: the resolution failed: " << failure << '\n';
      }
    }
    [[nodiscard]] bool hasConverged() const noexcept override {
      return this->backend->getReport().converged;
    }
    [[nodiscard]] std::optional<int> getNumberOfIterations()
        const noexcept override {
      if (!this->iterative) {
        return {};
      }
      return this->backend->getReport().iterations;
    }
    [[nodiscard]] std::optional<real> getRelativeTolerance()
        const noexcept override {
      try {
        const auto tol = this->backend->getRelativeTolerance();
        if (tol >= 0) {
          return tol;
        }
      } catch (...) {
      }
      return {};
    }
    [[nodiscard]] bool setRelativeTolerance(const real tol) noexcept override {
      try {
        return this->backend->setRelativeTolerance(tol);
      } catch (...) {
      }
      return false;
    }

   private:
    //! \brief solver or preconditioner of Trilinos
    std::unique_ptr<trilinos::Operator> backend;
    //! \brief preconditioner of a Krylov solver of Belos (not owned)
    mfem::Solver* preconditioner = nullptr;
    //! \brief number of components of the unknowns
    const int number_of_components;
    //! \brief if true, the unknowns are ordered by nodes
    const bool ordered_by_nodes;
    //! \brief if true, the solver is iterative
    const bool iterative;
  };  // end of struct TrilinosSolver

  /*!
   * \return true if the solvers of Trilinos compute on the GPU: Trilinos is
   * built for a GPU and MFEM computes on the GPU. Otherwise, they compute
   * on the host, when Tpetra is built for a node of the host (see
   * `trilinos::hasGPUNode`).
   */
  [[nodiscard]] static bool computeTrilinosOnTheGPU() {
    return trilinos::hasGPUNode() &&
           mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                                mfem::Backend::HIP_MASK);
  }  // end of computeTrilinosOnTheGPU

  /*!
   * \return the preconditioner of mfem-mgis given to Belos: its vectors are
   * the ones of MFEM, given to Trilinos on the GPU when MFEM and Tpetra
   * compute there, and on the host otherwise
   * \param[in] p: preconditioner of mfem-mgis
   * \param[in] gpu: if true, the solver of Belos computes on the GPU
   */
  [[nodiscard]] static trilinos::ExternalPreconditioner
  makeExternalPreconditioner(mfem::Solver& p, const bool gpu) {
    const auto device = gpu && trilinos::hasGPUNode() &&
                        mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                                             mfem::Backend::HIP_MASK);
    // right hand side and result
    const auto v = std::make_shared<std::array<mfem::Vector, 2>>();
    auto e = trilinos::ExternalPreconditioner{};
    e.device = device;
    e.getRightHandSide = [&p, v, device] {
      auto& b = (*v)[0];
      b.UseDevice(true);
      b.SetSize(p.Height());
      return device ? b.Write() : b.HostWrite();
    };
    e.apply = [&p, v, device] {
      auto& [b, x] = *v;
      x.UseDevice(true);
      x.SetSize(p.Height());
      x = 0.;
      p.Mult(b, x);
      const auto* const r = device ? x.Read() : x.HostRead();
      if (device) {
        // Tpetra may compute on another stream than MFEM
        MFEM_DEVICE_SYNC;
      }
      return r;
    };
    return e;
  }  // end of makeExternalPreconditioner

  /*!
   * \return the space dimension and, for a finite element space whose
   * number of components is the space dimension (2 or 3), the coordinates
   * of its nodes (see `NodesCoordinates.hxx`)
   * \param[in] fespace: finite element space
   */
  [[nodiscard]] static trilinos::NodesCoordinates getTrilinosNodesCoordinates(
      FiniteElementSpace<true>& fespace) {
    auto c = trilinos::NodesCoordinates{};
    c.dimension = fespace.GetParMesh()->SpaceDimension();
    if (!hasRigidBodyModes(fespace)) {
      return c;
    }
    const auto nc = getNodesCoordinates(fespace);
    c.values.assign(nc.values.begin(), nc.values.end());
    return c;
  }  // end of getTrilinosNodesCoordinates

  /*!
   * \return the options of a solver of Trilinos given by the parameter
   * `Parameters`, or an empty value on error
   * \param[in, out] ctx: execution context
   * \param[in] params: parameters
   */
  [[nodiscard]] static std::optional<trilinos::Options> getTrilinosOptions(
      Context& ctx, const Parameters& params) {
    const auto o = get_if<Parameters>(ctx, params, "Parameters", Parameters{});
    if (isInvalid(o)) {
      return {};
    }
    return convertToTrilinosOptions(ctx, *o);
  }  // end of getTrilinosOptions

  /*!
   * \brief create a solver of Trilinos, Kokkos being initialized if
   * required, the errors being reported in the execution context
   * \return the solver, or an empty value on error
   * \param[in, out] ctx: execution context
   * \param[in] n: name of the solver, for the error messages
   * \param[in] f: function creating the solver
   */
  template <typename Function>
  [[nodiscard]] static std::optional<std::unique_ptr<trilinos::Operator>>
  callTrilinos(Context& ctx, const std::string& n, const Function& f) {
    try {
      initializeKokkosOnDemand();
      return f();
    } catch (std::exception& e) {
      std::ignore = ctx.registerErrorMessage(n + ": " + e.what());
    }
    return {};
  }  // end of callTrilinos

  /*!
   * \return the placement of the computations of Trilinos: where MFEM
   * computes when Trilinos is built for a GPU and for a node of the host,
   * on the GPU when it is built for the GPU only, on the host otherwise
   */
  [[nodiscard]] static LinearSolverPlacement getTrilinosPlacement() {
    if (!trilinos::hasGPUNode()) {
      return LinearSolverPlacement::HOST;
    }
    return trilinos::hasHostNode() ? LinearSolverPlacement::MFEM
                                   : LinearSolverPlacement::DEVICE;
  }  // end of getTrilinosPlacement

  //! \return the options of the Krylov solvers of Belos
  [[nodiscard]] static std::vector<std::string> getBelosSolverOptions() {
    auto options = getIterativeSolverParametersList();
    options.insert(options.end(), {"Solver", "Parameters", "Preconditioner"});
    return options;
  }  // end of getBelosSolverOptions

  /*!
   * \return the generator of the Krylov solvers of Belos
   * \param[in] d: description of the solvers
   */
  [[nodiscard]] static LinearSolverFactory<true>::Generator
  buildBelosSolverGenerator(const LinearSolverDescription d) {
    return [d](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      if (!checkParameters(ctx, params, d.options)) {
        return {};
      }
      const auto osolver =
          get_if<std::string>(ctx, params, "Solver", std::string{"GMRES"});
      const auto ooptions = getTrilinosOptions(ctx, params);
      if (isInvalid(osolver) || isInvalid(ooptions)) {
        return {};
      }
      // the settings of mfem-mgis, overridden by the parameters of Belos
      auto settings = trilinos::BelosSettings{};
      if (contains(params, Problem::SolverRelativeTolerance)) {
        const auto ov =
            get<double>(ctx, params, Problem::SolverRelativeTolerance);
        if (isInvalid(ov)) {
          return {};
        }
        settings.relative_tolerance =
            static_cast<double>(boundRelativeTolerance(*ov));
      }
      if (contains(params, Problem::SolverAbsoluteTolerance)) {
        const auto ov =
            get<double>(ctx, params, Problem::SolverAbsoluteTolerance);
        if (isInvalid(ov)) {
          return {};
        }
        if ((*ov < 0) || (*ov > 0)) {
          std::ignore = ctx.registerErrorMessage(
              "BelosSolver: the solvers of Belos only have a relative "
              "tolerance (the absolute tolerance must be null)");
          return {};
        }
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        const auto ov =
            get<int>(ctx, params, Problem::SolverMaximumNumberOfIterations);
        if (isInvalid(ov)) {
          return {};
        }
        settings.maximum_number_of_iterations = *ov;
      }
      if (contains(params, Problem::SolverVerbosityLevel)) {
        const auto ov = get<int>(ctx, params, Problem::SolverVerbosityLevel);
        if (isInvalid(ov)) {
          return {};
        }
        settings.verbosity_level = *ov;
      }
      auto oprec =
          makePreconditioner<true>(ctx, fespace, "BelosSolver", d, params);
      if (isInvalid(oprec)) {
        return {};
      }
      // a preconditioner of Trilinos is applied to the vectors of Tpetra,
      // the other ones to vectors of MFEM
      const auto gpu = computeTrilinosOnTheGPU();
      const auto* tp = static_cast<const trilinos::Operator*>(nullptr);
      auto external = trilinos::ExternalPreconditioner{};
      if (*oprec != nullptr) {
        const auto* const t = dynamic_cast<const TrilinosSolver*>(oprec->get());
        if (t != nullptr) {
          tp = &(t->getBackend());
        } else {
          external = makeExternalPreconditioner(**oprec, gpu);
        }
      }
      auto obackend = callTrilinos(ctx, "BelosSolver", [&] {
        return trilinos::makeBelosSolver(*osolver, *ooptions, settings, tp,
                                         std::move(external), gpu);
      });
      if (isInvalid(obackend)) {
        return {};
      }
      auto s =
          std::make_unique<TrilinosSolver>(std::move(*obackend), fespace, true);
      if (*oprec != nullptr) {
        s->setPreconditioner(**oprec);
      }
      return LinearSolverHandler{std::move(s), std::move(*oprec)};
    };
  }  // end of buildBelosSolverGenerator

  /*!
   * \return the direct solver of Amesos2
   * \param[in, out] ctx: execution context
   * \param[in] params: parameters (`Solver`, `Parameters`)
   */
  [[nodiscard]] static std::optional<std::unique_ptr<trilinos::Operator>>
  makeAmesos2(Context& ctx, const Parameters& params) {
    if (!checkParameters(ctx, params,
                         std::vector<std::string>{"Solver", "Parameters"})) {
      return {};
    }
    const auto osolver =
        get_if<std::string>(ctx, params, "Solver", std::string{"KLU2"});
    const auto ooptions = getTrilinosOptions(ctx, params);
    if (isInvalid(osolver) || isInvalid(ooptions)) {
      return {};
    }
    return callTrilinos(ctx, "Amesos2", [&] {
      if (!trilinos::isAmesos2SolverAvailable(*osolver)) {
        mgis::raise("the solver '" + *osolver + "' is not available");
      }
      return trilinos::makeAmesos2Solver(*osolver, *ooptions,
                                         computeTrilinosOnTheGPU());
    });
  }  // end of makeAmesos2

  void declareTrilinosLinearSolvers(LinearSolverFactory<true>& f) {
    const auto belos = LinearSolverDescription{
        .library = "Trilinos",
        .kind = LinearSolverKind::ITERATIVE,
        .placement = getTrilinosPlacement(),
        .preconditioners = AcceptedPreconditioners::ANY,
        .options = getBelosSolverOptions(),
        .description =
            "Krylov solver of Belos (Trilinos), selected by the 'Solver' "
            "parameter (GMRES by default, CG, BICGSTAB, GCRODR, MINRES, "
            "TFQMR...), with the parameters of Belos given by the "
            "'Parameters' parameter, preconditioned on the right"};
    declare(f, "BelosSolver", belos, buildBelosSolverGenerator(belos));
    declare(f, "Amesos2Solver",
            LinearSolverDescription{
                .library = "Trilinos",
                .kind = LinearSolverKind::DIRECT,
                .placement = getTrilinosPlacement(),
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Solver", "Parameters"},
                .description = "direct solver of Amesos2 (Trilinos), "
                               "selected by the 'Solver' parameter (KLU2 by "
                               "default, Basker, Tacho...)"},
            [](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& params) -> LinearSolverHandler {
              auto o = makeAmesos2(ctx, params);
              if (isInvalid(o)) {
                return {};
              }
              return LinearSolverHandler{std::make_unique<TrilinosSolver>(
                                             std::move(*o), fespace, false),
                                         nullptr};
            });
  }  // end of declareTrilinosLinearSolvers

  void declareTrilinosPreconditioners(PreconditionerFactory<true>& f) {
    const auto placement = getTrilinosPlacement();
    declare(f, "Ifpack2",
            LinearSolverDescription{
                .library = "Trilinos",
                .kind = LinearSolverKind::PRECONDITIONER,
                .placement = placement,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Type", "Parameters"},
                .description =
                    "preconditioner of Ifpack2 (Trilinos), selected by the "
                    "'Type' parameter (RELAXATION by default, CHEBYSHEV, "
                    "RILUK, ILUT, SCHWARZ...)"},
            [](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& params)
                -> std::unique_ptr<LinearSolverPreconditioner> {
              if (!checkParameters(
                      ctx, params,
                      std::vector<std::string>{"Type", "Parameters"})) {
                return {};
              }
              const auto otype = get_if<std::string>(ctx, params, "Type",
                                                     std::string{"RELAXATION"});
              const auto ooptions = getTrilinosOptions(ctx, params);
              if (isInvalid(otype) || isInvalid(ooptions)) {
                return {};
              }
              auto o = callTrilinos(ctx, "Ifpack2", [&] {
                return trilinos::makeIfpack2Preconditioner(
                    *otype, *ooptions, computeTrilinosOnTheGPU());
              });
              if (isInvalid(o)) {
                return {};
              }
              return std::make_unique<TrilinosSolver>(std::move(*o), fespace,
                                                      false);
            });
    declare(f, "MueLu",
            LinearSolverDescription{
                .library = "Trilinos",
                .kind = LinearSolverKind::PRECONDITIONER,
                .placement = placement,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Parameters"},
                .description =
                    "smoothed aggregation algebraic multigrid of MueLu "
                    "(Trilinos), one V-cycle, with the rigid body modes as "
                    "near null space when the number of components of the "
                    "unknowns is the space dimension"},
            [](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& params)
                -> std::unique_ptr<LinearSolverPreconditioner> {
              if (!checkParameters(ctx, params,
                                   std::vector<std::string>{"Parameters"})) {
                return {};
              }
              const auto ooptions = getTrilinosOptions(ctx, params);
              if (isInvalid(ooptions)) {
                return {};
              }
              // collective computation, before the calls to Trilinos,
              // which may fail on one process only
              auto c = getTrilinosNodesCoordinates(fespace);
              auto o = callTrilinos(ctx, "MueLu", [&] {
                return trilinos::makeMueLuPreconditioner(
                    *ooptions, std::move(c), computeTrilinosOnTheGPU());
              });
              if (isInvalid(o)) {
                return {};
              }
              return std::make_unique<TrilinosSolver>(std::move(*o), fespace,
                                                      false);
            });
    // FROSch computes with the default node of Tpetra
    declare(f, "FROSch",
            LinearSolverDescription{
                .library = "Trilinos",
                .kind = LinearSolverKind::PRECONDITIONER,
                .placement = trilinos::hasGPUNode()
                                 ? LinearSolverPlacement::DEVICE
                                 : LinearSolverPlacement::HOST,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Parameters"},
                .description =
                    "two level overlapping Schwarz preconditioner of FROSch "
                    "(Trilinos), with the null space of the linear "
                    "elasticity when the number of components of the "
                    "unknowns is the space dimension"},
            [](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& params)
                -> std::unique_ptr<LinearSolverPreconditioner> {
              if (!checkParameters(ctx, params,
                                   std::vector<std::string>{"Parameters"})) {
                return {};
              }
              const auto ooptions = getTrilinosOptions(ctx, params);
              if (isInvalid(ooptions)) {
                return {};
              }
              // collective computation, before the calls to Trilinos,
              // which may fail on one process only
              auto c = getTrilinosNodesCoordinates(fespace);
              auto o = callTrilinos(ctx, "FROSch", [&] {
                return trilinos::makeFROSchPreconditioner(
                    *ooptions, std::move(c), computeTrilinosOnTheGPU());
              });
              if (isInvalid(o)) {
                return {};
              }
              return std::make_unique<TrilinosSolver>(std::move(*o), fespace,
                                                      false);
            });
    declare(f, "Amesos2Preconditioner",
            LinearSolverDescription{
                .library = "Trilinos",
                .kind = LinearSolverKind::PRECONDITIONER,
                .placement = placement,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Solver", "Parameters"},
                .description = "factorization of a direct solver of "
                               "Amesos2 (Trilinos), applied once"},
            [](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& params)
                -> std::unique_ptr<LinearSolverPreconditioner> {
              auto o = makeAmesos2(ctx, params);
              if (isInvalid(o)) {
                return {};
              }
              return std::make_unique<TrilinosSolver>(std::move(*o), fespace,
                                                      false);
            });
  }  // end of declareTrilinosPreconditioners

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_TRILINOS) */
