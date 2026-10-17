/*!
 * \file   src/PETScLinearSolvers.cxx
 * \brief  This file implements the module of the linear solvers and of the
 * preconditioners of PETSc (see `LinearSolverModules.hxx`): a Krylov
 * solver (`KSP`) and a preconditioner (`PC`) of PETSc, configured by the
 * options database of PETSc.
 *
 * When MFEM computes on a GPU and PETSc is built for it (CUDA or HIP), the
 * solvers and the preconditioners of PETSc compute on the GPU, unless the
 * execution settings place the linear solvers on the host (see
 * `usePETScOnDevice`): the matrices are converted to the AIJ format of the
 * GPU (see `makeDeviceMatrix`), and the vectors of PETSc of the GPU wrap the
 * memory of the GPU of the vectors of MFEM, without copy (see
 * `createDeviceVector`; the shell preconditioners are applied to vectors of
 * MFEM copied on the GPU, see `applyOnDevice`). Otherwise, they compute on the
 * host: the matrices are copied in the AIJ format of the host (see
 * `makeHostMatrix`, which does not initialize the GPU), and, when MFEM
 * computes on a GPU, the vectors are copied from and to the vectors of MFEM
 * through the arrays of the host of PETSc.
 *
 * The vectors of PETSc of MFEM 4.10 are not used on the GPU: when the device
 * memory of MFEM is managed by Umpire (the default when MFEM is built with
 * Umpire), they are vectors of the host, which the matrices and the
 * preconditioners of PETSc of the GPU write on the GPU without the vectors
 * seeing it, and the Krylov solvers of PETSc saw null preconditioned
 * residuals (wrong solutions).
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC)

#include <memory>
#include <algorithm>
#include <tuple>
#include <string>
#include <vector>
#include <optional>
#include <utility>
#include <string_view>
#include "petscsys.h"
#include "petscksp.h"
#include "mfem/general/device.hpp"
#include "mfem/linalg/petsc.hpp"
#include "mfem/fem/pfespace.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverConvergenceReport.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "PETScSupport.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief Krylov solver of PETSc (`KSP`) configured by the options database
   * of PETSc. The operators are converted to the AIJ format of PETSc at each
   * call to `SetOperator`, so that the preconditioners of PETSc can be built
   * (they are thus built again at each call). A preconditioner of mfem-mgis
   * is seen by PETSc as a shell preconditioner, a preconditioner of PETSc
   * (`PETScPreconditioner`) becomes the preconditioner of the Krylov solver.
   *
   * The `iterative_mode` member of MFEM is given to PETSc at the
   * construction only (the solvers of mfem-mgis start from a null initial
   * guess). A resolution stopped by the maximum number of iterations
   * without norm (option `ksp_norm_type none`) is reported as converged, as
   * PETSc does.
   */
  struct PETScKrylovSolver final : public PETScOptions,
                                   public mfem::PetscLinearSolver,
                                   public LinearSolverConvergenceReport {
    /*!
     * \param[in] comm: communicator
     */
    explicit PETScKrylovSolver(MPI_Comm comm)
        : PETScOptions(comm, "ksp"),
          mfem::PetscLinearSolver(comm, this->prefix, false) {}
    /*!
     * \brief set a preconditioner of mfem-mgis: a preconditioner of PETSc
     * (`PETScPreconditioner`) becomes the preconditioner of the Krylov
     * solver, the others are seen by PETSc as shell preconditioners. This
     * method hides the one of MFEM, which is not virtual.
     *
     * The shell preconditioners are the ones of mfem-mgis, applied to
     * vectors of MFEM copied from and to the vectors of PETSc, on the GPU
     * when the solvers of PETSc compute on the GPU (see `applyOnDevice`),
     * through the arrays of the host of the vectors of PETSc otherwise.
     * \param[in] p: preconditioner
     */
    void SetPreconditioner(mfem::Solver& p) {
      if (dynamic_cast<mfem::PetscPreconditioner*>(&p) != nullptr) {
        mfem::PetscLinearSolver::SetPreconditioner(p);
        return;
      }
      this->shell_preconditioner = &p;
      auto pc = PC{};
      if ((KSPGetPC(static_cast<KSP>(*this), &pc) != 0) ||
          (PCSetType(pc, PCSHELL) != 0) || (PCShellSetContext(pc, this) != 0) ||
          (PCShellSetApply(pc, applyShellPreconditioner) != 0) ||
          (PCShellSetName(pc, "preconditioner of mfem-mgis") != 0)) {
        mgis::raise(
            "PETScKrylovSolver::SetPreconditioner: the shell preconditioner "
            "can't be set");
      }
    }
    void SetOperator(const mfem::Operator& op) override {
      // unlike the Krylov solvers of MFEM, the solver of PETSc does not give
      // the operator to a shell preconditioner (the operator of a
      // preconditioner of PETSc is set by PETSc)
      if (this->shell_preconditioner != nullptr) {
        this->shell_preconditioner->SetOperator(op);
      }
      // the operator is copied in the AIJ format of the GPU or of the host
      // (see the description of this file), referenced by the Krylov solver
      const auto a = usePETScOnDevice() ? makeDeviceMatrix(this->GetComm(), op)
                                        : makeHostMatrix(this->GetComm(), op);
      this->rigid_body_modes.attach(*a);
      mfem::PetscLinearSolver::SetOperator(*a);
    }
    /*!
     * \brief compute the rigid body modes given to the preconditioners of
     * PETSc (see `PETScRigidBodyModes`)
     * \param[in] fespace: finite element space
     */
    void initializeRigidBodyModes(FiniteElementSpace<true>& fespace) {
      this->rigid_body_modes.initialize(fespace);
    }
    /*!
     * \brief set the rigid body modes given to the preconditioners of PETSc
     * \param[in] m: rigid body modes
     */
    void setRigidBodyModes(const PETScRigidBodyModes& m) {
      this->rigid_body_modes = m;
    }
    /*!
     * \brief compute the fields given to a field split preconditioner of
     * PETSc (see `PETScComponentFields`)
     * \param[in] fespace: finite element space
     */
    void initializeComponentFields(FiniteElementSpace<true>& fespace) {
      this->component_fields.initialize(fespace);
    }
    /*!
     * \brief set the fields given to a field split preconditioner of PETSc
     * \param[in] f: fields
     */
    void setComponentFields(const PETScComponentFields& f) {
      this->component_fields = f;
    }
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override {
      this->markAsUsed();
      // the options of PETSc are applied before the fields are given to
      // a field split preconditioner, which must precede its setup
      this->Customize();
      auto pc = PC{};
      if (KSPGetPC(static_cast<KSP>(*this), &pc) != 0) {
        mgis::raise("PETScKrylovSolver::Mult: KSPGetPC failed");
      }
      this->component_fields.attach(pc);
      if (!mfem::Device::IsEnabled()) {
        mfem::PetscLinearSolver::Mult(b, x);
        return;
      }
      // see the description of this file
      auto* const ksp = static_cast<KSP>(*this);
      auto guess = PetscBool{};
      if (KSPGetInitialGuessNonzero(ksp, &guess) != 0) {
        mgis::raise(
            "PETScKrylovSolver::Mult: KSPGetInitialGuessNonzero failed");
      }
      auto vx = PETScVector{};
      auto vb = PETScVector{};
      if (usePETScOnDevice()) {
        createDeviceVector(vb, this->GetComm(), b.Read(), b.Size());
        createDeviceVector(vx, this->GetComm(),
                           guess == PETSC_TRUE ? x.ReadWrite() : x.Write(),
                           x.Size());
        // the right hand side is not modified
        if ((VecLockReadPush(vb.v) != 0) ||
            (KSPSolve(ksp, vb.v, vx.v) != 0) || (VecLockReadPop(vb.v) != 0)) {
          mgis::raise("PETScKrylovSolver::Mult: KSPSolve failed");
        }
        synchronizeDeviceVector(vx.v);
        this->finalizePETScAtExitOnce();
        return;
      }
      createHostVector(vb, this->GetComm(), b.Size());
      createHostVector(vx, this->GetComm(), x.Size());
      copyToPETScVector(vb.v, b);
      if (guess == PETSC_TRUE) {
        copyToPETScVector(vx.v, x);
      }
      if (KSPSolve(ksp, vb.v, vx.v) != 0) {
        mgis::raise("PETScKrylovSolver::Mult: KSPSolve failed");
      }
      copyFromPETScVector(x, vx.v);
    }
    [[nodiscard]] bool hasConverged() const noexcept override {
      // the accessors of MFEM are not constant
      return const_cast<PETScKrylovSolver*>(this)->GetConverged() > 0;
    }
    [[nodiscard]] std::optional<int> getNumberOfIterations()
        const noexcept override {
      return const_cast<PETScKrylovSolver*>(this)->GetNumIterations();
    }
    [[nodiscard]] std::optional<real> getRelativeTolerance()
        const noexcept override {
      auto rtol = PetscReal{};
      if (KSPGetTolerances(KSP(*this), &rtol, nullptr, nullptr, nullptr) != 0) {
        return {};
      }
      return static_cast<real>(rtol);
    }
    [[nodiscard]] bool setRelativeTolerance(const real tol) noexcept override {
      this->SetRelTol(tol);
      return true;
    }

   private:
    /*!
     * \brief apply the shell preconditioner
     * \return 0 on success, an error of PETSc otherwise (the exceptions do
     * not cross PETSc)
     * \param[in] pc: preconditioner of PETSc
     * \param[in] x: vector
     * \param[out] y: result
     */
    static PetscErrorCode applyShellPreconditioner(PC pc,
                                                   Vec x,
                                                   Vec y) noexcept {
      try {
        auto* d = static_cast<void*>(nullptr);
        if (PCShellGetContext(pc, &d) != 0) {
          return PETSC_ERR_LIB;
        }
        const auto& s = *static_cast<const PETScKrylovSolver*>(d);
        if (usePETScOnDevice()) {
          applyOnDevice(*(s.shell_preconditioner), s.shell_rhs,
                        s.shell_solution, x, y);
          return PETSC_SUCCESS;
        }
        auto n = PetscInt{};
        if (VecGetLocalSize(x, &n) != 0) {
          return PETSC_ERR_LIB;
        }
        s.shell_rhs.UseDevice(true);
        s.shell_rhs.SetSize(static_cast<int>(n));
        s.shell_solution.UseDevice(true);
        s.shell_solution.SetSize(static_cast<int>(n));
        copyFromPETScVector(s.shell_rhs, x);
        s.shell_preconditioner->Mult(s.shell_rhs, s.shell_solution);
        copyToPETScVector(y, s.shell_solution);
        return PETSC_SUCCESS;
      } catch (...) {
      }
      return PETSC_ERR_LIB;
    }  // end of applyShellPreconditioner
    //! \brief shell preconditioner, if any (not owned)
    mfem::Solver* shell_preconditioner = nullptr;
    //! \brief right hand side of the shell preconditioner
    mutable mfem::Vector shell_rhs;
    //! \brief result of the shell preconditioner
    mutable mfem::Vector shell_solution;
    //! \brief rigid body modes given to the preconditioners of PETSc
    PETScRigidBodyModes rigid_body_modes;
    //! \brief fields given to a field split preconditioner of PETSc
    PETScComponentFields component_fields;
  };  // end of struct PETScKrylovSolver

  //! \brief preconditioner of PETSc (`PC`) configured by the options database
  struct PETScPreconditioner final : public PETScOptions,
                                     public mfem::PetscPreconditioner {
    /*!
     * \param[in] comm: communicator
     */
    explicit PETScPreconditioner(MPI_Comm comm)
        : PETScOptions(comm, "pc"),
          mfem::PetscPreconditioner(comm, this->prefix) {}
    void SetOperator(const mfem::Operator& op) override {
      // the operator is copied in the AIJ format of the GPU or of the host
      // (see the description of this file), referenced by the preconditioner
      const auto a = usePETScOnDevice() ? makeDeviceMatrix(this->GetComm(), op)
                                        : makeHostMatrix(this->GetComm(), op);
      this->rigid_body_modes.attach(*a);
      mfem::PetscPreconditioner::SetOperator(*a);
    }
    /*!
     * \brief compute the rigid body modes given to the preconditioner of
     * PETSc (see `PETScRigidBodyModes`)
     * \param[in] fespace: finite element space
     */
    void initializeRigidBodyModes(FiniteElementSpace<true>& fespace) {
      this->rigid_body_modes.initialize(fespace);
    }
    /*!
     * \return the rigid body modes given to the preconditioner, used by a
     * Krylov solver of PETSc which it preconditions (the Krylov solver sets
     * the operators of its preconditioner itself)
     */
    [[nodiscard]] const PETScRigidBodyModes& getRigidBodyModes()
        const noexcept {
      return this->rigid_body_modes;
    }
    /*!
     * \brief compute the fields given to the preconditioner of PETSc, if
     * it is a field split preconditioner (see `PETScComponentFields`)
     * \param[in] fespace: finite element space
     */
    void initializeComponentFields(FiniteElementSpace<true>& fespace) {
      this->component_fields.initialize(fespace);
    }
    /*!
     * \return the fields given to the preconditioner, used by a Krylov
     * solver of PETSc which it preconditions
     */
    [[nodiscard]] const PETScComponentFields& getComponentFields()
        const noexcept {
      return this->component_fields;
    }
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override {
      this->markAsUsed();
      auto* const pc = static_cast<PC>(*this);
      // the options of PETSc are applied before the fields are given to
      // a field split preconditioner, which must precede its setup
      this->Customize();
      this->component_fields.attach(pc);
      if (!mfem::Device::IsEnabled()) {
        mfem::PetscPreconditioner::Mult(b, x);
        return;
      }
      // see the description of this file
      auto vx = PETScVector{};
      auto vb = PETScVector{};
      if (usePETScOnDevice()) {
        createDeviceVector(vb, this->GetComm(), b.Read(), b.Size());
        createDeviceVector(vx, this->GetComm(), x.Write(), x.Size());
        if ((VecLockReadPush(vb.v) != 0) || (PCSetUp(pc) != 0) ||
            (PCApply(pc, vb.v, vx.v) != 0) || (VecLockReadPop(vb.v) != 0)) {
          mgis::raise("PETScPreconditioner::Mult: PCApply failed");
        }
        synchronizeDeviceVector(vx.v);
        this->finalizePETScAtExitOnce();
        return;
      }
      createHostVector(vb, this->GetComm(), b.Size());
      createHostVector(vx, this->GetComm(), x.Size());
      copyToPETScVector(vb.v, b);
      if ((PCSetUp(pc) != 0) || (PCApply(pc, vb.v, vx.v) != 0)) {
        mgis::raise("PETScPreconditioner::Mult: PCApply failed");
      }
      copyFromPETScVector(x, vx.v);
    }

   private:
    //! \brief rigid body modes given to the preconditioner
    PETScRigidBodyModes rigid_body_modes;
    //! \brief fields given to a field split preconditioner
    PETScComponentFields component_fields;
  };  // end of struct PETScPreconditioner

  /*!
   * \brief where the solvers and the preconditioners of PETSc compute: where
   * the execution settings place the linear solvers if PETSc is built for the
   * GPU of MFEM (see `usePETScOnDevice`), on the host otherwise
   */
  static constexpr auto petsc_placement =
#if (defined(MFEM_USE_CUDA) && defined(PETSC_HAVE_CUDA)) || \
    (defined(MFEM_USE_HIP) && defined(PETSC_HAVE_HIP))
      LinearSolverPlacement::REQUESTED;
#else
      LinearSolverPlacement::HOST;
#endif

  //! \return the options of the Krylov solver of PETSc
  [[nodiscard]] static std::vector<std::string> getPETScKrylovSolverOptions() {
    auto options = getIterativeSolverParametersList();
    options.insert(options.end(),
                   {"Options", "Preconditioner", "RigidBodyModes", "Fields"});
    return options;
  }  // end of getPETScKrylovSolverOptions

  /*!
   * \return if the fields of the components of the unknowns are given to a
   * field split preconditioner of PETSc (option `Fields`: `Components`, or
   * `None`, the default), see `PETScComponentFields`
   * \param[in, out] ctx: execution context
   * \param[in] params: parameters
   */
  [[nodiscard]] static std::optional<bool> getComponentFieldsOption(
      Context& ctx, const Parameters& params) {
    const auto of =
        get_if<std::string>(ctx, params, "Fields", std::string{"None"});
    if (isInvalid(of)) {
      return {};
    }
    if (*of == "Components") {
      return true;
    }
    if (*of == "None") {
      return false;
    }
    return ctx.registerErrorMessage("invalid value '" + *of +
                                    "' of the option 'Fields' (expected "
                                    "'Components' or 'None')");
  }  // end of getComponentFieldsOption

  /*!
   * \return the generator of the Krylov solver of PETSc
   * \param[in] d: description of the solver
   */
  [[nodiscard]] static LinearSolverFactory<true>::Generator
  buildPETScKrylovSolverGenerator(const LinearSolverDescription d) {
    return [d](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      if (!checkParameters(ctx, params, d.options)) {
        return {};
      }
      const auto oopts =
          get_if<Parameters>(ctx, params, "Options", Parameters{});
      const auto orbm = get_if<bool>(ctx, params, "RigidBodyModes", true);
      const auto ofields = getComponentFieldsOption(ctx, params);
      if (isInvalid(oopts) || isInvalid(orbm) || isInvalid(ofields)) {
        return {};
      }
      initializePETSc();
      auto s = std::make_unique<PETScKrylovSolver>(fespace.GetComm());
      if (!s->setOptions(ctx, *oopts)) {
        return {};
      }
      auto& ksp = *s;
      // the options of mfem-mgis, overridden by the ones of PETSc
      if (contains(params, Problem::SolverRelativeTolerance)) {
        const auto ov =
            get<double>(ctx, params, Problem::SolverRelativeTolerance);
        if (isInvalid(ov)) {
          return {};
        }
        ksp.SetRelTol(boundRelativeTolerance(*ov));
      }
      if (contains(params, Problem::SolverAbsoluteTolerance)) {
        const auto ov =
            get<double>(ctx, params, Problem::SolverAbsoluteTolerance);
        if (isInvalid(ov)) {
          return {};
        }
        ksp.SetAbsTol(*ov);
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        const auto ov =
            get<int>(ctx, params, Problem::SolverMaximumNumberOfIterations);
        if (isInvalid(ov)) {
          return {};
        }
        ksp.SetMaxIter(*ov);
      }
      if (contains(params, Problem::SolverVerbosityLevel)) {
        const auto ov = get<int>(ctx, params, Problem::SolverVerbosityLevel);
        if (isInvalid(ov)) {
          return {};
        }
        ksp.SetPrintLevel(*ov);
      }
      // a preconditioner of mfem-mgis is seen by PETSc as a shell
      // preconditioner
      auto oprec = makePreconditioner<true>(ctx, fespace, "PETScLinearSolver",
                                            d, params);
      if (isInvalid(oprec)) {
        return {};
      }
      if (*oprec != nullptr) {
        // PETSc would replace the preconditioner by the one given by these
        // options
        for (const auto& [n, p] : *oopts) {
          if (n.starts_with("pc_")) {
            std::ignore = ctx.registerErrorMessage(
                "the option '" + n +
                "' of the preconditioner of PETSc can't be given to a solver "
                "with a preconditioner of mfem-mgis (the options of the "
                "'PETScPreconditioner' preconditioner are given by its "
                "'Options' option)");
            return {};
          }
        }
        s->SetPreconditioner(**oprec);
      }
      // rigid body modes and fields: a preconditioner of PETSc shares its
      // object of PETSc with the solver, which sets its operators: its own
      // options apply
      const auto* const pp =
          dynamic_cast<const PETScPreconditioner*>(oprec->get());
      if (pp != nullptr) {
        s->setRigidBodyModes(pp->getRigidBodyModes());
        s->setComponentFields(pp->getComponentFields());
      } else {
        if (*orbm) {
          s->initializeRigidBodyModes(fespace);
        }
        if (*ofields) {
          s->initializeComponentFields(fespace);
        }
      }
      return LinearSolverHandler{std::move(s), std::move(*oprec)};
    };
  }  // end of buildPETScKrylovSolverGenerator

  /*!
   * \return the preconditioner of PETSc
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options (`Options`, the options of PETSc;
   * `RigidBodyModes`, true by default, see `PETScRigidBodyModes`; `Fields`,
   * see `getComponentFieldsOption`)
   */
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  makePETScPreconditioner(Context& ctx,
                          FiniteElementSpace<true>& fespace,
                          const Parameters& opts) {
    if (!checkParameters(
            ctx, opts,
            std::vector<std::string>{"Options", "RigidBodyModes", "Fields"})) {
      return {};
    }
    const auto oopts = get_if<Parameters>(ctx, opts, "Options", Parameters{});
    const auto orbm = get_if<bool>(ctx, opts, "RigidBodyModes", true);
    const auto ofields = getComponentFieldsOption(ctx, opts);
    if (isInvalid(oopts) || isInvalid(orbm) || isInvalid(ofields)) {
      return {};
    }
    initializePETSc();
    auto p = std::make_unique<PETScPreconditioner>(fespace.GetComm());
    if (!p->setOptions(ctx, *oopts)) {
      return {};
    }
    if (*orbm) {
      p->initializeRigidBodyModes(fespace);
    }
    if (*ofields) {
      p->initializeComponentFields(fespace);
    }
    return p;
  }  // end of makePETScPreconditioner

  void declarePETScLinearSolvers(LinearSolverFactory<true>& f) {
    const auto d = LinearSolverDescription{
        .library = "PETSc",
        .kind = LinearSolverKind::ITERATIVE,
        .placement = petsc_placement,
        .preconditioners = AcceptedPreconditioners::ANY,
        .options = getPETScKrylovSolverOptions(),
        .description =
            "Krylov solver of PETSc (KSP), configured by the "
            "options of PETSc given by the 'Options' parameter"};
    declare(f, "PETScLinearSolver", d, buildPETScKrylovSolverGenerator(d));
  }  // end of declarePETScLinearSolvers

  void declarePETScPreconditioners(PreconditionerFactory<true>& f) {
    declare(f, "PETScPreconditioner",
            LinearSolverDescription{
                .library = "PETSc",
                .kind = LinearSolverKind::PRECONDITIONER,
                .placement = petsc_placement,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Options", "RigidBodyModes", "Fields"},
                .description = "preconditioner of PETSc (PC), configured by "
                               "the options of PETSc given by the 'Options' "
                               "parameter"},
            makePETScPreconditioner);
  }  // end of declarePETScPreconditioners

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC) */
