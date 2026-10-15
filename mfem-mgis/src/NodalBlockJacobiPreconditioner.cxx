/*!
 * \file   src/NodalBlockJacobiPreconditioner.cxx
 * \brief  This file implements the Jacobi preconditioner by nodal blocks,
 * whose blocks are inverted and applied by the batched dense linear algebra
 * of MFEM (see `LinearSolverModules.hxx`)
 * \date   15/10/2026
 */

#include <array>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <optional>
#include <string_view>
#include <type_traits>
#include "mfem/general/backends.hpp"
#include "mfem/general/device.hpp"
#include "mfem/linalg/densemat.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/linalg/batched/batched.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/MatrixFreeTangentOperator.hxx"
#include "LinearSolverModules.hxx"
#include "NodalBlockJacobiKernels.hxx"

namespace mfem_mgis::internals {

  //! \brief settings of the Jacobi preconditioner by nodal blocks
  struct NodalBlockJacobiSettings {
    //! \brief backend of the batched linear algebra of MFEM
    mfem::BatchedLinAlg::Backend backend = mfem::BatchedLinAlg::NATIVE;
    //! \brief size of the blocks
    int block_size = 1;
    //! \brief if true, the unknowns are ordered by nodes (`byNODES`)
    bool ordered_by_nodes = true;
    /*!
     * \brief if true, the blocks are factorized (LU factors applied by
     * batched solves), otherwise inverted (inverses applied by batched
     * products)
     */
    bool factorize = false;
    //! \brief damping factor, if any
    std::optional<real> damping;
  };

  /*!
   * \brief Jacobi preconditioner by nodal blocks: the blocks coupling the
   * components of each node, extracted from the matrix (from the diagonal
   * block of each process in parallel), are inverted or factorized by the
   * batched dense linear algebra of MFEM (native kernels, cuBLAS or
   * hipBLAS, or MAGMA), which computes where MFEM computes. The vectors are
   * gathered by nodes before the batched products or solves, and scattered
   * back, when the unknowns are ordered by nodes.
   */
  struct NodalBlockJacobiPreconditioner final : public mfem::Solver {
    //! \param[in] s: settings
    explicit NodalBlockJacobiPreconditioner(const NodalBlockJacobiSettings& s)
        : settings(s) {
      // the vectors ordered by nodes are gathered and scattered where MFEM
      // computes
      this->rhs.UseDevice(true);
      this->solution.UseDevice(true);
      this->columns.UseDevice(true);
    }
    void SetOperator(const mfem::Operator& op) override {
      if (const auto* const a = dynamic_cast<const mfem::SparseMatrix*>(&op);
          a != nullptr) {
        this->setup(*a);
        return;
      }
#ifdef MFEM_USE_MPI
      if (const auto* const a = dynamic_cast<const mfem::HypreParMatrix*>(&op);
          a != nullptr) {
        this->setup(*a);
        return;
      }
#endif /* MFEM_USE_MPI */
      if (const auto* const p = dynamic_cast<const NodalBlocksProvider*>(&op);
          p != nullptr) {
        this->setup(*p, op.Height());
        return;
      }
      mgis::raise(
          "NodalBlockJacobiPreconditioner::SetOperator: the operator is "
          "neither an assembled matrix nor an operator giving its nodal "
          "blocks");
    }
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override {
      if (this->blocks.TotalSize() == 0) {
        mgis::raise(
            "NodalBlockJacobiPreconditioner::Mult: the operator is not set");
      }
      if ((b.Size() != this->Width()) || (x.Size() != this->Height())) {
        mgis::raise("NodalBlockJacobiPreconditioner::Mult: invalid sizes");
      }
      const auto& batched = mfem::BatchedLinAlg::Get(this->settings.backend);
      const auto gathered =
          this->settings.ordered_by_nodes && (this->settings.block_size != 1);
      // vectors ordered by nodes
      if (gathered) {
        b.GetSubVector(this->nodes, this->rhs);
      }
      const auto& rb = gathered ? this->rhs : b;
      auto& xb = gathered ? this->solution : x;
      if (gathered) {
        xb.SetSize(b.Size());
      }
      if (this->settings.factorize) {
        xb = rb;
        this->synchronize();
        batched.LUSolve(this->blocks, this->pivots, xb);
      } else {
        this->synchronize();
        batched.Mult(this->blocks, rb, xb);
      }
      this->synchronize();
      if (gathered) {
        x.SetSubVector(this->nodes, xb);
      }
      if (this->settings.damping.has_value()) {
        x *= *(this->settings.damping);
      }
    }

   private:
    /*!
     * \brief extract the blocks of a matrix, and invert or factorize them
     * \param[in] a: matrix
     */
    void setup(const mfem::SparseMatrix& a) {
      const auto use_device = mfem::Device::IsEnabled();
      this->setup(a.Height(), a.Width(), a.ReadI(use_device),
                  a.ReadJ(use_device), a.ReadData(use_device), use_device);
    }
#ifdef MFEM_USE_MPI
    /*!
     * \brief extract the blocks of the diagonal block of the process of a
     * matrix of hypre, read where hypre stores it
     * \param[in] a: matrix
     */
    void setup(const mfem::HypreParMatrix& a) {
      static_assert(std::is_same_v<HYPRE_Int, int>,
                    "unsupported integer type of hypre");
      static_assert(std::is_same_v<HYPRE_Real, real>,
                    "unsupported real type of hypre");
      const auto use_device =
          mfem::Device::IsEnabled() &&
          (mfem::GetHypreMemoryClass() != mfem::MemoryClass::HOST);
      if (use_device) {
        a.HypreRead();
      } else {
        a.HostRead();
      }
      auto* const diagonal =
          hypre_ParCSRMatrixDiag(static_cast<hypre_ParCSRMatrix*>(a));
      this->setup(hypre_CSRMatrixNumRows(diagonal),
                  hypre_CSRMatrixNumCols(diagonal), hypre_CSRMatrixI(diagonal),
                  hypre_CSRMatrixJ(diagonal), hypre_CSRMatrixData(diagonal),
                  use_device);
      if (!use_device) {
        a.HypreRead();
      }
    }
#endif /* MFEM_USE_MPI */
    /*!
     * \brief extract the blocks of a matrix in the compressed row format
     * \param[in] n: number of rows
     * \param[in] m: number of columns
     * \param[in] I: offsets of the rows
     * \param[in] J: columns
     * \param[in] A: values
     * \param[in] use_device: if true, the arrays are read on the device
     */
    void setup(const int n,
               const int m,
               const int* const I,
               const int* const J,
               const real* const A,
               const bool use_device) {
      const auto bs = this->settings.block_size;
      if ((m != n) || (n % bs != 0)) {
        mgis::raise(
            "NodalBlockJacobiPreconditioner::SetOperator: the size of the "
            "matrix is not a multiple of the size of the blocks");
      }
      const auto l = this->getLayout(n);
      this->blocks.SetSize(bs, bs, l.nb);
      launchNodalBlocksExtractionKernel(
          l, I, J, A, this->blocks.Write(use_device), use_device);
      this->factorize(n);
    }
    /*!
     * \brief gather the blocks of an operator applied without assembly, and
     * invert or factorize them
     * \param[in] p: operator giving its nodal blocks
     * \param[in] n: number of unknowns
     */
    void setup(const NodalBlocksProvider& p, const int n) {
      const auto bs = this->settings.block_size;
      p.AssembleNodalBlocks(this->columns);
      if ((n % bs != 0) || (this->columns.Size() != bs * n)) {
        mgis::raise(
            "NodalBlockJacobiPreconditioner::SetOperator: the size of the "
            "blocks must be the number of components of the unknowns for an "
            "operator applied without assembly");
      }
      const auto l = this->getLayout(n);
      const auto use_device = mfem::Device::IsEnabled();
      this->blocks.SetSize(bs, bs, l.nb);
      launchNodalBlocksGatherKernel(l, this->columns.Read(use_device),
                                    this->blocks.Write(use_device), use_device);
      this->factorize(n);
    }
    /*!
     * \return the layout of the unknowns
     * \param[in] n: number of unknowns
     */
    [[nodiscard]] NodalBlocksLayout getLayout(const int n) const noexcept {
      const auto bs = this->settings.block_size;
      return NodalBlocksLayout{
          .nb = n / bs, .bs = bs, .by_nodes = this->settings.ordered_by_nodes};
    }
    /*!
     * \brief invert or factorize the blocks, and prepare the vectors
     * \param[in] n: number of unknowns
     */
    void factorize(const int n) {
      const auto bs = this->settings.block_size;
      if (this->settings.ordered_by_nodes && (bs != 1) &&
          (this->nodes.Size() != n)) {
        // position of the unknowns in the vectors ordered by nodes
        const auto nb = n / bs;
        this->nodes.SetSize(n);
        this->nodes.UseDevice(true);
        auto* const nodes_ = this->nodes.HostWrite();
        for (int i = 0; i != n; ++i) {
          nodes_[(i % nb) * bs + i / nb] = i;
        }
      }
      const auto& batched = mfem::BatchedLinAlg::Get(this->settings.backend);
      this->synchronize();
      if (this->settings.factorize) {
        batched.LUFactor(this->blocks, this->pivots);
      } else {
        batched.Invert(this->blocks);
      }
      this->synchronize();
      this->height = this->width = n;
    }
    /*!
     * \brief synchronize the device around the operations of the backends
     * cuBLAS, hipBLAS and MAGMA, which compute on their own queue: the
     * kernels of MFEM writing their arguments are completed before them,
     * and their results before the next kernels of MFEM. Without it, with
     * the backend MAGMA, the GMRES solver of MFEM in parallel on the GPU
     * gave arguments still being computed to MAGMA, and diverged.
     */
    void synchronize() const {
      if (this->settings.backend != mfem::BatchedLinAlg::NATIVE) {
        MFEM_DEVICE_SYNC;
      }
    }
    //! \brief settings
    const NodalBlockJacobiSettings settings;
    //! \brief inverses or LU factors of the blocks
    mfem::DenseTensor blocks;
    //! \brief pivots of the LU factors
    mfem::Array<int> pivots;
    //! \brief unknowns ordered by nodes
    mfem::Array<int> nodes;
    //! \brief columns of the blocks of an operator applied without assembly
    mfem::Vector columns;
    //! \brief right hand side ordered by nodes
    mutable mfem::Vector rhs;
    //! \brief solution ordered by nodes
    mutable mfem::Vector solution;
  };  // end of struct NodalBlockJacobiPreconditioner

  /*!
   * \return the backend of the batched linear algebra of MFEM given by its
   * name, or an empty value if it is unknown or not available
   * \param[in, out] ctx: execution context
   * \param[in] n: name
   */
  [[nodiscard]] static std::optional<mfem::BatchedLinAlg::Backend>
  getBatchedLinearAlgebraBackend(Context& ctx, const std::string& n) {
    using Batched = mfem::BatchedLinAlg;
    const auto backends =
        std::array<std::pair<const char*, Batched::Backend>, 3>{
            {{"Native", Batched::NATIVE},
             {"GPUBLAS", Batched::GPU_BLAS},
             {"MAGMA", Batched::MAGMA}}};
    for (const auto& [name, backend] : backends) {
      if (n != name) {
        continue;
      }
      if (!Batched::IsAvailable(backend)) {
        return ctx.registerErrorMessage(
            "the backend '" + n +
            "' of the batched linear algebra of MFEM is not available (it "
            "requires MFEM built with CUDA or HIP, and MAGMA for the backend "
            "'MAGMA', and the device of MFEM on the GPU)");
      }
      return backend;
    }
    return ctx.registerErrorMessage(
        "invalid backend '" + n +
        "' of the batched linear algebra of MFEM (expected 'Native', "
        "'GPUBLAS' or 'MAGMA')");
  }  // end of getBatchedLinearAlgebraBackend

  /*!
   * \return the Jacobi preconditioner by nodal blocks
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options (`Backend`: `Native`, `GPUBLAS` or `MAGMA`,
   * the active backend of MFEM by default; `BlockSize`, the number of
   * components of the unknowns by default; `Factorization`: `Inverse`, by
   * default, or `LU`; `Damping`, one by default)
   */
  template <bool parallel>
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  makeNodalBlockJacobiPreconditioner(Context& ctx,
                                     FiniteElementSpace<parallel>& fespace,
                                     const Parameters& opts) {
    using Batched = mfem::BatchedLinAlg;
    if (!checkParameters(
            ctx, opts,
            std::vector<std::string>{"Backend", "BlockSize", "Factorization",
                                     "Damping"})) {
      return {};
    }
    // the active backend of MFEM, if available (MAGMA, then cuBLAS or
    // hipBLAS, then the native kernels)
    const auto active = Batched::GetActiveBackend();
    auto settings = NodalBlockJacobiSettings{};
    settings.backend = Batched::IsAvailable(active) ? active : Batched::NATIVE;
    if (contains(opts, "Backend")) {
      const auto oname = get<std::string>(ctx, opts, "Backend");
      if (isInvalid(oname)) {
        return {};
      }
      const auto obackend = getBatchedLinearAlgebraBackend(ctx, *oname);
      if (isInvalid(obackend)) {
        return {};
      }
      settings.backend = *obackend;
    }
    const auto osize = get_if<int>(ctx, opts, "BlockSize", fespace.GetVDim());
    const auto ofactorization =
        get_if<std::string>(ctx, opts, "Factorization", std::string{"Inverse"});
    const auto odamping = get_if<double>(ctx, opts, "Damping", 1.);
    const auto damped = contains(opts, "Damping");
    if (isInvalid(osize) || isInvalid(ofactorization) || isInvalid(odamping)) {
      return {};
    }
    if ((*ofactorization != "Inverse") && (*ofactorization != "LU")) {
      return ctx.registerErrorMessage(
          "invalid factorization '" + *ofactorization +
          "' of the Jacobi preconditioner by nodal blocks (expected "
          "'Inverse' or 'LU')");
    }
    if (*osize < 1) {
      return ctx.registerErrorMessage(
          "invalid block size of the Jacobi preconditioner by nodal blocks");
    }
    if (*odamping <= 0) {
      return ctx.registerErrorMessage(
          "invalid damping factor of the Jacobi preconditioner by nodal "
          "blocks");
    }
    settings.block_size = *osize;
    settings.ordered_by_nodes =
        fespace.GetOrdering() == mfem::Ordering::byNODES;
    settings.factorize = *ofactorization == "LU";
    if (damped) {
      settings.damping = static_cast<real>(*odamping);
    }
    return std::make_unique<NodalBlockJacobiPreconditioner>(settings);
  }  // end of makeNodalBlockJacobiPreconditioner

  template <bool parallel>
  void declareNodalBlockJacobiPreconditioner(
      PreconditionerFactory<parallel>& f) {
    declare(f, "NodalBlockJacobi",
            LinearSolverDescription{
                .library = "MFEM",
                .kind = LinearSolverKind::PRECONDITIONER,
                .placement = LinearSolverPlacement::MFEM,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Backend", "BlockSize", "Factorization", "Damping"},
                .description =
                    "Jacobi preconditioner by nodal blocks, inverted by the "
                    "batched dense linear algebra of MFEM, applied to the "
                    "diagonal block of each process in parallel"},
            makeNodalBlockJacobiPreconditioner<parallel>);
  }  // end of declareNodalBlockJacobiPreconditioner

#ifdef MFEM_USE_MPI
  template void declareNodalBlockJacobiPreconditioner<true>(
      PreconditionerFactory<true>&);
#endif /* MFEM_USE_MPI */
  template void declareNodalBlockJacobiPreconditioner<false>(
      PreconditionerFactory<false>&);

}  // end of namespace mfem_mgis::internals
