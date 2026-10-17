/*!
 * \file   src/PETScSupport.cxx
 * \brief  This file implements the functions and the classes supporting the
 * solvers and the preconditioners of PETSc
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC)

#include <array>
#include <vector>
#include <memory>
#include <tuple>
#include <atomic>
#include <string>
#include <iomanip>
#include <sstream>
#include <optional>
#include <utility>
#include <algorithm>
#include <type_traits>
#include "petscsys.h"
#include "petscvec.h"
#include "petscdevice.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscis.h"
#include "mfem/general/device.hpp"
#include "mfem/general/backends.hpp"
#include "mfem/general/mem_manager.hpp"
#include "mfem/fem/pfespace.hpp"
#include "mfem/linalg/hypre.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/ExecutionSettings.hxx"
#include "ExecutionSettings.hxx"
#include "NodesCoordinates.hxx"
#include "PETScSupport.hxx"

namespace mfem_mgis::internals {

  void initializePETSc() {
    auto initialized = PetscBool{};
    if ((PetscInitialized(&initialized) == 0) && (initialized == PETSC_TRUE)) {
      return;
    }
    initializePETScOnDemand();
  }  // end of initializePETSc

  /*!
   * \return a prefix of the options database of PETSc which is unique in
   * the process, so that the options of two solvers do not interfere
   * \param[in] kind: kind of object (`ksp`, `pc`, `snes`)
   */
  [[nodiscard]] static std::string makePETScOptionsPrefix(
      const std::string_view kind) {
    static auto count = std::atomic<std::size_t>{};
    return "mfem_mgis_" + std::string{kind} + std::to_string(count++) + "_";
  }  // end of makePETScOptionsPrefix

  /*!
   * \return the value of an option of PETSc given by a parameter, or an
   * empty value on error: a boolean (`true` or `false`), an integer, a
   * floating-point value (written with all its digits) or a string (an empty
   * string being an option without value)
   * \param[in, out] ctx: execution context
   * \param[in] n: name of the option
   * \param[in] p: parameter
   */
  [[nodiscard]] static std::optional<std::string> getPETScOptionValue(
      Context& ctx, const std::string& n, const Parameter& p) {
    if (mfem_mgis::is<bool>(p)) {
      const auto ov = mfem_mgis::get<bool>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      return std::string{*ov ? "true" : "false"};
    }
    if (mfem_mgis::is<int>(p)) {
      const auto ov = mfem_mgis::get<int>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      return std::to_string(*ov);
    }
    if (mfem_mgis::is<double>(p)) {
      const auto ov = mfem_mgis::get<double>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      auto os = std::ostringstream{};
      os << std::setprecision(17) << *ov;
      return os.str();
    }
    if (mfem_mgis::is<std::string>(p)) {
      const auto ov = mfem_mgis::get<std::string>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      return std::string{*ov};
    }
    std::ignore = ctx.registerErrorMessage(
        "invalid value of the option '" + n +
        "' of PETSc (a boolean, an integer, a floating-point value or a "
        "string is expected)");
    return {};
  }  // end of getPETScOptionValue

  PETScOptions::PETScOptions(MPI_Comm c, const std::string_view kind)
      : prefix(makePETScOptionsPrefix(kind)), communicator(c) {}

  bool PETScOptions::setOptions(Context& ctx, const Parameters& opts) {
    for (const auto& [n, p] : opts) {
      if (n.empty() || n.starts_with('-')) {
        return ctx.registerErrorMessage(
            "invalid option '" + n +
            "' of PETSc (the options are given without the leading dash)");
      }
      const auto ovalue = getPETScOptionValue(ctx, n, p);
      if (isInvalid(ovalue)) {
        return false;
      }
      auto option = "-" + this->prefix + n;
      const auto* const value = ovalue->empty() ? nullptr : ovalue->c_str();
      if (PetscOptionsSetValue(nullptr, option.c_str(), value) != 0) {
        return ctx.registerErrorMessage("the option '" + n +
                                        "' can't be given to PETSc");
      }
      this->options.push_back(std::move(option));
    }
    return true;
  }  // end of setOptions

  void PETScOptions::markAsUsed() const noexcept {
    this->used = true;
  }  // end of markAsUsed

  void PETScOptions::finalizePETScAtExitOnce() const {
    if (!this->finalized_at_exit) {
      finalizePETScAtExit();
      this->finalized_at_exit = true;
    }
  }  // end of finalizePETScAtExitOnce

  PETScOptions::~PETScOptions() {
    auto finalized = PetscBool{};
    if ((PetscFinalized(&finalized) != 0) || (finalized == PETSC_TRUE)) {
      return;
    }
    auto rank = int{};
    MPI_Comm_rank(this->communicator, &rank);
    for (const auto& o : this->options) {
      auto read = PetscBool{};
      if ((this->used) && (rank == 0) &&
          (PetscOptionsUsed(nullptr, o.c_str() + 1, &read) == 0) &&
          (read == PETSC_FALSE)) {
        getErrorStream() << "warning: the option '"
                         << o.substr(this->prefix.size() + 1)
                         << "' was never read by PETSc\n";
      }
      std::ignore = PetscOptionsClearValue(nullptr, o.c_str());
    }
  }  // end of ~PETScOptions

  PETScVector::~PETScVector() {
    if (this->v != nullptr) {
      std::ignore = VecDestroy(&(this->v));
    }
  }  // end of ~PETScVector

  void copyToPETScVector(Vec dst, const mfem::Vector& src) {
    auto* a = static_cast<PetscScalar*>(nullptr);
    if (VecGetArrayWrite(dst, &a) != 0) {
      mgis::raise("copyToPETScVector: VecGetArrayWrite failed");
    }
    const auto* const s = src.HostRead();
    std::copy(s, s + src.Size(), a);
    if (VecRestoreArrayWrite(dst, &a) != 0) {
      mgis::raise("copyToPETScVector: VecRestoreArrayWrite failed");
    }
  }  // end of copyToPETScVector

  void copyFromPETScVector(mfem::Vector& dst, Vec src) {
    const auto* a = static_cast<const PetscScalar*>(nullptr);
    if (VecGetArrayRead(src, &a) != 0) {
      mgis::raise("copyFromPETScVector: VecGetArrayRead failed");
    }
    std::copy(a, a + dst.Size(), dst.HostWrite());
    if (VecRestoreArrayRead(src, &a) != 0) {
      mgis::raise("copyFromPETScVector: VecRestoreArrayRead failed");
    }
  }  // end of copyFromPETScVector

  void createHostVector(PETScVector& v, MPI_Comm c, const int n) {
    if ((VecCreate(c, &(v.v)) != 0) ||
        (VecSetSizes(v.v, static_cast<PetscInt>(n), PETSC_DECIDE) != 0) ||
        (VecSetType(v.v, VECSTANDARD) != 0)) {
      mgis::raise("createHostVector: the vector can't be created");
    }
  }  // end of createHostVector

  void convertToHostMatrix(Mat a) {
    if (a == nullptr) {
      return;
    }
    auto device = PetscBool{};
    if (PetscObjectTypeCompareAny(
            reinterpret_cast<PetscObject>(a), &device, MATSEQAIJCUSPARSE,
            MATMPIAIJCUSPARSE, MATAIJCUSPARSE, MATSEQAIJHIPSPARSE,
            MATMPIAIJHIPSPARSE, MATAIJHIPSPARSE, "") != 0) {
      mgis::raise("convertToHostMatrix: PetscObjectTypeCompareAny failed");
    }
    if (device == PETSC_FALSE) {
      return;
    }
    if (MatConvert(a, MATAIJ, MAT_INPLACE_MATRIX, &a) != 0) {
      mgis::raise("convertToHostMatrix: MatConvert failed");
    }
  }  // end of convertToHostMatrix

  /*!
   * \return a matrix of PETSc in the AIJ format of the host, copy of a
   * matrix of hypre, on its communicator and with its distribution
   * \param[in] a: matrix of hypre
   */
  static Mat copyHypreMatrix(const mfem::HypreParMatrix& a) {
    // the arrays of hypre are read on the host, and given back to hypre,
    // whose computations may be on the GPU, even on error
    struct HypreLocationGuard {
      ~HypreLocationGuard() { this->matrix.HypreRead(); }
      const mfem::HypreParMatrix& matrix;
    };
    a.HostRead();
    const auto guard = HypreLocationGuard{a};
    auto* const h =
        static_cast<hypre_ParCSRMatrix*>(const_cast<mfem::HypreParMatrix&>(a));
    auto* const diag = hypre_ParCSRMatrixDiag(h);
    auto* const offd = hypre_ParCSRMatrixOffd(h);
    const auto* const cmap = hypre_ParCSRMatrixColMapOffd(h);
    const auto* const di = hypre_CSRMatrixI(diag);
    const auto* const dj = hypre_CSRMatrixJ(diag);
    const auto* const dv = hypre_CSRMatrixData(diag);
    // the off-diagonal block may have no row offsets when it is empty
    const auto* const oi = hypre_CSRMatrixI(offd);
    const auto* const oj = hypre_CSRMatrixJ(offd);
    const auto* const ov = hypre_CSRMatrixData(offd);
    const auto nrows = static_cast<PetscInt>(hypre_CSRMatrixNumRows(diag));
    const auto ncolumns = static_cast<PetscInt>(hypre_CSRMatrixNumCols(diag));
    const auto first_row =
        static_cast<PetscInt>(hypre_ParCSRMatrixFirstRowIndex(h));
    const auto first_column =
        static_cast<PetscInt>(hypre_ParCSRMatrixFirstColDiag(h));
    auto offd_size = [oi](const PetscInt r) -> PetscInt {
      return oi == nullptr ? 0 : static_cast<PetscInt>(oi[r + 1] - oi[r]);
    };
    auto d_nnz = std::vector<PetscInt>(static_cast<std::size_t>(nrows));
    auto o_nnz = std::vector<PetscInt>(static_cast<std::size_t>(nrows));
    for (PetscInt r = 0; r != nrows; ++r) {
      d_nnz[r] = static_cast<PetscInt>(di[r + 1] - di[r]);
      o_nnz[r] = offd_size(r);
    }
    auto m = Mat{};
    auto check = [&m](const PetscErrorCode e) {
      if (e != 0) {
        std::ignore = MatDestroy(&m);
        mgis::raise("makeHostMatrix: the matrix of hypre can't be copied");
      }
    };
    check(MatCreate(a.GetComm(), &m));
    check(MatSetSizes(m, nrows, ncolumns, PETSC_DETERMINE, PETSC_DETERMINE));
    check(MatSetType(m, MATAIJ));
    // only the preallocation of the actual type is used
    check(MatSeqAIJSetPreallocation(m, 0, d_nnz.data()));
    check(MatMPIAIJSetPreallocation(m, 0, d_nnz.data(), 0, o_nnz.data()));
    // each process sets its own rows only (an entry of another process is
    // an error)
    check(MatSetOption(m, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE));
    auto columns = std::vector<PetscInt>{};
    auto values = std::vector<PetscScalar>{};
    for (PetscInt r = 0; r != nrows; ++r) {
      columns.clear();
      values.clear();
      for (auto k = di[r]; k != di[r + 1]; ++k) {
        columns.push_back(first_column + static_cast<PetscInt>(dj[k]));
        values.push_back(static_cast<PetscScalar>(dv[k]));
      }
      if (offd_size(r) != 0) {
        for (auto k = oi[r]; k != oi[r + 1]; ++k) {
          columns.push_back(static_cast<PetscInt>(cmap[oj[k]]));
          values.push_back(static_cast<PetscScalar>(ov[k]));
        }
      }
      const auto row = first_row + r;
      check(MatSetValues(m, 1, &row, static_cast<PetscInt>(columns.size()),
                         columns.data(), values.data(), INSERT_VALUES));
    }
    check(MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY));
    check(MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY));
    return m;
  }  // end of copyHypreMatrix

  std::unique_ptr<mfem::PetscParMatrix> makeHostMatrix(
      MPI_Comm c, const mfem::Operator& op) {
    const auto* const a = dynamic_cast<const mfem::HypreParMatrix*>(&op);
    if (a != nullptr) {
      return std::make_unique<mfem::PetscParMatrix>(copyHypreMatrix(*a), false);
    }
    // a matrix of PETSc is referenced, and converted in place to the format
    // of the host if needed
    const auto* const p = dynamic_cast<const mfem::PetscParMatrix*>(&op);
    if (p != nullptr) {
      auto m =
          std::make_unique<mfem::PetscParMatrix>(static_cast<Mat>(*p), true);
      convertToHostMatrix(*m);
      return m;
    }
    auto m = std::make_unique<mfem::PetscParMatrix>(
        c, &op, mfem::Operator::PETSC_MATAIJ);
    convertToHostMatrix(*m);
    return m;
  }  // end of makeHostMatrix

  bool usePETScOnDevice() noexcept {
#if defined(MFEM_USE_CUDA) && defined(PETSC_HAVE_CUDA)
    const auto device = mfem::Device::Allows(mfem::Backend::CUDA_MASK);
#elif defined(MFEM_USE_HIP) && defined(PETSC_HAVE_HIP)
    const auto device = mfem::Device::Allows(mfem::Backend::HIP_MASK);
#else
    const auto device = false;
#endif
    if ((!device) || (!isExecutionModeSpecified())) {
      return device;
    }
    return getExecutionSettings().linear_solver == ExecutionLocation::GPU;
  }  // end of usePETScOnDevice

  /*!
   * \brief convert a matrix of PETSc in the AIJ format of the host to the AIJ
   * format of the GPU (CUDA or HIP), in place
   * \param[in, out] a: matrix
   */
  static void convertToDeviceMatrix(Mat a) {
    if (a == nullptr) {
      return;
    }
    auto host = PetscBool{};
    if (PetscObjectTypeCompareAny(reinterpret_cast<PetscObject>(a), &host,
                                  MATSEQAIJ, MATMPIAIJ, "") != 0) {
      mgis::raise("convertToDeviceMatrix: PetscObjectTypeCompareAny failed");
    }
    if (host == PETSC_FALSE) {
      return;
    }
#if defined(MFEM_USE_CUDA) && defined(PETSC_HAVE_CUDA)
    if (MatConvert(a, MATAIJCUSPARSE, MAT_INPLACE_MATRIX, &a) != 0) {
      mgis::raise("convertToDeviceMatrix: MatConvert failed");
    }
#elif defined(MFEM_USE_HIP) && defined(PETSC_HAVE_HIP)
    if (MatConvert(a, MATAIJHIPSPARSE, MAT_INPLACE_MATRIX, &a) != 0) {
      mgis::raise("convertToDeviceMatrix: MatConvert failed");
    }
#else
    mgis::raise("convertToDeviceMatrix: PETSc is not built for the GPU");
#endif
  }  // end of convertToDeviceMatrix

  std::unique_ptr<mfem::PetscParMatrix> makeDeviceMatrix(
      MPI_Comm c, const mfem::Operator& op) {
    const auto* const a = dynamic_cast<const mfem::HypreParMatrix*>(&op);
    if (a != nullptr) {
      // the arrays of hypre are where hypre computes
      a->HypreRead();
      auto m =
          std::make_unique<mfem::PetscParMatrix>(a, mfem::Operator::PETSC_MATAIJ);
      convertToDeviceMatrix(*m);
      return m;
    }
    const auto* const p = dynamic_cast<const mfem::PetscParMatrix*>(&op);
    if (p != nullptr) {
      auto m =
          std::make_unique<mfem::PetscParMatrix>(static_cast<Mat>(*p), true);
      convertToDeviceMatrix(*m);
      return m;
    }
    auto m = std::make_unique<mfem::PetscParMatrix>(
        c, &op, mfem::Operator::PETSC_MATAIJ);
    convertToDeviceMatrix(*m);
    return m;
  }  // end of makeDeviceMatrix

  void createDeviceVector(PETScVector& v,
                          MPI_Comm c,
                          const real* const a,
                          const int n) {
    static_assert(std::is_same_v<real, PetscScalar>,
                  "the vectors of MFEM and of PETSc shall have the same "
                  "numerical type");
    auto size = int{};
    if (MPI_Comm_size(c, &size) != MPI_SUCCESS) {
      mgis::raise("createDeviceVector: MPI_Comm_size failed");
    }
    const auto ln = static_cast<PetscInt>(n);
#if defined(MFEM_USE_CUDA) && defined(PETSC_HAVE_CUDA)
    const auto e = (size == 1)
                       ? VecCreateSeqCUDAWithArray(c, 1, ln, a, &(v.v))
                       : VecCreateMPICUDAWithArray(c, 1, ln, PETSC_DECIDE, a,
                                                   &(v.v));
#elif defined(MFEM_USE_HIP) && defined(PETSC_HAVE_HIP)
    const auto e = (size == 1)
                       ? VecCreateSeqHIPWithArray(c, 1, ln, a, &(v.v))
                       : VecCreateMPIHIPWithArray(c, 1, ln, PETSC_DECIDE, a,
                                                  &(v.v));
#else
    std::ignore = ln;
    std::ignore = a;
    const auto e = PETSC_ERR_SUP;
#endif
    if (e != 0) {
      mgis::raise("createDeviceVector: the vector can't be created");
    }
  }  // end of createDeviceVector

  void synchronizeDeviceVector(Vec v) {
    // the array of the GPU is given by PETSc after the copy of the values of
    // the host, if needed
    const auto* a = static_cast<const PetscScalar*>(nullptr);
    auto t = PetscMemType{};
    if ((VecGetArrayReadAndMemType(v, &a, &t) != 0) ||
        (VecRestoreArrayReadAndMemType(v, &a) != 0)) {
      mgis::raise("synchronizeDeviceVector: the array can't be read");
    }
  }  // end of synchronizeDeviceVector

  /*!
   * \brief copy an array of the GPU in another one, after the end of the
   * computations of PETSc and of MFEM on the GPU, the copy being finished on
   * return (PETSc computes in its own streams; the copies of PETSc only
   * accept the arrays allocated by PETSc)
   * \param[out] dst: destination
   * \param[in] src: source
   * \param[in] n: size
   */
  static void copyDeviceArray(PetscScalar* const dst,
                              const PetscScalar* const src,
                              const PetscInt n) {
    if (n == 0) {
      return;
    }
    auto dctx = PetscDeviceContext{};
    if ((PetscDeviceContextGetCurrentContext(&dctx) != 0) ||
        (PetscDeviceContextSynchronize(dctx) != 0)) {
      mgis::raise("copyDeviceArray: the context of PETSc can't be synchronized");
    }
    const auto bytes = static_cast<std::size_t>(n) * sizeof(PetscScalar);
#if defined(MFEM_USE_CUDA)
    mfem::CuMemcpyDtoD(dst, src, bytes);
#elif defined(MFEM_USE_HIP)
    mfem::HipMemcpyDtoD(dst, src, bytes);
#else
    std::ignore = bytes;
    mgis::raise("copyDeviceArray: MFEM is built without GPU support");
#endif
    MFEM_DEVICE_SYNC;
  }  // end of copyDeviceArray

  void applyOnDevice(const mfem::Operator& op,
                     mfem::Vector& r,
                     mfem::Vector& z,
                     Vec x,
                     Vec y) {
    // arrays of the GPU of the vectors of PETSc, given back to PETSc on exit
    // (the values of the array of x being up to date on the GPU, y being
    // written on the GPU)
    struct DeviceArrays {
      ~DeviceArrays() {
        if (this->ya != nullptr) {
          std::ignore = VecRestoreArrayWriteAndMemType(this->y, &(this->ya));
        }
        if (this->xa != nullptr) {
          std::ignore = VecRestoreArrayReadAndMemType(this->x, &(this->xa));
        }
      }
      Vec x;
      Vec y;
      const PetscScalar* xa = nullptr;
      PetscScalar* ya = nullptr;
    };
    auto nx = PetscInt{};
    auto ny = PetscInt{};
    if ((VecGetLocalSize(x, &nx) != 0) || (VecGetLocalSize(y, &ny) != 0)) {
      mgis::raise("applyOnDevice: VecGetLocalSize failed");
    }
    auto arrays = DeviceArrays{.x = x, .y = y};
    auto tx = PetscMemType{};
    auto ty = PetscMemType{};
    if ((VecGetArrayReadAndMemType(x, &(arrays.xa), &tx) != 0) ||
        (VecGetArrayWriteAndMemType(y, &(arrays.ya), &ty) != 0)) {
      mgis::raise("applyOnDevice: the arrays of the vectors can't be read");
    }
    if (((nx != 0) && (!PetscMemTypeDevice(tx))) ||
        ((ny != 0) && (!PetscMemTypeDevice(ty)))) {
      mgis::raise("applyOnDevice: the vectors are not vectors of the GPU");
    }
    r.UseDevice(true);
    r.SetSize(static_cast<int>(nx));
    z.UseDevice(true);
    z.SetSize(static_cast<int>(ny));
    copyDeviceArray(r.Write(), arrays.xa, nx);
    op.Mult(r, z);
    // `Read` copies the values computed on the host by the operator, if any,
    // on the GPU
    copyDeviceArray(arrays.ya, z.Read(), ny);
  }  // end of applyOnDevice

  void PETScRigidBodyModes::initialize(FiniteElementSpace<true>& fespace) {
    this->dimension = 0;
    this->coordinates.clear();
    if ((!hasRigidBodyModes(fespace)) ||
        (fespace.GetOrdering() != mfem::Ordering::byVDIM)) {
      return;
    }
    const auto c = getNodesCoordinates(fespace);
    const auto d = static_cast<std::size_t>(c.dimension);
    const auto n = c.values.size() / d;
    // coordinates relative to the centroid of the nodes of all the
    // processes: PETSc orthogonalizes the rotations against the
    // translations in a single pass, which loses accuracy far from the
    // origin
    auto sums = std::array<double, 4>{0, 0, 0, static_cast<double>(n)};
    for (std::size_t k = 0; k != d; ++k) {
      for (std::size_t node = 0; node != n; ++node) {
        sums[k] += static_cast<double>(c.values[k * n + node]);
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, sums.data(), 4, MPI_DOUBLE, MPI_SUM,
                  fespace.GetComm());
    this->coordinates.resize(c.values.size());
    for (std::size_t node = 0; node != n; ++node) {
      for (std::size_t k = 0; k != d; ++k) {
        this->coordinates[node * d + k] = static_cast<PetscScalar>(
            static_cast<double>(c.values[k * n + node]) -
            sums[k] / std::max(sums[3], 1.));
      }
    }
    this->dimension = c.dimension;
  }  // end of initialize

  void PETScRigidBodyModes::attach(Mat m) const {
    if ((this->dimension == 0) || (m == nullptr)) {
      return;
    }
    // only the assembled matrices (not the shell matrices of the operators
    // applied without assembly)
    auto aij = PetscBool{};
    if (PetscObjectBaseTypeCompareAny(reinterpret_cast<PetscObject>(m), &aij,
                                      MATSEQAIJ, MATMPIAIJ, "") != 0) {
      mgis::raise(
          "PETScRigidBodyModes::attach: PetscObjectBaseTypeCompareAny failed");
    }
    if (aij == PETSC_FALSE) {
      return;
    }
    auto comm = MPI_Comm{};
    auto nrows = PetscInt{};
    auto ncolumns = PetscInt{};
    if ((PetscObjectGetComm(reinterpret_cast<PetscObject>(m), &comm) != 0) ||
        (MatGetLocalSize(m, &nrows, &ncolumns) != 0)) {
      mgis::raise("PETScRigidBodyModes::attach: the matrix can't be described");
    }
    // the modes are only given to the matrices of the finite element space
    // (a matrix of another size or distribution is left unchanged), the
    // decision being taken on all the processes
    auto mismatch = static_cast<int>(static_cast<std::size_t>(nrows) !=
                                     this->coordinates.size());
    MPI_Allreduce(MPI_IN_PLACE, &mismatch, 1, MPI_INT, MPI_LOR, comm);
    if (mismatch != 0) {
      return;
    }
    auto x = PETScVector{};
    const auto n = static_cast<PetscInt>(this->coordinates.size());
    auto* values = static_cast<PetscScalar*>(nullptr);
    if ((MatSetBlockSize(m, this->dimension) != 0) ||
        (VecCreate(comm, &(x.v)) != 0) ||
        (VecSetSizes(x.v, n, PETSC_DETERMINE) != 0) ||
        (VecSetBlockSize(x.v, this->dimension) != 0) ||
        (VecSetType(x.v, VECSTANDARD) != 0) ||
        (VecGetArray(x.v, &values) != 0)) {
      mgis::raise(
          "PETScRigidBodyModes::attach: the coordinates can't be given to "
          "PETSc");
    }
    std::copy(this->coordinates.begin(), this->coordinates.end(), values);
    auto nsp = MatNullSpace{};
    if ((VecRestoreArray(x.v, &values) != 0) ||
        (MatNullSpaceCreateRigidBody(x.v, &nsp) != 0)) {
      mgis::raise(
          "PETScRigidBodyModes::attach: MatNullSpaceCreateRigidBody failed");
    }
    const auto status = MatSetNearNullSpace(m, nsp);
    std::ignore = MatNullSpaceDestroy(&nsp);
    if (status != 0) {
      mgis::raise("PETScRigidBodyModes::attach: MatSetNearNullSpace failed");
    }
  }  // end of attach

  void PETScComponentFields::initialize(FiniteElementSpace<true>& fespace) {
    this->rows.clear();
    this->size = static_cast<PetscInt>(fespace.GetTrueVSize());
    const auto vdim = fespace.GetVDim();
    if (vdim < 2) {
      return;
    }
    const auto ndofs = fespace.GetNDofs();
    const auto by_nodes = fespace.GetOrdering() == mfem::Ordering::byNODES;
    const auto offset = static_cast<PetscInt>(fespace.GetMyTDofOffset());
    // component of each true degree of freedom of the process, given by its
    // local degree of freedom (the true degrees of freedom are a subset of
    // the local ones)
    auto components =
        std::vector<int>(static_cast<std::size_t>(this->size), -1);
    for (int ldof = 0; ldof != fespace.GetVSize(); ++ldof) {
      const auto tdof = fespace.GetLocalTDofNumber(ldof);
      if (tdof < 0) {
        continue;
      }
      components[static_cast<std::size_t>(tdof)] =
          by_nodes ? ldof / ndofs : ldof % vdim;
    }
    this->rows.resize(static_cast<std::size_t>(vdim));
    for (std::size_t tdof = 0; tdof != components.size(); ++tdof) {
      const auto c = components[tdof];
      if ((c < 0) || (c >= vdim)) {
        mgis::raise(
            "PETScComponentFields::initialize: a true degree of freedom has "
            "no component");
      }
      this->rows[static_cast<std::size_t>(c)].push_back(
          offset + static_cast<PetscInt>(tdof));
    }
  }  // end of initialize

  void PETScComponentFields::attach(PC pc) const {
    if ((this->rows.empty()) || (pc == nullptr)) {
      return;
    }
    auto fieldsplit = PetscBool{};
    if (PetscObjectTypeCompare(reinterpret_cast<PetscObject>(pc), PCFIELDSPLIT,
                               &fieldsplit) != 0) {
      mgis::raise(
          "PETScComponentFields::attach: PetscObjectTypeCompare failed");
    }
    if (fieldsplit == PETSC_FALSE) {
      return;
    }
    // fields already given (by a previous call)
    auto first = IS{};
    if (PCFieldSplitGetIS(pc, "0", &first) != 0) {
      mgis::raise("PETScComponentFields::attach: PCFieldSplitGetIS failed");
    }
    if (first != nullptr) {
      return;
    }
    auto a = Mat{};
    auto p = Mat{};
    auto comm = MPI_Comm{};
    auto nrows = PetscInt{};
    auto ncolumns = PetscInt{};
    if ((PCGetOperators(pc, &a, &p) != 0) ||
        (PetscObjectGetComm(reinterpret_cast<PetscObject>(pc), &comm) != 0) ||
        (MatGetLocalSize(p, &nrows, &ncolumns) != 0)) {
      mgis::raise(
          "PETScComponentFields::attach: the operator of the preconditioner "
          "can't be described");
    }
    // the fields are only given to the operators of the finite element
    // space, the decision being taken on all the processes
    auto mismatch = static_cast<int>(nrows != this->size);
    MPI_Allreduce(MPI_IN_PLACE, &mismatch, 1, MPI_INT, MPI_LOR, comm);
    if (mismatch != 0) {
      return;
    }
    for (std::size_t c = 0; c != this->rows.size(); ++c) {
      const auto& r = this->rows[c];
      auto is = IS{};
      if (ISCreateGeneral(comm, static_cast<PetscInt>(r.size()), r.data(),
                          PETSC_COPY_VALUES, &is) != 0) {
        mgis::raise("PETScComponentFields::attach: ISCreateGeneral failed");
      }
      // the preconditioner keeps its own reference to the set
      const auto name = std::to_string(c);
      const auto status = PCFieldSplitSetIS(pc, name.c_str(), is);
      ISDestroy(&is);
      if (status != 0) {
        mgis::raise("PETScComponentFields::attach: PCFieldSplitSetIS failed");
      }
    }
  }  // end of attach

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC) */
