/*!
 * \file   src/PETScSupport.cxx
 * \brief  This file implements the functions and the classes supporting the
 * solvers and the preconditioners of PETSc
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC)

#include <array>
#include <tuple>
#include <atomic>
#include <string>
#include <iomanip>
#include <sstream>
#include <optional>
#include <utility>
#include <algorithm>
#include "petscsys.h"
#include "petscvec.h"
#include "petscmat.h"
#include "mfem/fem/pfespace.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
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

  bool PETScOptions::setOptions(Context& ctx, const Parameters& options) {
    for (const auto& [n, p] : options) {
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

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC) */
