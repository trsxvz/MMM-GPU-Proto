/*!
 * \file   src/PETScSupport.cxx
 * \brief  This file implements the functions and the classes supporting the
 * solvers and the preconditioners of PETSc
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC)

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
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
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

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC) */
