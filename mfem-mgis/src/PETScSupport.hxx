/*!
 * \file   src/PETScSupport.hxx
 * \brief  This file declares the functions and the classes supporting the
 * solvers and the preconditioners of PETSc (internal header)
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_PETSCSUPPORT_HXX
#define LIB_MFEM_MGIS_SRC_PETSCSUPPORT_HXX

#include "mfem/config/config.hpp"
#if defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC)
#include <string>
#include <vector>
#include <string_view>
#include <mpi.h>
#include "petscvec.h"
#include "petscmat.h"
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Parameters.hxx"
#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC) */
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief initialize PETSc, if not done yet by mfem-mgis, for its solvers
   * and its preconditioners, PETSc being finalized at the end of the
   * execution of mfem-mgis. This does not select the non linear solver of
   * PETSc (see `usePETSc`). This function shall not be called if the program
   * initialized PETSc itself.
   */
  void initializePETScOnDemand();

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC)

  /*!
   * \brief initialize PETSc, unless it is already initialized (by
   * mfem-mgis, or by the program, which then finalizes it)
   */
  void initializePETSc();

  /*!
   * \brief options of a solver of PETSc, written in the options database of
   * PETSc with a prefix which is unique to the solver, and removed from it
   * when the solver is destroyed. This class is the first base class (or
   * the first member) of the solvers of PETSc, so that the prefix is known
   * when the object of PETSc is created and that the options are removed
   * after its destruction.
   *
   * PETSc ignores the options it does not know (for instance a misspelt
   * option): the options which PETSc never read are reported when a solver
   * which has been used is destroyed.
   */
  struct PETScOptions {
    /*!
     * \param[in] c: communicator
     * \param[in] kind: kind of object (`ksp`, `pc`, `snes`)
     */
    PETScOptions(MPI_Comm c, const std::string_view kind);
    PETScOptions(PETScOptions&&) = delete;
    PETScOptions(const PETScOptions&) = delete;
    PETScOptions& operator=(PETScOptions&&) = delete;
    PETScOptions& operator=(const PETScOptions&) = delete;
    /*!
     * \brief write options in the options database of PETSc
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] options: options, given without the leading dash and the
     * prefix (for instance `{"ksp_type", "cg"}`)
     */
    [[nodiscard]] bool setOptions(Context& ctx, const Parameters& options);
    /*!
     * \brief destructor, reporting the options which PETSc never read, if
     * the solver has been used, and removing the options from the options
     * database (unless PETSc has been finalized)
     */
    ~PETScOptions();
    //! \brief mark the solver as used (called by the `Mult` methods)
    void markAsUsed() const noexcept;
    //! \brief prefix of the options
    const std::string prefix;

   private:
    //! \brief communicator of the solver
    const MPI_Comm communicator;
    //! \brief options written in the options database, with their prefix
    std::vector<std::string> options;
    //! \brief if true, the solver has been used
    mutable bool used = false;
  };  // end of struct PETScOptions

  //! \brief vector of PETSc, destroyed with this object
  struct PETScVector {
    PETScVector() = default;
    PETScVector(PETScVector&&) = delete;
    PETScVector(const PETScVector&) = delete;
    PETScVector& operator=(PETScVector&&) = delete;
    PETScVector& operator=(const PETScVector&) = delete;
    //! \brief destructor
    ~PETScVector();
    //! \brief vector
    Vec v = nullptr;
  };  // end of struct PETScVector

  /*!
   * \brief copy a vector of MFEM in a vector of PETSc, through the array of
   * the host of PETSc (PETSc then synchronizes its arrays, whatever the type
   * of the vector)
   * \param[out] dst: vector of PETSc
   * \param[in] src: vector of MFEM
   */
  void copyToPETScVector(Vec dst, const mfem::Vector& src);
  /*!
   * \brief copy a vector of PETSc in a vector of MFEM, through the array of
   * the host of PETSc
   * \param[out] dst: vector of MFEM
   * \param[in] src: vector of PETSc
   */
  void copyFromPETScVector(mfem::Vector& dst, Vec src);
  /*!
   * \brief create a vector of PETSc of the host
   * \param[out] v: vector
   * \param[in] c: communicator
   * \param[in] n: local size
   */
  void createHostVector(PETScVector& v, MPI_Comm c, const int n);
  /*!
   * \brief convert a matrix of PETSc of the GPU (AIJ formats of CUDA or HIP)
   * to the AIJ format of the host, in place
   * \param[in, out] a: matrix
   */
  void convertToHostMatrix(Mat a);

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC) */

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_PETSCSUPPORT_HXX */
