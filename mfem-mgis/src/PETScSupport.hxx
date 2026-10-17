/*!
 * \file   src/PETScSupport.hxx
 * \brief  This file declares the functions and the classes supporting the
 * solvers and the preconditioners of PETSc (internal header)
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_PETSCSUPPORT_HXX
#define LIB_MFEM_MGIS_SRC_PETSCSUPPORT_HXX

#include "mfem/config/config.hpp"
#if defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC)
#include <memory>
#include <string>
#include <vector>
#include <string_view>
#include <mpi.h>
#include "petscvec.h"
#include "petscmat.h"
#include "petscpc.h"
#include "mfem/linalg/vector.hpp"
#include "mfem/linalg/petsc.hpp"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Parameters.hxx"
#endif /* defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC) */
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

#if defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC)

  /*!
   * \brief initialize PETSc, unless it is already initialized (by
   * mfem-mgis, or by the program, which then finalizes it)
   */
  void initializePETSc();
  /*!
   * \brief finalize PETSc at exit, before the destruction of the objects of
   * PETSc built since the start of the program, if PETSc has been
   * initialized by mfem-mgis (it is otherwise finalized with MPI, by
   * `finalize` or at the destruction of the static objects of mfem-mgis).
   *
   * The solvers and the preconditioners of PETSc call this function after
   * their first resolution on the GPU: PETSc then holds pools of events and
   * of memory of the GPU in static objects of functions, built at their
   * first use, after the static objects of mfem-mgis, and thus destroyed
   * before them at exit. `PetscFinalize`, if called afterwards, wrote in
   * these destroyed pools (heap corruption) and freed the memory of the
   * GPU with invalid streams (errors of CUDA), and the program aborted.
   */
  void finalizePETScAtExit();

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
     * \param[in] opts: options, given without the leading dash and the
     * prefix (for instance `{"ksp_type", "cg"}`)
     */
    [[nodiscard]] bool setOptions(Context& ctx, const Parameters& opts);
    /*!
     * \brief destructor, reporting the options which PETSc never read, if
     * the solver has been used, and removing the options from the options
     * database (unless PETSc has been finalized)
     */
    ~PETScOptions();
    //! \brief mark the solver as used (called by the `Mult` methods)
    void markAsUsed() const noexcept;
    /*!
     * \brief finalize PETSc at exit (see `finalizePETScAtExit`), once per
     * solver, after its first resolution on the GPU
     */
    void finalizePETScAtExitOnce() const;
    //! \brief prefix of the options
    const std::string prefix;

   private:
    //! \brief communicator of the solver
    const MPI_Comm communicator;
    //! \brief options written in the options database, with their prefix
    std::vector<std::string> options;
    //! \brief if true, the solver has been used
    mutable bool used = false;
    //! \brief if true, `finalizePETScAtExit` has been called by the solver
    mutable bool finalized_at_exit = false;
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
  /*!
   * \return a matrix of PETSc in the AIJ format of the host, copy of an
   * operator. A matrix of hypre is copied row by row, on its communicator
   * and with its distribution: MFEM converts it
   * through the `MATHYPRE` format of PETSc, whose creation initializes the
   * GPU when PETSc and hypre are built for it (a context per process, about
   * 260 MiB of the GPU and 0.4 second each with CUDA 13), although the
   * solvers of PETSc compute on the host in mfem-mgis. A matrix of PETSc is
   * referenced, and converted in place to the format of the host if needed.
   * The other operators are converted by MFEM, then to the format of the
   * host.
   * \param[in] c: communicator of the operators which are not matrices of
   * hypre
   * \param[in] op: operator
   */
  [[nodiscard]] std::unique_ptr<mfem::PetscParMatrix> makeHostMatrix(
      MPI_Comm c, const mfem::Operator& op);
  /*!
   * \return if the solvers and the preconditioners of PETSc compute on the
   * GPU: MFEM computes on a GPU with CUDA (or HIP), PETSc is built with CUDA
   * (or HIP), and the execution settings, if given, place the linear solvers
   * on the GPU. Otherwise, when MFEM computes on a GPU, they compute on the
   * host (see `makeHostMatrix` and `createHostVector`).
   */
  [[nodiscard]] bool usePETScOnDevice() noexcept;
  /*!
   * \return a matrix of PETSc in the AIJ format of the GPU (CUDA or HIP),
   * copy of an operator, when the solvers of PETSc compute on the GPU (see
   * `usePETScOnDevice`). A matrix of hypre is converted by MFEM, through
   * the `MATHYPRE` format of PETSc, where hypre stores it (on the GPU when
   * hypre computes on the GPU). A matrix of PETSc is referenced, and
   * converted in place to the format of the GPU if needed. The other
   * operators are converted by MFEM.
   * \param[in] c: communicator of the operators which are not matrices of
   * hypre
   * \param[in] op: operator
   */
  [[nodiscard]] std::unique_ptr<mfem::PetscParMatrix> makeDeviceMatrix(
      MPI_Comm c, const mfem::Operator& op);
  /*!
   * \brief create a vector of PETSc of the GPU (CUDA or HIP) whose array of
   * the GPU is the one of a vector of MFEM, without copy. The vector of PETSc
   * does not own this array, and shall be destroyed before the vector of
   * MFEM is used again, after `synchronizeDeviceVector` if PETSc wrote it.
   *
   * The vectors of PETSc of MFEM 4.10 are not used: when the device memory
   * of MFEM is managed by Umpire (the default when MFEM is built with
   * Umpire), they are vectors of the host, which the matrices and the
   * preconditioners of PETSc of the GPU write on the GPU without the
   * vectors seeing it (wrong solutions).
   * \param[out] v: vector
   * \param[in] c: communicator
   * \param[in] a: array of the GPU (`Read`, `Write` or `ReadWrite` of the
   * vector of MFEM on the device), null if the local size is null
   * \param[in] n: local size
   */
  void createDeviceVector(PETScVector& v,
                          MPI_Comm c,
                          const real* const a,
                          const int n);
  /*!
   * \brief make the array of the GPU of a vector of PETSc of the GPU up to
   * date: PETSc may have left the values in its array of the host (for
   * instance after a direct solver computing on the host)
   * \param[in] v: vector
   */
  void synchronizeDeviceVector(Vec v);
  /*!
   * \brief apply an operator of MFEM to vectors of PETSc of the GPU (shell
   * preconditioners), copied on the GPU from and to vectors of MFEM of the
   * GPU. These vectors are kept by the caller: the operators of MFEM may keep
   * aliases of the memory of the vectors given to `Mult` (the solvers of
   * hypre do, in their own vectors), which shall thus outlive the call (a
   * vector of MFEM wrapping the array of a vector of PETSc, destroyed after
   * the call, would leave dangling aliases in the memory manager of MFEM,
   * which erases a memory without checking its aliases).
   * \param[in] op: operator
   * \param[in, out] r: vector of MFEM receiving `x`
   * \param[in, out] z: vector of MFEM receiving the result
   * \param[in] x: vector
   * \param[out] y: result
   */
  void applyOnDevice(const mfem::Operator& op,
                     mfem::Vector& r,
                     mfem::Vector& z,
                     Vec x,
                     Vec y);

  /*!
   * \brief rigid body modes of the unknowns of a finite element space,
   * given to the assembled matrices of PETSc as near null space, the number
   * of components being their block size: the algebraic multigrid of PETSc
   * (GAMG) aggregates the nodes and builds its coarse spaces from them. The
   * modes are only defined when the number of components of the unknowns
   * is the space dimension (2 or 3), and given when the unknowns are ordered
   * by vector dimension, the block size of PETSc grouping consecutive rows.
   */
  struct PETScRigidBodyModes {
    /*!
     * \brief compute the coordinates of the nodes of a finite element space
     * (collective)
     * \param[in] fespace: finite element space
     */
    void initialize(FiniteElementSpace<true>& fespace);
    /*!
     * \brief set the block size and the near null space of an assembled
     * matrix (AIJ formats), if the modes are given and if the matrix has
     * the size and the distribution of the finite element space (checked on
     * all the processes)
     * \param[in] m: matrix
     */
    void attach(Mat m) const;

   private:
    //! \brief space dimension, null if the modes are not given
    int dimension = 0;
    /*!
     * \brief coordinates of the nodes, node by node, relative to the
     * centroid of the nodes
     */
    std::vector<PetscScalar> coordinates;
  };  // end of struct PETScRigidBodyModes

  /*!
   * \brief fields of the field split preconditioner of PETSc
   * (`PCFIELDSPLIT`), one per component of the unknowns of a finite element
   * space, named after the index of the component (`0`, `1`, ...): each
   * field is the set of the rows of the true degrees of freedom of a
   * component, whatever the ordering of the unknowns (`byNODES` or
   * `byVDIM`). The options of PETSc then configure each field by the prefix
   * `fieldsplit_<index>_` of the preconditioner (for instance
   * `fieldsplit_0_pc_type`).
   */
  struct PETScComponentFields {
    /*!
     * \brief compute the rows of each component owned by the process
     * \param[in] fespace: finite element space
     */
    void initialize(FiniteElementSpace<true>& fespace);
    /*!
     * \brief give the fields to a preconditioner of PETSc, if they are
     * defined (see `initialize`), if the preconditioner is a field split
     * preconditioner without fields, and if its operator has the size and
     * the distribution of the finite element space (checked on all the
     * processes). This method shall be called after the options of PETSc
     * are applied and before the setup of the preconditioner.
     * \param[in] pc: preconditioner
     */
    void attach(PC pc) const;

   private:
    /*!
     * \brief global rows of the true degrees of freedom of the process, by
     * component (empty if the fields are not defined)
     */
    std::vector<std::vector<PetscInt>> rows;
    //! \brief number of true degrees of freedom of the process
    PetscInt size = 0;
  };  // end of struct PETScComponentFields

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC) */

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_PETSCSUPPORT_HXX */
