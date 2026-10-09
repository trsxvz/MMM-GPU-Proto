/*!
 * \file   include/MFEMMGIS/MaterialArrays.hxx
 * \brief  This file declares the `MaterialArrays` class
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_MATERIALARRAYS_HXX
#define LIB_MFEM_MGIS_MATERIALARRAYS_HXX

#include <array>
#include "mfem/linalg/vector.hpp"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Behaviour.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"

namespace mfem_mgis {

  /*!
   * \brief arrays of the data of a material at its integration points, in
   * the memory of MFEM (`mfem::Memory`, with a copy on the host and, with a
   * device, a copy on the device): gradients, thermodynamic forces,
   * internal state variables and energies at the beginning and at the end
   * of the time step, and blocks of the tangent operator.
   *
   * The material data manager of MGIS uses these arrays without copy (see
   * `getInitializer`). It is given the pointers of their copies on the
   * host, at the `HOST` location, in the layout chosen at construction
   * (point by point by default, or component by component, which gives
   * coalesced accesses on a GPU). The arrays are never reallocated: they
   * are private, and only pointers to their values are given.
   *
   * Validity of the copies: MFEM keeps track of the valid copies of each
   * array, but MGIS (and the views returned by the accessors of the
   * material) read and write the arrays through the pointers of the host,
   * which MFEM does not see. Hence:
   *
   * - a phase of accesses on the host (integration of the behaviours by
   *   MGIS, update, revert, accessors of the material, post-processings...)
   *   must start with `synchronizeOnHost`, which copies the arrays back to
   *   the host if a kernel modified them on the device, and invalidates
   *   their copies on the device, or with `synchronizeOnHostForReading` for
   *   read-only accesses, which keeps the copies on the device valid;
   * - no access on the device may happen during such a phase: a kernel
   *   which reads an array (`Read`) after the synchronization makes its
   *   copy on the device valid again, and a later write on the host would
   *   leave this copy valid but outdated. The views returned by the
   *   accessors are only valid until the next access on the device;
   * - a kernel accessing the arrays on the device gets its pointers with
   *   `Read`, `Write` or `ReadWrite`, which copy the arrays to the device
   *   if needed.
   *
   * Without device (`cpu`, `omp`), both steps do nothing. With the `debug`
   * device of MFEM, the copies are made between the memory of the host and
   * a separate memory, which detects the violations of these rules.
   */
  struct MFEM_MGIS_EXPORT MaterialArrays {
    //! \brief arrays of a state of the material
    enum struct StateArray {
      GRADIENTS,
      THERMODYNAMIC_FORCES,
      INTERNAL_STATE_VARIABLES,
      STORED_ENERGIES,
      DISSIPATED_ENERGIES
    };
    /*!
     * \brief constructor
     * \param[in] b: behaviour
     * \param[in] n: number of integration points
     * \param[in] l: layout of the arrays
     */
    MaterialArrays(const Behaviour&,
                   const size_type,
                   const mgis::FieldLayout = mgis::FieldLayout::POINT_BY_POINT);
    /*!
     * \return the initializer of the material data manager of MGIS using
     * these arrays
     * \note the host copies of the arrays are made valid, and the copies on
     * the device invalid (see `synchronizeOnHost`).
     */
    mgis::behaviour::MaterialDataManagerInitializer getInitializer();
    /*!
     * \brief makes the copies of the arrays on the host valid (copying them
     * from the device if needed) and invalidates the copies on the device,
     * at the beginning of a phase of accesses on the host.
     */
    void synchronizeOnHost();
    /*!
     * \brief makes the copies of the arrays on the host valid (copying them
     * from the device if needed), keeping the copies on the device valid,
     * at the beginning of a phase of read-only accesses on the host.
     */
    void synchronizeOnHostForReading() const;
    /*!
     * \brief makes the copy of one array of a state on the host valid, and
     * invalidates its copy on the device (see `synchronizeOnHost`)
     * \param[in] a: array
     * \param[in] s: state
     */
    void synchronizeOnHost(const StateArray, const TimeStepStage);
    /*!
     * \brief makes the copy of one array of a state on the host valid,
     * keeping its copy on the device valid (see
     * `synchronizeOnHostForReading`)
     * \param[in] a: array
     * \param[in] s: state
     */
    void synchronizeOnHostForReading(const StateArray,
                                     const TimeStepStage) const;
    /*!
     * \brief makes the copy of the blocks of the tangent operator on the
     * host valid, keeping their copy on the device valid
     */
    void synchronizeTangentOperatorBlocksOnHostForReading() const;
    /*!
     * \return an array of a state, for read-only accesses through `Read`
     * or `HostRead`, which keep track of the valid copies. The pointer
     * returned by `GetData` must not be used: it is not tracked. The array
     * is empty if the behaviour does not use it.
     * \param[in] a: array
     * \param[in] s: state
     */
    const mfem::Vector& getArray(const StateArray, const TimeStepStage) const;
    //! \return the blocks of the tangent operator, for read-only accesses
    const mfem::Vector& getTangentOperatorBlocks() const;
    /*!
     * \return a pointer to the values of an array of a state, for reading
     * and writing them on the device (or on the host): see
     * `mfem::Vector::ReadWrite`
     * \param[in] a: array
     * \param[in] s: state
     * \param[in] d: if true, the values are accessed on the device (if a
     * device is enabled)
     */
    real* readWrite(const StateArray, const TimeStepStage, const bool = true);
    /*!
     * \return a pointer to the values of an array of a state, for writing
     * them on the device (or on the host): see `mfem::Vector::Write`
     * \param[in] a: array
     * \param[in] s: state
     * \param[in] d: if true, the values are accessed on the device (if a
     * device is enabled)
     */
    real* write(const StateArray, const TimeStepStage, const bool = true);
    /*!
     * \return a pointer to the blocks of the tangent operator, for reading
     * and writing them on the device (or on the host)
     * \param[in] d: if true, the values are accessed on the device (if a
     * device is enabled)
     */
    real* readWriteTangentOperatorBlocks(const bool = true);
    /*!
     * \return a pointer to the blocks of the tangent operator, for writing
     * them on the device (or on the host)
     * \param[in] d: if true, the values are accessed on the device (if a
     * device is enabled)
     */
    real* writeTangentOperatorBlocks(const bool = true);
    //! \brief destructor
    ~MaterialArrays();
    //! \brief layout of the arrays
    const mgis::FieldLayout array_layout;

   private:
    //! \return an array of a state
    mfem::Vector& getMutableArray(const StateArray, const TimeStepStage);
    //! \brief arrays of the two states, indexed by `StateArray`
    std::array<std::array<mfem::Vector, 2u>, 5u> state_arrays;
    //! \brief blocks of the tangent operator
    mfem::Vector tangent_operator_blocks;
    //! \brief copy constructor (disabled)
    MaterialArrays(const MaterialArrays&) = delete;
    //! \brief move constructor (disabled)
    MaterialArrays(MaterialArrays&&) = delete;
    //! \brief standard assignment (disabled)
    MaterialArrays& operator=(const MaterialArrays&) = delete;
    //! \brief move assignment (disabled)
    MaterialArrays& operator=(MaterialArrays&&) = delete;
  };  // end of struct MaterialArrays

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MATERIALARRAYS_HXX */
