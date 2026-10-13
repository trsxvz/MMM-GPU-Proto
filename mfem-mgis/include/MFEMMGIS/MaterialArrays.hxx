/*!
 * \file   include/MFEMMGIS/MaterialArrays.hxx
 * \brief  This file declares the `MaterialArrays` class
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_MATERIALARRAYS_HXX
#define LIB_MFEM_MGIS_MATERIALARRAYS_HXX

#include <span>
#include <array>
#include <string>
#include "mfem/linalg/vector.hpp"
#include "MGIS/Raise.hxx"
#include "MGIS/FieldView.hxx"
#include "MGIS/MemoryLocation.hxx"
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
   * which MFEM does not see, except for the behaviours and the point-wise
   * models integrated on the device, whose views are built from pointers
   * obtained by `Read`, `Write` or `ReadWrite` (see
   * `Material::makeDeviceView`).
   * Hence:
   *
   * - a phase of accesses on the host (integration of the behaviours by
   *   MGIS on the host, update and revert of these behaviours, accessors of
   *   the material, post-processings...) must start with
   *   `synchronizeOnHost`, which copies the arrays back to the host if a
   *   kernel modified them on the device, and invalidates their copies on
   *   the device, or with `synchronizeOnHostForReading` for read-only
   *   accesses, which keeps the copies on the device valid;
   * - no access on the device may happen during such a phase: a kernel
   *   which reads an array (`Read`) after the synchronization makes its
   *   copy on the device valid again, and a later write on the host would
   *   leave this copy valid but outdated. The views returned by the
   *   accessors are only valid until the next access on the device;
   * - a kernel accessing the arrays on the device gets its pointers with
   *   `Read`, `Write` or `ReadWrite`, which copy the arrays to the device
   *   if needed;
   * - the update and the revert of a behaviour integrated on the device are
   *   made on the device, then the arrays are synchronized on the host (see
   *   `BehaviourIntegratorBase::update`); the ones of a point-wise model
   *   are made on the host (see `PointWiseModel::update`), its arrays being
   *   made valid on the host for reading after its integration on the
   *   device.
   *
   * Without device (`cpu`, `omp`), the synchronizations do nothing, the
   * copy on the host being the only one. With the `debug` device of MFEM,
   * the copies are made between the memory of the host and a separate
   * memory, which detects the violations of these rules.
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
    MaterialArrays(
        const Behaviour& b,
        const size_type n,
        const mgis::FieldLayout l = mgis::FieldLayout::POINT_BY_POINT);
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
    void synchronizeOnHost(const StateArray a, const TimeStepStage s);
    /*!
     * \brief makes the copy of one array of a state on the host valid,
     * keeping its copy on the device valid (see
     * `synchronizeOnHostForReading`)
     * \param[in] a: array
     * \param[in] s: state
     */
    void synchronizeOnHostForReading(const StateArray a,
                                     const TimeStepStage s) const;
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
    const mfem::Vector& getArray(const StateArray a,
                                 const TimeStepStage s) const;
    //! \return the blocks of the tangent operator, for read-only accesses
    const mfem::Vector& getTangentOperatorBlocks() const;
    /*!
     * \return a pointer to the values of an array of a state, for reading
     * them on the device (or on the host): see `mfem::Vector::Read`. The
     * copy made, if any, is counted (see `getMaterialArraysTransfers`).
     * \param[in] a: array
     * \param[in] s: state
     * \param[in] d: if true, the values are accessed on the device (if a
     * device is enabled)
     */
    const real* read(const StateArray a,
                     const TimeStepStage s,
                     const bool d = true) const;
    /*!
     * \return a pointer to the blocks of the tangent operator, for reading
     * them on the device (or on the host), the copy made being counted
     * \param[in] d: if true, the values are accessed on the device (if a
     * device is enabled)
     */
    const real* readTangentOperatorBlocks(const bool d = true) const;
    /*!
     * \return a pointer to the values of an array of a state, for reading
     * and writing them on the device (or on the host): see
     * `mfem::Vector::ReadWrite`
     * \param[in] a: array
     * \param[in] s: state
     * \param[in] d: if true, the values are accessed on the device (if a
     * device is enabled)
     */
    real* readWrite(const StateArray a,
                    const TimeStepStage s,
                    const bool d = true);
    /*!
     * \return a pointer to the values of an array of a state, for writing
     * them on the device (or on the host): see `mfem::Vector::Write`
     * \param[in] a: array
     * \param[in] s: state
     * \param[in] d: if true, the values are accessed on the device (if a
     * device is enabled)
     */
    real* write(const StateArray a, const TimeStepStage s, const bool d = true);
    /*!
     * \return a pointer to the blocks of the tangent operator, for reading
     * and writing them on the device (or on the host)
     * \param[in] d: if true, the values are accessed on the device (if a
     * device is enabled)
     */
    real* readWriteTangentOperatorBlocks(const bool d = true);
    /*!
     * \return a pointer to the blocks of the tangent operator, for writing
     * them on the device (or on the host)
     * \param[in] d: if true, the values are accessed on the device (if a
     * device is enabled)
     */
    real* writeTangentOperatorBlocks(const bool d = true);
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

  /*!
   * \brief copies of the arrays of the materials between the host and the
   * device, made by the accesses of `MaterialArrays` (reads, reads and
   * writes, synchronizations on the host) when the accessed copy was not
   * valid. The copies made by the reads through `getArray` or
   * `getTangentOperatorBlocks` are not counted, nor the copies to the device
   * of the material properties, of the external state variables and of the
   * mass density made at each integration on the device (see
   * `Material::makeDeviceView`), nor the ones of the rotation matrices.
   */
  struct MaterialArraysTransfers {
    //! \brief number of copies from the host to the device
    size_type host_to_device_copies = 0;
    //! \brief number of bytes copied from the host to the device
    std::size_t host_to_device_bytes = 0;
    //! \brief time spent in the copies from the host to the device, in s
    double host_to_device_time = 0;
    //! \brief number of copies from the device to the host
    size_type device_to_host_copies = 0;
    //! \brief number of bytes copied from the device to the host
    std::size_t device_to_host_bytes = 0;
    //! \brief time spent in the copies from the device to the host, in s
    double device_to_host_time = 0;
  };

  /*!
   * \return the copies of the arrays of the materials between the host and
   * the device made by the current process since the start of the program.
   * They measure the transfers of the placements of the computations on
   * both processors (see `ExecutionSettings`), for instance the finite
   * elements on the GPU and the behaviours on the CPU.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] MaterialArraysTransfers
  getMaterialArraysTransfers() noexcept;

  /*!
   * \return the location, for MGIS, of the copies of the arrays on the device
   * enabled in MFEM (`Read`, `Write`, `ReadWrite`): the memory of the host
   * for the `debug` device of MFEM, whose copies on the device are in the
   * memory of the host, the memory of the device otherwise
   * \pre a device is enabled in MFEM (`DEVICE` is also returned otherwise)
   */
  MFEM_MGIS_EXPORT [[nodiscard]] mgis::MemoryLocation
  getDeviceMemoryLocation() noexcept;

  /*!
   * \return a view of the values of a view covering an array, in a copy of
   * this array (for instance its copy on the device), with the same layout
   * \param[in] v: view of the array (an empty view is returned unchanged)
   * \param[in] p: address of the copy of the array
   * \param[in] n: size of the copy, which must be the size of the array
   */
  template <typename ValueType>
  [[nodiscard]] mgis::FieldView<ValueType> relocate(
      const mgis::FieldView<ValueType>& v, ValueType* const p, const int n) {
    if (v.empty()) {
      return v;
    }
    if (static_cast<std::size_t>(n) != v.getValues().size()) {
      mgis::raise("relocate: unexpected size of an array (" +
                  std::to_string(n) + " values instead of " +
                  std::to_string(v.getValues().size()) + ")");
    }
    return mgis::FieldView<ValueType>(
        std::span<ValueType>(p, v.getValues().size()), v.getNumberOfPoints(),
        v.getNumberOfComponents(), v.getPointStride(), v.getComponentStride());
  }  // end of relocate

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MATERIALARRAYS_HXX */
