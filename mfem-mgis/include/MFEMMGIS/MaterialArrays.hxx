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
   * GPU device, a copy on the device): gradients, thermodynamic forces,
   * internal state variables and energies at the beginning and at the end
   * of the time step, and blocks of the tangent operator.
   *
   * The material data manager of MGIS uses these arrays without copy (see
   * `getInitializer`). It is given the pointers of their copies on the
   * host, at the `HOST` location, in the point by point layout.
   *
   * Validity of the copies: MFEM keeps track of the valid copies of each
   * array. MGIS reads and writes the arrays through the pointers given at
   * construction, which MFEM does not see. Hence:
   *
   * - the code accessing the arrays through MGIS on the host (integration
   *   of the behaviours, accessors of the material, post-processings...)
   *   must first call `synchronizeOnHost`, which copies the arrays back to
   *   the host if a kernel modified them on the device, and invalidates
   *   their copies on the device;
   * - a kernel accessing the arrays on the device must get its pointers
   *   with `Read`, `Write` or `ReadWrite`, which copy the arrays to the
   *   device if needed.
   *
   * Without GPU device, both steps do nothing.
   */
  struct MFEM_MGIS_EXPORT MaterialArrays {
    /*!
     * \brief constructor
     * \param[in] b: behaviour
     * \param[in] n: number of integration points
     */
    MaterialArrays(const Behaviour&, const size_type);
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
     * before an access of the arrays by MGIS on the host.
     */
    void synchronizeOnHost();
    /*!
     * \brief makes the copies of the arrays on the host valid (copying them
     * from the device if needed), keeping the copies on the device valid,
     * before a read-only access of the arrays on the host.
     */
    void synchronizeOnHostForReading() const;
    //! \brief destructor
    ~MaterialArrays();
    //! \brief gradients, at the beginning and at the end of the time step
    std::array<mfem::Vector, 2u> gradients;
    //! \brief thermodynamic forces
    std::array<mfem::Vector, 2u> thermodynamic_forces;
    //! \brief internal state variables
    std::array<mfem::Vector, 2u> internal_state_variables;
    //! \brief stored energies (empty if the behaviour does not compute them)
    std::array<mfem::Vector, 2u> stored_energies;
    /*!
     * \brief dissipated energies (empty if the behaviour does not compute
     * them)
     */
    std::array<mfem::Vector, 2u> dissipated_energies;
    //! \brief blocks of the tangent operator
    mfem::Vector tangent_operator_blocks;

   private:
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
