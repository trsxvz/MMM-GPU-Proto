/*!
 * \file   src/NodesCoordinates.hxx
 * \brief  This file declares the computation of the coordinates of the nodes
 * of the true degrees of freedom of a finite element space, used to build
 * the rigid body modes given to the algebraic multigrids (internal header)
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_NODESCOORDINATES_HXX
#define LIB_MFEM_MGIS_SRC_NODESCOORDINATES_HXX

#include "mfem/config/config.hpp"

#ifdef MFEM_USE_MPI

#include <vector>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis::internals {

  //! \brief coordinates of the nodes of the true degrees of freedom
  struct NodesCoordinates {
    //! \brief space dimension
    int dimension = 0;
    /*!
     * \brief coordinates of the nodes of the process, component by component
     * (`dimension` arrays of the number of nodes)
     */
    std::vector<real> values;
  };

  /*!
   * \return the coordinates of the nodes of the true degrees of freedom of
   * the process, interpolated on the scalar finite element space of the
   * same collection: the node `k` carries the true degree of freedom `k` of
   * a scalar space, and the true degrees of freedom `c * n + k` (ordering
   * by nodes) or `k * vdim + c` (ordering by vector dimension) of a space
   * with `vdim` components, `n` being the number of nodes
   * \param[in] fespace: finite element space
   * \note this function is collective
   */
  [[nodiscard]] NodesCoordinates getNodesCoordinates(
      FiniteElementSpace<true>& fespace);

  /*!
   * \return true if the rigid body modes of the unknowns of a finite element
   * space are defined: the number of components of the unknowns is the space
   * dimension, 2 or 3
   * \param[in] fespace: finite element space
   */
  [[nodiscard]] bool hasRigidBodyModes(const FiniteElementSpace<true>& fespace);

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_USE_MPI */

#endif /* LIB_MFEM_MGIS_SRC_NODESCOORDINATES_HXX */
