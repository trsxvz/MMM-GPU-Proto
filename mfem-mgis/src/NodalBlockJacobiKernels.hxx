/*!
 * \file   src/NodalBlockJacobiKernels.hxx
 * \brief  This file declares the kernels building the nodal blocks of the
 * Jacobi preconditioner by nodal blocks (internal header)
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_NODALBLOCKJACOBIKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_NODALBLOCKJACOBIKERNELS_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief layout of the unknowns of the nodal blocks: the component `c` of
   * the node `i` is the unknown `c * nb + i` when the unknowns are ordered
   * by nodes (`byNODES`), `i * bs + c` otherwise (`byVDIM`), `nb` being the
   * number of nodes and `bs` the size of the blocks
   */
  struct NodalBlocksLayout {
    //! \brief number of nodes
    int nb = 0;
    //! \brief size of the blocks
    int bs = 1;
    //! \brief if true, the unknowns are ordered by nodes
    bool by_nodes = true;
  };

  /*!
   * \brief launches the kernel extracting the nodal blocks of a sparse
   * matrix in the compressed row format: one thread per row
   * \param[in] l: layout of the unknowns
   * \param[in] I: offsets of the rows
   * \param[in] J: columns
   * \param[in] A: values
   * \param[out] blocks: blocks, stored column by column, one after the other
   * \param[in] use_device: if true, the kernel runs on the device
   */
  void launchNodalBlocksExtractionKernel(const NodalBlocksLayout& l,
                                         const int* const I,
                                         const int* const J,
                                         const real* const A,
                                         real* const blocks,
                                         const bool use_device);
  /*!
   * \brief launches the kernel building the nodal blocks from their columns:
   * the vector of the column `c2` of the blocks gives, at the unknown of
   * the component `c1` of a node, the entry `(c1, c2)` of the block of this
   * node. One thread per entry of the blocks.
   * \param[in] l: layout of the unknowns
   * \param[in] columns: columns of the blocks, one vector of the unknowns
   * after the other
   * \param[out] blocks: blocks, stored column by column, one after the other
   * \param[in] use_device: if true, the kernel runs on the device
   */
  void launchNodalBlocksGatherKernel(const NodalBlocksLayout& l,
                                     const real* const columns,
                                     real* const blocks,
                                     const bool use_device);

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_NODALBLOCKJACOBIKERNELS_HXX */
