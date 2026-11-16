/*!
 * \file   src/Batch/GlobalStiffnessMatrixKernels.hxx
 * \brief  This file declares the kernels assembling the values of the
 * global matrix of a parallel finite element space from the local matrix of
 * the process and from the values received from the other processes
 * (internal header)
 * \date   15/11/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_GLOBALSTIFFNESSMATRIXKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_GLOBALSTIFFNESSMATRIXKERNELS_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::batch::internals {

  /*!
   * \brief arguments of the kernel computing the values of the part of the
   * global matrix owned by the process (see `GlobalStiffnessMatrix`)
   */
  struct GlobalStiffnessMatrixKernelArguments {
    /*!
     * \brief offsets of the contributions to each entry of the global
     * matrix (entries of the diagonal block, then of the off diagonal
     * block)
     */
    const int* offsets = nullptr;
    /*!
     * \brief contributions to the entries of the global matrix: indices of
     * entries of the local matrix, then, from the number of entries of the
     * local matrix, indices of received values
     */
    const int* contributions = nullptr;
    //! \brief values of the local matrix
    const real* local_values = nullptr;
    //! \brief values received from the other processes
    const real* received_values = nullptr;
    //! \brief values of the diagonal block of the global matrix
    real* diagonal_block = nullptr;
    //! \brief values of the off diagonal block of the global matrix
    real* off_diagonal_block = nullptr;
    //! \brief number of entries of the local matrix
    int local_entries = 0;
    //! \brief number of entries of the diagonal block
    int diagonal_block_entries = 0;
    //! \brief number of entries of the off diagonal block
    int off_diagonal_block_entries = 0;
  };

  /*!
   * \brief sets values to one (matrix giving the pattern of the global
   * matrix)
   * \param[out] v: values
   * \param[in] n: number of values
   * \param[in] use_device: if true, the kernel is run on the device and
   * the pointer is located on the device
   */
  void launchGlobalStiffnessMatrixOnesKernel(real* const,
                                             const int,
                                             const bool);

  /*!
   * \brief copies the values of the local matrix sent to the other
   * processes
   * \param[out] s: values sent
   * \param[in] entries: entries of the local matrix sent
   * \param[in] a: values of the local matrix
   * \param[in] n: number of values sent
   * \param[in] use_device: if true, the kernel is run on the device and
   * the pointers are located on the device
   */
  void launchGlobalStiffnessMatrixSentValuesKernel(
      real* const, const int* const, const real* const, const int, const bool);

  /*!
   * \brief computes each value of the part of the global matrix owned by
   * the process as the sum of its contributions, added in their order
   * \param[in] a: arguments
   * \param[in] use_device: if true, the kernel is run on the device and
   * the pointers are located on the device
   */
  void launchGlobalStiffnessMatrixKernel(
      const GlobalStiffnessMatrixKernelArguments&, const bool);

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_GLOBALSTIFFNESSMATRIXKERNELS_HXX */
