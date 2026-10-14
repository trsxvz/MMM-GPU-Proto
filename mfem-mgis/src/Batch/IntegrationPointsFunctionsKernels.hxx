/*!
 * \file   src/Batch/IntegrationPointsFunctionsKernels.hxx
 * \brief  This file declares the kernels computing integrals, right-hand
 * sides of L2 projections and local projections of functions defined at
 * the integration points of a material, and mass matrices of the elements
 * (internal header)
 * \date   11/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_INTEGRATIONPOINTSFUNCTIONSKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_INTEGRATIONPOINTSFUNCTIONSKERNELS_HXX

#include "mfem/general/array.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::batch::internals {

  //! \brief geometric data of an integration rule used by the kernels
  struct IntegrationRuleKernelArguments {
    /*!
     * \brief values of the shape functions on the reference element (`nq x
     * nd`)
     */
    const real* shape_functions = nullptr;
    //! \brief determinants of the jacobians (`nq x ne`)
    const real* determinants = nullptr;
    //! \brief weights of the integration points on the reference element
    const real* weights = nullptr;
    //! \brief number of integration points per element
    int nq = 0;
  };

  //! \brief arguments of the kernels on functions at the integration points
  struct IntegrationPointsFunctionsKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    //! \brief geometric data of the integration rule of the material
    IntegrationRuleKernelArguments points;
    /*!
     * \brief values of the function (`nc x nq x number_of_elements`), or
     * null for the measure
     */
    const real* values = nullptr;
    //! \brief number of components of the function
    int nc = 1;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of degrees of freedom per element and per component
    int nd = 0;
    /*!
     * \brief if true, the geometric data are the ones of the elements of the
     * material only, the `k`-th element of the material being the `k`-th
     * element of the geometric data
     */
    bool compact_geometry = false;
  };

  /*!
   * \return the integral of a scalar function over the elements of a
   * material, or their measure if no values are given. The contributions of
   * the points are accumulated in double, whatever the numerical type of
   * mfem-mgis: in float, a sum over all the points would lose about the
   * machine epsilon times the number of points.
   * \param[in] a: arguments
   * \param[in, out] w: workspace of the reduction
   * \param[in] use_device: if true, the kernel runs on the device (if a device
   * is enabled), and the pointers of the arguments are located on the device
   */
  [[nodiscard]] double launchIntegralKernel(
      const IntegrationPointsFunctionsKernelArguments& a,
      mfem::Array<double>& w,
      const bool use_device);

  /*!
   * \brief launches the kernel adding the right-hand side of the L2
   * projection of a component of a function to values element by element
   * (`nd x ne`)
   * \param[in] a: arguments
   * \param[in] c: component
   * \param[in, out] r: values element by element
   * \param[in] use_device: if true, the kernel runs on the device (if a device
   * is enabled), and the pointers of the arguments are located on the device
   */
  void launchL2ProjectionRightHandSideKernel(
      const IntegrationPointsFunctionsKernelArguments& a,
      const int c,
      real* const r,
      const bool use_device);

  /*!
   * \brief launches the kernel computing the values at the nodes of the
   * local projections of a function (`nd x nc x ne`), the number of degrees
   * of freedom per element being at most
   * `maximumNumberOfNodesOfLocalProjections`
   * \param[in] a: arguments
   * \param[in] m: geometric data of the integration rule of the mass
   * matrices
   * \param[in, out] r: values element by element
   * \param[in] use_device: if true, the kernel runs on the device (if a device
   * is enabled), and the pointers of the arguments are located on the device
   */
  void launchLocalProjectionsKernel(
      const IntegrationPointsFunctionsKernelArguments& a,
      const IntegrationRuleKernelArguments& m,
      real* const r,
      const bool use_device);

  /*!
   * \brief launches the kernel computing, or adding, the mass matrices of
   * all the elements (`nd x nd x ne`)
   * \param[in] m: geometric data of the integration rule
   * \param[in] p: native index of each degree of freedom of the mass
   * matrices (lexicographic ordering), or null for the native ordering
   * \param[in] nd: number of degrees of freedom per element
   * \param[in] ne: number of elements
   * \param[in, out] r: mass matrices
   * \param[in] add: if true, the mass matrices are added to `r`
   * \param[in] use_device: if true, the kernel runs on the device (if a device
   * is enabled), and the pointers of the arguments are located on the device
   */
  void launchMassMatricesKernel(const IntegrationRuleKernelArguments& m,
                                const int* const p,
                                const int nd,
                                const int ne,
                                real* const r,
                                const bool add,
                                const bool use_device);

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_INTEGRATIONPOINTSFUNCTIONSKERNELS_HXX */
