/*!
 * \file   include/MFEMMGIS/BehaviourIntegratorTraits.hxx
 * \brief  This file declares the `BehaviourIntegratorTraits` class
 * \author Thomas Helfer
 * \date   06/04/2021
 */

#ifndef LIB_MFEMMGIS_BEHAVIOURINTEGRATORTRAITS_HXX
#define LIB_MFEMMGIS_BEHAVIOURINTEGRATORTRAITS_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief a traits class aimed at specializing the
   * `StandardBehaviourIntegratorCRTPBase` and
   * `FBarBehaviourIntegratorCRTPBase` classes.
   *
   * \note The computations by kernels of the batch path
   * (`StandardBehaviourIntegratorCRTPBase::supportsKernels`) are chosen from
   * these traits and from the structure of the behaviour (see
   * `batch::getGradientsKind`): they assume the standard computations of
   * the gradients, of the inner forces and of the stiffness matrices for
   * each kind of gradients, and the weights of the integration points of
   * the cartesian hypotheses. A behaviour integrator whose computations
   * differ (axisymmetrical weights, other gradients...) must override
   * `supportsKernels`.
   */
  template <typename BehaviourIntegrator>
  struct BehaviourIntegratorTraits {
    //! \brief number of components of the unknowns
    static constexpr size_type unknownsSize = 0u;
    //! \brief if the computation of the gradients requires the shape functions
    static constexpr bool gradientsComputationRequiresShapeFunctions = false;
    /*!
     * \brief if the computation of the gradients requires the derivatives of
     * the shape functions
     */
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = false;
    //! \brief if the external state variables are updated from the unknowns
    static constexpr bool updateExternalStateVariablesFromUnknownsValues =
        false;
  };  // end of struct BehaviourIntegratorTraits

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_BEHAVIOURINTEGRATORTRAITS_HXX */
