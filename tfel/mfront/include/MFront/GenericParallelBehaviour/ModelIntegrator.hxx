/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/ModelIntegrator.hxx
 * \brief  This file implements the integration strategy of the standalone
 * models treated by the `generic-parallel` interface.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_MODELINTEGRATOR_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_MODELINTEGRATOR_HXX

#include "MFront/GenericBehaviour/Integrate.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb {

  /*!
   * \brief integration strategy of the standalone models, which are treated
   * as behaviours without gradients: the outputs of the model are the
   * internal state variables and its inputs are the external state
   * variables. The model is built from the behaviour data of a point and
   * from the values of its parameters, then its functions are executed.
   */
  struct ModelIntegrator {
    /*!
     * \brief the model is executed and the state at the end of the time
     * step is updated whatever the computation options, as in the `generic`
     * interface of the models (see `updatesStateOnPrediction`)
     */
    static constexpr bool updates_state_on_prediction = true;
    /*!
     * \return the sizes of the data of a point
     * \tparam PointTraits: description of the data of a point
     */
    template <typename PointTraits>
    TFEL_HOST_DEVICE static constexpr PointDataSizes getPointDataSizes(
        const ComputationOptions&) noexcept {
      return PointDataSizes{PointTraits::value};
    }  // end of getPointDataSizes
    /*!
     * \brief execute the model at one point
     * \return 1 on success, -1 otherwise. On failure, the time step scaling
     * factor is set to 0.1, as in the `generic` interface for models. On the
     * host, the exceptions thrown by the model are reported as failures.
     *
     * A model has no prediction operator: as in the `generic` interface of
     * the models, the model is also executed if a prediction is requested.
     *
     * \tparam Model: class implementing the model
     * \tparam ParametersType: type of the structure holding the parameters
     * \param[in,out] d: behaviour data
     * \param[in] o: computation options, not used: the model is executed
     * whatever the type of computation
     * \param[in] p: out of bounds policy, which is the default policy of the
     * model, applied by the model itself
     * \param[in] parameters: values of the parameters
     */
    template <typename Model, typename ParametersType>
    TFEL_HOST_DEVICE static int integrate(
        mfront::gb::BehaviourData& d,
        const ComputationOptions&,
        const tfel::material::OutOfBoundsPolicy,
        const ParametersType& parameters) {
      // on the host, the exceptions thrown by the model (bounds checks,
      // material laws) are reported as failures
      const auto status =
          mfront::gb::callAndReportExceptions(d, [&d, &parameters] {
            const auto m = Model(d, parameters);
            return m.execute(d) ? 1 : -1;
          });
      if (status == -1) {
        *(d.rdt) = real{0.1};
      }
      return status;
    }  // end of integrate
  };

}  // end of namespace mfront::gpb

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_MODELINTEGRATOR_HXX */
