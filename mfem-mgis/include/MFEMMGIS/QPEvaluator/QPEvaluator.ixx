/*!
 * \file   MFEMMGIS/QPEvaluator/QPEvaluator.ixx
 * \brief  This file implements the inline functions declared in
 * `MFEMMGIS/QPEvaluator/QPEvaluator.hxx`
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MFEM_MGIS_QPEVALUATOR_IXX
#define LIB_MFEM_MGIS_QPEVALUATOR_IXX

#include <span>
#include <cstddef>
#include <iterator>
#include <exception>
#include <algorithm>

namespace mfem_mgis::algorithm {

  /*!
   * \brief copy N values
   * \tparam N: number of values to be copied
   * \tparam InputIterator: input iterator type
   * \tparam OutputIterator: output iterator type
   * \param[in] p: iterator to the beginning of the values to be copied
   * \param[in] pe: iterator past the end of the values to be copied
   * \param[in] po: iterator to the output values
   */
  template <size_type N, typename InputIterator, typename OutputIterator>
  void copy(const InputIterator p,
            const InputIterator pe,
            const OutputIterator po) requires(N > 0) {
    if constexpr ((std::random_access_iterator<InputIterator>)&&  //
                  (std::random_access_iterator<OutputIterator>)) {
      if constexpr (N > 9) {
        std::copy(p, pe, po);
      } else if constexpr (N == 1) {
        *po = *p;
      } else if constexpr (N == 2) {
        po[0] = p[0];
        po[1] = p[1];
      } else if constexpr (N == 3) {
        po[0] = p[0];
        po[1] = p[1];
        po[2] = p[2];
      } else if constexpr (N == 4) {
        po[0] = p[0];
        po[1] = p[1];
        po[2] = p[2];
        po[3] = p[3];
      } else {
        *po = *p;
        copy<N - 1>(std::next(p), pe, std::next(po));
      }
    } else {
      std::copy(p, pe, po);
    }
  }  // end of copy

}  // end of namespace mfem_mgis::algorithm

namespace mfem_mgis::internals {

  /*!
   * \return the values of an integration point
   * \tparam N: number of components, or `dynamic_extent`
   * \param[in] values: values of all the integration points, stored point by
   * point, valid on the host
   * \param[in] i: offset of the integration point
   * \param[in] nc: number of components
   */
  template <size_type N>
  auto getIntegrationPointValues(const mfem::Vector& values,
                                 const size_type i,
                                 const size_type nc) noexcept {
    const auto o = static_cast<std::size_t>(i) * static_cast<std::size_t>(nc);
    const auto* const v = values.HostRead() + o;
    if constexpr (N == dynamic_extent) {
      return std::span<const real>(v, static_cast<std::size_t>(nc));
    } else {
      return std::span<const real, static_cast<std::size_t>(N)>(
          v, static_cast<std::size_t>(N));
    }
  }  // end of getIntegrationPointValues

}  // end of namespace mfem_mgis::internals

namespace mfem_mgis {

  inline auto RotationMatrixQPEvaluator::operator()(const size_type i) const {
    return this->material.getRotationMatrixAtIntegrationPoint(i);
  }

  inline const PartialQuadratureSpace& getSpace(
      const RotationMatrixQPEvaluator& e) {
    return e.getPartialQuadratureSpace();
  }  // end of getSpace

  inline bool check(AbstractErrorHandler& eh,
                    const RotationMatrixQPEvaluator& e) {
    return e.check(eh);
  }  // end of check

  constexpr mgis::size_type getNumberOfComponents(
      const RotationMatrixQPEvaluator&) noexcept {
    return 9u;
  }  // end of getNumberOfComponents

  template <size_type ThermodynamicForcesSize>
  RotatedThermodynamicForcesMatrixQPEvaluator<ThermodynamicForcesSize>::
      RotatedThermodynamicForcesMatrixQPEvaluator(
          const Material& m, const Material::StateSelection s) noexcept
      : material(m), stage(s) {}  // end of
                                  // RotatedThermodynamicForcesMatrixQPEvaluator

  template <size_type ThermodynamicForcesSize>
  bool
  RotatedThermodynamicForcesMatrixQPEvaluator<ThermodynamicForcesSize>::check(
      AbstractErrorHandler& ctx) const {
    if (this->material.b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      return ctx.registerErrorMessage("material is not orthotropic");
    }
    if constexpr (ThermodynamicForcesSize != dynamic_extent) {
      const auto& sm = getStateManager(this->material, this->stage);
      const auto thsize = sm.thermodynamic_forces_stride;
      if (ThermodynamicForcesSize != thsize) {
        return ctx.registerErrorMessage(
            "inconsistent number of components of the thermodynamic forces");
      }
    }
    try {
      this->material.computeThermodynamicForcesInGlobalFrame(this->values,
                                                             this->stage);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(e.what());
    } catch (...) {
      return ctx.registerErrorMessage("unknown exception");
    }
    return true;
  }

  template <size_type ThermodynamicForcesSize>
  const PartialQuadratureSpace& RotatedThermodynamicForcesMatrixQPEvaluator<
      ThermodynamicForcesSize>::getPartialQuadratureSpace() const {
    return this->material.getPartialQuadratureSpace();
  }

  template <size_type ThermodynamicForcesSize>
  size_type RotatedThermodynamicForcesMatrixQPEvaluator<
      ThermodynamicForcesSize>::getNumberOfComponents() const noexcept {
    if constexpr (ThermodynamicForcesSize == dynamic_extent) {
      const auto& sm = getStateManager(this->material, this->stage);
      return sm.thermodynamic_forces_stride;
    } else {
      return ThermodynamicForcesSize;
    }
  }

  template <size_type ThermodynamicForcesSize>
  auto RotatedThermodynamicForcesMatrixQPEvaluator<
      ThermodynamicForcesSize>::operator()(const size_type i) const {
    // values computed by `check`
    return internals::getIntegrationPointValues<ThermodynamicForcesSize>(
        this->values, i, this->getNumberOfComponents());
  }

  template <size_type ThermodynamicForcesSize>
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotatedThermodynamicForcesMatrixQPEvaluator<
          ThermodynamicForcesSize>& e) {
    return e.getPartialQuadratureSpace();
  }  // end of getSpace

  template <size_type ThermodynamicForcesSize>
  bool check(AbstractErrorHandler& eh,
             const RotatedThermodynamicForcesMatrixQPEvaluator<
                 ThermodynamicForcesSize>& e) {
    return e.check(eh);
  }  // end of check

  template <size_type ThermodynamicForcesSize>
  inline mgis::size_type getNumberOfComponents(
      const RotatedThermodynamicForcesMatrixQPEvaluator<
          ThermodynamicForcesSize>& e) noexcept {
    return e.getNumberOfComponents();
  }

  template <size_type GradientsSize>
  RotatedGradientsMatrixQPEvaluator<GradientsSize>::
      RotatedGradientsMatrixQPEvaluator(const Material& m,
                                        const Material::StateSelection s)
      : material(m), stage(s) {}  // end of RotatedGradientsMatrixQPEvaluator

  template <size_type GradientsSize>
  bool RotatedGradientsMatrixQPEvaluator<GradientsSize>::check(
      AbstractErrorHandler& ctx) const {
    if (this->material.b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      return ctx.registerErrorMessage("material is not orthotropic");
    }
    if constexpr (GradientsSize != dynamic_extent) {
      const auto& sm = getStateManager(this->material, this->stage);
      const auto gsize = sm.gradients_stride;
      if (GradientsSize != gsize) {
        return ctx.registerErrorMessage(
            "inconsistent number of components of the gradients");
      }
    }
    try {
      this->material.computeGradientsInGlobalFrame(this->values, this->stage);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(e.what());
    } catch (...) {
      return ctx.registerErrorMessage("unknown exception");
    }
    return true;
  }  // end of check

  template <size_type GradientsSize>
  const PartialQuadratureSpace& RotatedGradientsMatrixQPEvaluator<
      GradientsSize>::getPartialQuadratureSpace() const {
    return this->material.getPartialQuadratureSpace();
  }

  template <size_type GradientsSize>
  size_type RotatedGradientsMatrixQPEvaluator<
      GradientsSize>::getNumberOfComponents() const noexcept {
    if constexpr (GradientsSize == dynamic_extent) {
      const auto& sm = getStateManager(this->material, this->stage);
      return sm.gradients_stride;
    } else {
      return GradientsSize;
    }
  }

  template <size_type GradientsSize>
  auto RotatedGradientsMatrixQPEvaluator<GradientsSize>::operator()(
      const size_type i) const {
    // values computed by `check`
    return internals::getIntegrationPointValues<GradientsSize>(
        this->values, i, this->getNumberOfComponents());
  }

  template <size_type GradientsSize>
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotatedGradientsMatrixQPEvaluator<GradientsSize>& e) {
    return e.getPartialQuadratureSpace();
  }  // end of getSpace

  template <size_type GradientsSize>
  bool check(AbstractErrorHandler& eh,
             const RotatedGradientsMatrixQPEvaluator<GradientsSize>& e) {
    return e.check(eh);
  }  // end of check

  template <size_type GradientsSize>
  inline mgis::size_type getNumberOfComponents(
      const RotatedGradientsMatrixQPEvaluator<GradientsSize>& e) noexcept {
    return e.getNumberOfComponents();
  }

  template <QPEvaluatorConcept EvaluatorType1,
            QPEvaluatorConcept EvaluatorType2>
  bool checkMatchingQuadratureSpaces(Context& ctx,
                                     const EvaluatorType1& e1,
                                     const EvaluatorType2& e2) noexcept {
    const auto& qspace1 = getSpace(e1);
    const auto& qspace2 = getSpace(e2);
    if (&qspace1 != &qspace2) {
      return ctx.registerErrorMessage("unmatched quadrature spaces");
    }
    return true;
  }  // end of checkMatchingQuadratureSpaces

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_QPEVALUATOR_IXX */
