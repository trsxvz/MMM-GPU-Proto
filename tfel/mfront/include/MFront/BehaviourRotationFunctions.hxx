/*!
 * \file   mfront/include/MFront/BehaviourRotationFunctions.hxx
 * \brief  This file declares the description of the rotation functions of
 * orthotropic behaviours, shared by the `generic` and `generic-parallel`
 * interfaces.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_BEHAVIOURROTATIONFUNCTIONS_HXX
#define LIB_MFRONT_BEHAVIOURROTATIONFUNCTIONS_HXX

#include <string>
#include <vector>
#include <optional>
#include <cstddef>
#include <iosfwd>
#include "TFEL/Material/ModellingHypothesis.hxx"
#include "MFront/MFrontConfig.hxx"

namespace mfront {

  // forward declaration
  struct BehaviourDescription;

  /*!
   * \brief description of a function rotating the values of an integration
   * point (gradients, thermodynamic forces or tangent operator blocks) of an
   * orthotropic behaviour.
   *
   * The rotation matrix `rv` given to those functions is the rotation matrix
   * from the global frame to the material frame, stored in column-major
   * format. The
   * gradients are rotated from the global frame to the material frame, the
   * thermodynamic forces and the tangent operator blocks from the material
   * frame to the global frame.
   */
  struct MFRONT_VISIBILITY_EXPORT BehaviourRotationFunction {
    /*!
     * \brief a tangent operator block for which TFEL has no type which can
     * be rotated (derivative of a vector with respect to a symmetric tensor
     * or a tensor, and conversely, derivative of a tensor with respect to a
     * symmetric tensor), rotated by `mfront::gb::rotateTangentOperatorBlock`
     */
    struct Block {
      //! \brief type of the derivated object
      std::string type1;
      //! \brief type of the object with respect to which the derivative is
      //! taken
      std::string type2;
    };
    //! \brief a rotated object
    struct Object {
      /*!
       * \brief type of the object, empty for a scalar or for a block for
       * which TFEL has no type (see `block`, which is tested first)
       */
      std::string type;
      //! \brief offset of the object in the values of a point
      std::size_t offset;
      //! \brief description of a block for which TFEL has no type
      std::optional<Block> block = {};
    };
    //! \brief name of the function, appended to the name of the behaviour
    std::string name;
    /*!
     * \brief name of the function rotating an array of values, appended to
     * the name of the behaviour
     */
    std::string array_name;
    //! \brief description of the rotated values (used in the documentation)
    std::string description;
    //! \brief rotated objects
    std::vector<Object> objects;
    //! \brief number of values of a point
    std::size_t size = 0;
    //! \brief true if the values are rotated to the material frame
    bool to_material_frame = false;
  };

  /*!
   * \return the rotation functions of a behaviour for the given modelling
   * hypothesis (empty if the behaviour is not orthotropic)
   * \param[in] bd: behaviour description
   * \param[in] h: modelling hypothesis
   */
  MFRONT_VISIBILITY_EXPORT std::vector<BehaviourRotationFunction>
  getBehaviourRotationFunctions(
      const BehaviourDescription&,
      const tfel::material::ModellingHypothesis::Hypothesis);

  /*!
   * \brief write the statements rotating the values of one point. Those
   * statements use the `dest`, `src` and `rv` pointers, which have the
   * meaning of the arguments of the rotation functions of the `generic`
   * interface. `dest` and `src` may be equal. The generated code can be used
   * in device code.
   * \param[out] os: output stream
   * \param[in] f: rotation function
   */
  MFRONT_VISIBILITY_EXPORT void writePointRotation(
      std::ostream&, const BehaviourRotationFunction&);

}  // end of namespace mfront

#endif /* LIB_MFRONT_BEHAVIOURROTATIONFUNCTIONS_HXX */
