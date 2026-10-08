/*!
 * \file   mfront/src/BehaviourRotationFunctions.cxx
 * \brief  This file implements the description of the rotation functions of
 * orthotropic behaviours.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <ostream>
#include <optional>
#include "TFEL/Material/ModellingHypothesis.hxx"
#include "MFront/SupportedTypes.hxx"
#include "MFront/VariableDescription.hxx"
#include "MFront/BehaviourDescription.hxx"
#include "MFront/BehaviourRotationFunctions.hxx"

namespace mfront {

  using Hypothesis = tfel::material::ModellingHypothesis::Hypothesis;

  //! \return true if the behaviour is a finite strain behaviour
  static bool isFiniteStrainBehaviour(const BehaviourDescription& bd) {
    return (bd.getBehaviourType() ==
            BehaviourDescription::STANDARDFINITESTRAINBEHAVIOUR) ||
           ((bd.getBehaviourType() ==
             BehaviourDescription::STANDARDSTRAINBASEDBEHAVIOUR) &&
            (bd.isStrainMeasureDefined()) &&
            (bd.getStrainMeasure() != BehaviourDescription::LINEARISED));
  }  // end of isFiniteStrainBehaviour

  //! \return the type of a mathematical object in the generated code
  static std::string getObjectType(const std::string& t,
                                   const unsigned short N) {
    return "tfel::math::" + t + "<" + std::to_string(N) + ", mfront::gb::real>";
  }  // end of getObjectType

  //! \return the type of a variable, empty for a scalar
  static std::string getVariableType(const VariableDescription& v,
                                     const unsigned short N) {
    const auto f = SupportedTypes::getTypeFlag(v.type);
    if (f == SupportedTypes::TVECTOR) {
      return getObjectType("tvector", N);
    } else if (f == SupportedTypes::STENSOR) {
      return getObjectType("stensor", N);
    } else if (f == SupportedTypes::TENSOR) {
      return getObjectType("tensor", N);
    }
    return "";
  }  // end of getVariableType

  /*!
   * \brief add the objects associated with a tangent operator block. A
   * block associated with arrays of variables is made of one sub-block for
   * each pair of entries, stored contiguously. Each sub-block is a single
   * object of TFEL if TFEL has a type for it which can be rotated, or is
   * rotated by `mfront::gb::rotateTangentOperatorBlock` otherwise
   * (derivative of a vector with respect to a symmetric tensor or a tensor,
   * and conversely, derivative of a tensor with respect to a symmetric
   * tensor).
   * \param[in,out] f: rotation function
   * \param[in] v1: first variable (the derivated one)
   * \param[in] v2: second variable
   * \param[in] N: space dimension
   * \param[in] h: modelling hypothesis
   */
  static void addTangentOperatorBlock(BehaviourRotationFunction& f,
                                      const VariableDescription& v1,
                                      const VariableDescription& v2,
                                      const unsigned short N,
                                      const Hypothesis h) {
    using Flag = SupportedTypes::TypeFlag;
    const auto f1 = SupportedTypes::getTypeFlag(v1.type);
    const auto f2 = SupportedTypes::getTypeFlag(v2.type);
    const auto is = [f1, f2](const Flag a, const Flag b) {
      return (f1 == a) && (f2 == b);
    };
    const auto type = [&]() -> std::optional<std::string> {
      if (f1 == SupportedTypes::SCALAR) {
        return getVariableType(v2, N);
      } else if (f2 == SupportedTypes::SCALAR) {
        return getVariableType(v1, N);
      } else if (is(SupportedTypes::TVECTOR, SupportedTypes::TVECTOR)) {
        return "tfel::math::tmatrix<" + std::to_string(N) + ", " +
               std::to_string(N) + ", mfront::gb::real>";
      } else if (is(SupportedTypes::STENSOR, SupportedTypes::STENSOR)) {
        return getObjectType("st2tost2", N);
      } else if (is(SupportedTypes::STENSOR, SupportedTypes::TENSOR)) {
        return getObjectType("t2tost2", N);
      } else if (is(SupportedTypes::TENSOR, SupportedTypes::TENSOR)) {
        return getObjectType("t2tot2", N);
      }
      return {};
    }();
    const auto size = [h](const VariableDescription& v) {
      return static_cast<std::size_t>(SupportedTypes::getTypeSize(v.type, 1u)
                                          .getValueForModellingHypothesis(h));
    };
    const auto s = size(v1) * size(v2);
    for (unsigned short i1 = 0; i1 != v1.arraySize; ++i1) {
      for (unsigned short i2 = 0; i2 != v2.arraySize; ++i2) {
        if (type.has_value()) {
          f.objects.push_back({.type = *type, .offset = f.size});
        } else {
          f.objects.push_back({.type = "",
                               .offset = f.size,
                               .block = BehaviourRotationFunction::Block{
                                   .type1 = getVariableType(v1, N),
                                   .type2 = getVariableType(v2, N)}});
        }
        f.size += s;
      }
    }
  }  // end of addTangentOperatorBlock

  /*!
   * \brief add the objects associated with a variable or an array of
   * variables
   * \param[in,out] f: rotation function
   * \param[in] v: variable
   * \param[in] N: space dimension
   * \param[in] h: modelling hypothesis
   */
  static void addVariable(BehaviourRotationFunction& f,
                          const VariableDescription& v,
                          const unsigned short N,
                          const Hypothesis h) {
    const auto type = getVariableType(v, N);
    const auto s =
        static_cast<std::size_t>(SupportedTypes::getTypeSize(v.type, 1u)
                                     .getValueForModellingHypothesis(h));
    for (unsigned short i = 0; i != v.arraySize; ++i) {
      f.objects.push_back({type, f.size});
      f.size += s;
    }
  }  // end of addVariable

  //! \return a rotation function of a single object
  static BehaviourRotationFunction makeRotationFunction(
      const std::string& name,
      const std::string& array_name,
      const std::string& description,
      const std::string& type,
      const std::size_t size) {
    auto f = BehaviourRotationFunction{};
    f.name = name;
    f.array_name = array_name;
    f.description = description;
    f.objects.push_back({type, 0});
    f.size = size;
    return f;
  }  // end of makeRotationFunction

  std::vector<BehaviourRotationFunction> getBehaviourRotationFunctions(
      const BehaviourDescription& bd, const Hypothesis h) {
    auto functions = std::vector<BehaviourRotationFunction>{};
    if (bd.getSymmetryType() != mfront::ORTHOTROPIC) {
      return functions;
    }
    const auto N = tfel::material::getSpaceDimension(h);
    // gradients
    auto gradients = BehaviourRotationFunction{};
    gradients.name = "rotateGradients";
    gradients.array_name = "rotateArrayOfGradients";
    gradients.description =
        "the gradients from the global frame to the material frame";
    gradients.to_material_frame = true;
    for (const auto& mv : bd.getMainVariables()) {
      addVariable(gradients, mv.first, N, h);
    }
    functions.push_back(std::move(gradients));
    if (isFiniteStrainBehaviour(bd)) {
      const auto S = tfel::material::getStensorSize(h);
      const auto T = tfel::material::getTensorSize(h);
      const auto stensor = getObjectType("stensor", N);
      const auto tensor = getObjectType("tensor", N);
      functions.push_back(makeRotationFunction(
          "rotateThermodynamicForces_CauchyStress",
          "rotateArrayOfThermodynamicForces_CauchyStress",
          "the Cauchy stress from the material frame to the global frame",
          stensor, S));
      functions.push_back(makeRotationFunction(
          "rotateThermodynamicForces_PK2Stress",
          "rotateArrayOfThermodynamicForces_PK2Stress",
          "the second Piola-Kirchhoff stress from the material frame to the "
          "global frame",
          stensor, S));
      functions.push_back(makeRotationFunction(
          "rotateThermodynamicForces_PK1Stress",
          "rotateArrayOfThermodynamicForces_PK1Stress",
          "the first Piola-Kirchhoff stress from the material frame to the "
          "global frame",
          tensor, T));
      functions.push_back(makeRotationFunction(
          "rotateTangentOperatorBlocks_dsig_dF",
          "rotateArrayOfTangentOperatorBlocks_dsig_dF",
          "the derivative of the Cauchy stress with respect to the "
          "deformation gradient from the material frame to the global frame",
          getObjectType("t2tost2", N), S * T));
      functions.push_back(makeRotationFunction(
          "rotateTangentOperatorBlocks_dPK1_dF",
          "rotateArrayOfTangentOperatorBlocks_dPK1_dF",
          "the derivative of the first Piola-Kirchhoff stress with respect to "
          "the deformation gradient from the material frame to the global "
          "frame",
          getObjectType("t2tot2", N), T * T));
      functions.push_back(makeRotationFunction(
          "rotateTangentOperatorBlocks_dPK2_dEGL",
          "rotateArrayOfTangentOperatorBlocks_dPK2_dEGL",
          "the derivative of the second Piola-Kirchhoff stress with respect "
          "to the Green-Lagrange strain from the material frame to the global "
          "frame",
          getObjectType("st2tost2", N), S * S));
      functions.push_back(makeRotationFunction(
          "rotateTangentOperatorBlocks_dtau_ddF",
          "rotateArrayOfTangentOperatorBlocks_dtau_ddF",
          "the derivative of the Kirchhoff stress with respect to the spatial "
          "increment of the deformation gradient from the material frame to "
          "the global frame",
          getObjectType("t2tost2", N), S * T));
      return functions;
    }
    // thermodynamic forces
    auto forces = BehaviourRotationFunction{};
    forces.name = "rotateThermodynamicForces";
    forces.array_name = "rotateArrayOfThermodynamicForces";
    forces.description =
        "the thermodynamic forces from the material frame to the global "
        "frame";
    for (const auto& mv : bd.getMainVariables()) {
      addVariable(forces, mv.second, N, h);
    }
    functions.push_back(std::move(forces));
    // tangent operator blocks
    auto blocks = BehaviourRotationFunction{};
    blocks.name = "rotateTangentOperatorBlocks";
    blocks.array_name = "rotateArrayOfTangentOperatorBlocks";
    blocks.description =
        "the tangent operator blocks from the material frame to the global "
        "frame";
    for (const auto& [v1, v2] : bd.getTangentOperatorBlocks()) {
      addTangentOperatorBlock(blocks, v1, v2, N, h);
    }
    functions.push_back(std::move(blocks));
    return functions;
  }  // end of getBehaviourRotationFunctions

  void writePointRotation(std::ostream& os,
                          const BehaviourRotationFunction& f) {
    // the components of the rotation matrix `rv` are read as the components
    // of the matrix `r` stored by rows, so that `r` is the transpose of the
    // rotation matrix from the global frame to the material frame, which is
    // stored in column-major format. `r` is transposed for the rotations to
    // the global frame.
    const auto transposed = !f.to_material_frame;
    os << "const auto r = tfel::math::rotation_matrix<mfront::gb::real>{";
    for (unsigned short i = 0; i != 3; ++i) {
      for (unsigned short j = 0; j != 3; ++j) {
        os << ((i == 0) && (j == 0) ? "" : ", ") << "rv["
           << (transposed ? 3 * j + i : 3 * i + j) << "]";
      }
    }
    os << "};\n";
    for (const auto& o : f.objects) {
      const auto offset = std::to_string(o.offset);
      if (o.block.has_value()) {
        os << "mfront::gb::rotateTangentOperatorBlock<" << o.block->type1
           << ", " << o.block->type2 << ">(dest + " << offset << ", src + "
           << offset << ", r);\n";
      } else if (o.type.empty()) {
        os << "dest[" << offset << "] = src[" << offset << "];\n";
      } else {
        // the source is copied before the destination is written, so that
        // the destination may be the source
        os << "tfel::math::map<" << o.type << ">(dest + " << offset
           << ") = tfel::math::change_basis(" << o.type << "{src + " << offset
           << "}, r);\n";
      }
    }
  }  // end of writePointRotation

}  // end of namespace mfront
