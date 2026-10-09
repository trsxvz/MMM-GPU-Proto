/*!
 * \file   src/BehaviourIntegratorBase.cxx
 * \brief  This file implements the `BehaviourIntegratorBase` class
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include <span>
#include <algorithm>
#include <limits>
#include <utility>
#include <variant>
#include "mfem/general/device.hpp"
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"
#include "MGIS/Behaviour/MaterialDataView.hxx"
#include "MFEMMGIS/Batch/MFEMFieldOperations.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/QPEvaluator/AbstractQPEvaluator.hxx"
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  [[nodiscard]] static std::string prependBehaviourDescriptionToErrorMessage(
      const Behaviour& b, const std::string& e) {
    return "error while treating behaviour '" + b.behaviour +
           "' implemented by '" + b.function + "' in library '" + b.library +
           "': " + e;
  }  // end of prependBehaviourDescriptionToErrorMessage

  //! \return a textual representation of a list of sizes
  [[nodiscard]] static std::string convertSizesToString(
      const std::vector<size_type>& sizes) {
    auto r = std::string{};
    for (const auto s : sizes) {
      r += r.empty() ? "" : ", ";
      r += std::to_string(s);
    }
    return "(" + r + ")";
  }  // end of convertSizesToString

  //! \return a textual representation of a list of sizes of tangent blocks
  [[nodiscard]] static std::string convertSizesOfTangentOperatorBlocksToString(
      const std::vector<std::pair<size_type, size_type>>& sizes) {
    auto r = std::string{};
    for (const auto& [s1, s2] : sizes) {
      r += r.empty() ? "" : ", ";
      r += std::to_string(s1) + "x" + std::to_string(s2);
    }
    return "(" + r + ")";
  }  // end of convertSizesOfTangentOperatorBlocksToString

  BehaviourIntegratorBase::BehaviourIntegratorBase(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::unique_ptr<const Behaviour> b_ptr)
      : Material(s, std::move(b_ptr)) {
    // The following arrays are storing material properties and external
    // state variables. They can be allocated a single time
    this->wks.mps.resize(getArraySize(this->b.mps, this->b.hypothesis));
    this->wks.esvs0.resize(getArraySize(this->b.esvs, this->b.hypothesis));
    this->wks.esvs1.resize(getArraySize(this->b.esvs, this->b.hypothesis));
  }  // end of BehaviourIntegratorBase

  void BehaviourIntegratorBase::checkHypothesis(attributes::Throwing,
                                                const Hypothesis h) const {
    using namespace mgis::behaviour;
    if (this->b.hypothesis != h) {
      const auto h1 = std::string(toString(this->b.hypothesis));
      const auto h2 = std::string(toString(h));
      raise(prependBehaviourDescriptionToErrorMessage(
          this->b, "the behaviour hypothesis (" + h1 + ") does not match the " +
                       "integrator hypothesis (" + h2 + ")"));
    }
  }  // end of BehaviourIntegratorBase::checkHypothesis

  void BehaviourIntegratorBase::checkIfAFiniteStrainBehaviourIsDeclared(
      attributes::Throwing) const {
    auto isDeformationGradient = [](const mgis::behaviour::Variable& g) {
      if (g.name != "DeformationGradient") {
        return false;
      }
      if (g.type != mgis::behaviour::Variable::TENSOR) {
        return false;
      }
      return true;
    };
    auto isPK1 = [](const mgis::behaviour::Variable& g) {
      if (g.name != "FirstPiolaKirchhoffStress") {
        return false;
      }
      if (g.type != mgis::behaviour::Variable::TENSOR) {
        return false;
      }
      return true;
    };
    if ((this->b.btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) &&
        (this->b.btype != Behaviour::GENERALBEHAVIOUR)) {
      this->throwInvalidBehaviourType(
          throwing,
          "expected a finite strain behaviour or a general behaviour");
    }
    if (this->b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      if (this->b.kinematic != Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY) {
        this->throwInvalidBehaviourKinematic(
            throwing, "expected a standard finite strain behaviour kinematic");
      }
    }
    if ((this->b.gradients.size() != 1) ||
        (!isDeformationGradient(this->b.gradients.at(0)))) {
      this->throwInvalidBehaviourType(
          throwing,
          "expected a finite strain behaviour or a general behaviour whose "
          "only gradient is the deformation gradient, a non symmetric "
          "tensor");
    }
    if ((this->b.thermodynamic_forces.size() != 1) ||
        (!isPK1(this->b.thermodynamic_forces.at(0)))) {
      this->throwInvalidBehaviourType(
          throwing,
          "expected a finite strain behaviour or a general behaviour whose "
          "only thermodynamic force is the first Piola-Kirchhoff stress, a non "
          "symmetric tensor");
    }
    if ((this->b.to_blocks.size() != 1) ||
        (!isPK1(this->b.to_blocks.at(0).first)) ||
        (!isDeformationGradient(this->b.to_blocks.at(0).second))) {
      this->throwInvalidBehaviourType(
          throwing,
          "expected a finite strain behaviour or a general behaviour whose "
          "only tangent operator block is the derivative of the first "
          "Piola-Kirchhoff stress with respect to the deformation gradient, "
          "both being non symmetric tensors");
    }
  }  // end of checkIfAFiniteStrainBehaviourIsDeclared

  void BehaviourIntegratorBase::checkBehaviourVariablesSizes(
      attributes::Throwing,
      const CheckBehaviourVariablesSizesArguments& args) const {
    const auto h = this->b.hypothesis;
    auto get_size = [h](const mgis::behaviour::Variable& v) {
      return static_cast<size_type>(getVariableSize(v, h));
    };
    auto gsizes = std::vector<size_type>{};
    for (const auto& g : this->b.gradients) {
      gsizes.push_back(get_size(g));
    }
    if (gsizes != args.gradients_sizes) {
      this->throwInvalidBehaviourType(
          throwing, "the sizes of the gradients " +
                        convertSizesToString(gsizes) +
                        " do not match the expected ones " +
                        convertSizesToString(args.gradients_sizes));
    }
    auto fsizes = std::vector<size_type>{};
    for (const auto& f : this->b.thermodynamic_forces) {
      fsizes.push_back(get_size(f));
    }
    if (fsizes != args.thermodynamic_forces_sizes) {
      this->throwInvalidBehaviourType(
          throwing, "the sizes of the thermodynamic forces " +
                        convertSizesToString(fsizes) +
                        " do not match the expected ones " +
                        convertSizesToString(args.thermodynamic_forces_sizes));
    }
    auto bsizes = std::vector<std::pair<size_type, size_type>>{};
    for (const auto& [v1, v2] : this->b.to_blocks) {
      bsizes.emplace_back(get_size(v1), get_size(v2));
    }
    if (bsizes != args.tangent_operator_blocks_sizes) {
      this->throwInvalidBehaviourType(
          throwing, "the sizes of the tangent operator blocks " +
                        convertSizesOfTangentOperatorBlocksToString(bsizes) +
                        " do not match the expected ones " +
                        convertSizesOfTangentOperatorBlocksToString(
                            args.tangent_operator_blocks_sizes));
    }
  }  // end of checkBehaviourVariablesSizes

  void BehaviourIntegratorBase::checkBehaviourSymmetry(
      attributes::Throwing, const Behaviour::Symmetry s) const {
    if (this->b.symmetry != s) {
      if (s == Behaviour::ORTHOTROPIC) {
        this->throwInvalidBehaviourSymmetry(
            throwing, "expected an orthotropic behaviour");
      } else {
        this->throwInvalidBehaviourSymmetry(  //
            throwing, "expected an isotropic behaviour");
      }
    }
  }  // end of checkBehaviourSymmetry

  void BehaviourIntegratorBase::throwInvalidBehaviourType(
      attributes::Throwing, const std::string& e) const {
    auto msg = std::string{"invalid behaviour type"};
    if (!e.empty()) {
      msg += ", " + e;
    }
    raise(prependBehaviourDescriptionToErrorMessage(this->b, msg));
  }  // end of throwInvalidBehaviourType

  void BehaviourIntegratorBase::throwInvalidBehaviourSymmetry(
      attributes::Throwing, const std::string& e) const {
    auto msg = std::string{"invalid behaviour symmetry"};
    if (!e.empty()) {
      msg += ", " + e;
    }
    raise(prependBehaviourDescriptionToErrorMessage(this->b, msg));
  }  // end of throwInvalidBehaviourSymmetry

  void BehaviourIntegratorBase::throwInvalidBehaviourKinematic(
      attributes::Throwing, const std::string& e) const {
    auto msg = std::string{"invalid behaviour kinematic"};
    if (!e.empty()) {
      msg += ", " + e;
    }
    raise(prependBehaviourDescriptionToErrorMessage(this->b, msg));
  }  // end of throwInvalidBehaviourKinematic

  static bool checkQuadratureFunctionEvaluator(
      Context& ctx,
      const AbstractBehaviourIntegrator& b,
      const AbstractQPEvaluator& e) {
    const auto& bqspace = b.getPartialQuadratureSpace();
    const auto& eqspace = *(e.getPartialQuadratureSpacePointer());
    const auto& bfed = bqspace.getFiniteElementDiscretization();
    const auto& efed = eqspace.getFiniteElementDiscretization();
    if (static_cast<const MeshDiscretization&>(bfed) !=
        static_cast<const MeshDiscretization&>(efed)) {
      return ctx.registerErrorMessage(
          "partial quadrature function evaluator is not defined on the mesh of "
          "the behaviour integrator");
    }
    if (bqspace.getId() != eqspace.getId()) {
      return ctx.registerErrorMessage(
          "partial quadrature function is not defined on the same material as "
          "the behaviour integrator");
    }
    if (bqspace.getNumberOfIntegrationPoints() !=
        eqspace.getNumberOfIntegrationPoints()) {
      return ctx.registerErrorMessage(
          "the partial quadrature space of the partial quadrature function "
          "does not have the same number of quadrature points than the partial "
          "quadrature space of the the behaviour integrator");
    }
    return true;
  }  // end of checkQuadratureFunctionEvaluator

  bool BehaviourIntegratorBase::setMaterialProperty(
      Context& ctx,
      std::string_view name,
      std::shared_ptr<const AbstractQPEvaluator> e,
      const TimeStepStage ts) noexcept {
    const auto omp = getVariable(ctx, this->b.mps, name);
    if (isInvalid(omp)) {
      return false;
    }
    if (e.get() == nullptr) {
      return ctx.registerErrorMessage(
          "invalid partial quadrature function evaluator given for material "
          "property '" +
          omp->name + "'");
    }
    if (!checkQuadratureFunctionEvaluator(ctx, *this, *e)) {
      return false;
    }
    if (omp->type != mgis::behaviour::Variable::SCALAR) {
      return ctx.registerErrorMessage(
          "invalid material property '" + omp->name +
          "' (only scalar material properties are supported)");
    }
    if (e->getNumberOfComponents() != 1) {
      return ctx.registerErrorMessage(
          "invalid number of components for partial quadrature function "
          "evaluator given for material property '" +
          omp->name + "' (" + std::to_string(e->getNumberOfComponents()) +
          ", expected 1)");
    }
    if (ts == TimeStepStage::BEGINNING_OF_TIME_STEP) {
      if (!unsetMaterialProperty(ctx, this->s0, name)) {
        return false;
      }
      this->material_properties_evaluators_bts.insert_or_assign(
          std::string{name}, e);
    } else {
      if (!unsetMaterialProperty(ctx, this->s1, name)) {
        return false;
      }
      this->material_properties_evaluators_ets.insert_or_assign(
          std::string{name}, e);
    }
    return true;
  }  // end of setMaterialProperty

  bool BehaviourIntegratorBase::setExternalStateVariable(
      Context& ctx,
      std::string_view name,
      std::shared_ptr<const AbstractQPEvaluator> e,
      const TimeStepStage ts) noexcept {
    const auto oesv = getVariable(ctx, this->b.esvs, name);
    if (isInvalid(oesv)) {
      return false;
    }
    if (e.get() == nullptr) {
      return ctx.registerErrorMessage(
          "invalid partial quadrature function evaluator given for external "
          "state variable '" +
          oesv->name + "'");
    }
    if (!checkQuadratureFunctionEvaluator(ctx, *this, *e)) {
      return false;
    }
    const auto osize = getVariableSize(ctx, *oesv, this->b.hypothesis);
    if (isInvalid(osize)) {
      return false;
    }
    if (e->getNumberOfComponents() != *osize) {
      return ctx.registerErrorMessage(
          "invalid number of components for partial quadrature function "
          "evaluator given for external state variable '" +
          oesv->name + "' (" + std::to_string(e->getNumberOfComponents()) +
          ", expected " + std::to_string(*osize) + ")");
    }
    if (ts == TimeStepStage::BEGINNING_OF_TIME_STEP) {
      if (!unsetExternalStateVariable(ctx, this->s0, name)) {
        return false;
      }
      this->external_state_variables_evaluators_bts.insert_or_assign(
          std::string{name}, e);
    } else {
      if (!unsetExternalStateVariable(ctx, this->s1, name)) {
        return false;
      }
      this->external_state_variables_evaluators_ets.insert_or_assign(
          std::string{name}, e);
    }
    return true;
  }  // end of setExternalStateVariable

  const PartialQuadratureSpace&
  BehaviourIntegratorBase::getPartialQuadratureSpace() const noexcept {
    return static_cast<const Material*>(this)->getPartialQuadratureSpace();
  }  // end of getPartialQuadratureSpace

  real BehaviourIntegratorBase::getTimeIncrement() const noexcept {
    return this->time_increment;
  }  // end of getTimeIncrement

  void BehaviourIntegratorBase::setTimeIncrement(const real dt) {
    this->time_increment = dt;
  }  // end of setTimeIncrement

  bool BehaviourIntegratorBase::hasMaterial() const noexcept { return true; }

  OptionalReference<Material> BehaviourIntegratorBase::getMaterial(
      Context&) noexcept {
    return OptionalReference<Material>{this};
  }  // end of getMaterial

  OptionalReference<const Material> BehaviourIntegratorBase::getMaterial(
      Context&) const noexcept {
    return OptionalReference<const Material>{this};
  }  // end of getMaterial

  bool BehaviourIntegratorBase::setup(Context& ctx,
                                      const real t,
                                      const real dt) noexcept {
    // note: synchronization among MPI processes is done later
    // in the calling nonlinear evolution problem
    this->wks.pqfcts_mps_bts.clear();
    this->wks.pqfcts_mps_ets.clear();
    this->wks.pqfcts_esvs_bts.clear();
    this->wks.pqfcts_esvs_ets.clear();
    for (const auto& [name, e] : this->material_properties_evaluators_bts) {
      if (e->isUniform()) {
        const auto ovalue = e->getUniformValue(ctx, t, dt);
        if (isInvalid(ovalue)) {
          return false;
        }
        if (std::holds_alternative<real>(*ovalue)) {
          if (!mgis::behaviour::setMaterialProperty(
                  ctx, this->s0, name, std::get<real>(*ovalue),
                  mgis::behaviour::MaterialStateManager::NOUPDATE)) {
            return false;
          }
        } else {
          if (!mgis::behaviour::setMaterialProperty(
                  ctx, this->s0, name, std::get<std::vector<real>>(*ovalue),
                  mgis::behaviour::MaterialStateManager::LOCAL_STORAGE,
                  mgis::behaviour::MaterialStateManager::NOUPDATE)) {
            return false;
          }
        }
      } else {
        auto of = evaluate(ctx, *e, t, dt, {.ensure_contiguous_storage = true});
        if (isInvalid(of)) {
          return false;
        }
        if (this->wks.pqfcts_mps_bts.contains(name)) {
          this->wks.pqfcts_mps_bts.erase(name);
        }
        this->wks.pqfcts_mps_bts.insert({name, std::move(*of)});
      }
    }
    for (const auto& [name, e] : this->material_properties_evaluators_ets) {
      if (e->isUniform()) {
        const auto ovalue = e->getUniformValue(ctx, t, dt);
        if (isInvalid(ovalue)) {
          return false;
        }
        if (std::holds_alternative<real>(*ovalue)) {
          if (!mgis::behaviour::setMaterialProperty(
                  ctx, this->s1, name, std::get<real>(*ovalue),
                  mgis::behaviour::MaterialStateManager::NOUPDATE)) {
            return false;
          }
        } else {
          if (!mgis::behaviour::setMaterialProperty(
                  ctx, this->s1, name, std::get<std::vector<real>>(*ovalue),
                  mgis::behaviour::MaterialStateManager::LOCAL_STORAGE,
                  mgis::behaviour::MaterialStateManager::NOUPDATE)) {
            return false;
          }
        }
      } else {
        auto of = evaluate(ctx, *e, t, dt, {.ensure_contiguous_storage = true});
        if (isInvalid(of)) {
          return false;
        }
        if (this->wks.pqfcts_mps_ets.contains(name)) {
          this->wks.pqfcts_mps_ets.erase(name);
        }
        this->wks.pqfcts_mps_ets.insert({name, std::move(*of)});
      }
    }
    for (const auto& [name, e] :
         this->external_state_variables_evaluators_bts) {
      if (e->isUniform()) {
        const auto ovalue = e->getUniformValue(ctx, t, dt);
        if (isInvalid(ovalue)) {
          return false;
        }
        if (std::holds_alternative<real>(*ovalue)) {
          if (!mgis::behaviour::setExternalStateVariable(
                  ctx, this->s0, name, std::get<real>(*ovalue),
                  mgis::behaviour::MaterialStateManager::NOUPDATE)) {
            return false;
          }
        } else {
          if (!mgis::behaviour::setExternalStateVariable(
                  ctx, this->s0, name, std::get<std::vector<real>>(*ovalue),
                  mgis::behaviour::MaterialStateManager::LOCAL_STORAGE,
                  mgis::behaviour::MaterialStateManager::NOUPDATE)) {
            return false;
          }
        }
      } else {
        auto of = evaluate(ctx, *e, t, dt, {.ensure_contiguous_storage = true});
        if (isInvalid(of)) {
          return false;
        }
        if (this->wks.pqfcts_esvs_bts.contains(name)) {
          this->wks.pqfcts_esvs_bts.erase(name);
        }
        this->wks.pqfcts_esvs_bts.insert({name, std::move(*of)});
      }
    }
    for (const auto& [name, e] :
         this->external_state_variables_evaluators_ets) {
      if (e->isUniform()) {
        const auto ovalue = e->getUniformValue(ctx, t, dt);
        if (isInvalid(ovalue)) {
          return false;
        }
        if (std::holds_alternative<real>(*ovalue)) {
          if (!mgis::behaviour::setExternalStateVariable(
                  ctx, this->s1, name, std::get<real>(*ovalue),
                  mgis::behaviour::MaterialStateManager::NOUPDATE)) {
            return false;
          }
        } else {
          if (!mgis::behaviour::setExternalStateVariable(
                  ctx, this->s1, name, std::get<std::vector<real>>(*ovalue),
                  mgis::behaviour::MaterialStateManager::LOCAL_STORAGE,
                  mgis::behaviour::MaterialStateManager::NOUPDATE)) {
            return false;
          }
        }
      } else {
        auto of = evaluate(ctx, *e, t, dt, {.ensure_contiguous_storage = true});
        if (isInvalid(of)) {
          return false;
        }
        if (this->wks.pqfcts_esvs_ets.contains(name)) {
          this->wks.pqfcts_esvs_ets.erase(name);
        }
        this->wks.pqfcts_esvs_ets.insert({name, std::move(*of)});
      }
    }
    //
    for (auto& [name, e] : this->wks.pqfcts_mps_bts) {
      if (!mgis::behaviour::setMaterialProperty(
              ctx, this->s0, name, e.getValues(),
              mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
              mgis::behaviour::MaterialStateManager::NOUPDATE)) {
        return false;
      }
    }
    for (auto& [name, e] : this->wks.pqfcts_mps_ets) {
      if (!mgis::behaviour::setMaterialProperty(
              ctx, this->s1, name, e.getValues(),
              mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
              mgis::behaviour::MaterialStateManager::NOUPDATE)) {
        return false;
      }
    }
    for (auto& [name, e] : this->wks.pqfcts_esvs_bts) {
      if (!mgis::behaviour::setExternalStateVariable(
              ctx, this->s0, name, e.getValues(),
              mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
              mgis::behaviour::MaterialStateManager::NOUPDATE)) {
        return false;
      }
    }
    for (auto& [name, e] : this->wks.pqfcts_esvs_ets) {
      if (!mgis::behaviour::setExternalStateVariable(
              ctx, this->s1, name, e.getValues(),
              mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
              mgis::behaviour::MaterialStateManager::NOUPDATE)) {
        return false;
      }
    }
    /*
     * \brief uniform values are treated immediately. For spatially variable
     * fields, we return the information needed to evaluate them
     */
    // This lambda function builds an *evaluator*. Its role is to fill up
    // the `v` vector
    auto dispatch =
        [this, &ctx](
            std::vector<real>& v,
            const std::map<std::string,
                           mgis::behaviour::MaterialStateManager::FieldHolder,
                           std::less<>>& field_holders,
            const std::vector<mgis::behaviour::Variable>& ds)
        -> std::optional<
            std::vector<std::tuple<size_type, size_type, const real*>>> {
      const auto h = this->b.hypothesis;
      if (v.size() != mgis::behaviour::getArraySize(ds, h)) {
        return ctx.registerErrorMessage("integrate: ill allocated memory");
      }
      // evaluators
      auto evs = std::vector<std::tuple<size_type, size_type, const real*>>{};
      auto i = mgis::size_type{};
      for (const auto& d : ds) {
        const auto vsize = mgis::behaviour::getVariableSize(d, h);
        auto p = field_holders.find(d.name);
        if (p == field_holders.end()) {
          auto msg = std::string{"integrate: no variable named '" + d.name +
                                 "' declared"};
          if (!field_holders.empty()) {
            msg += "\nThe following variables were declared: ";
            for (const auto& vs : field_holders) {
              msg += "\n- " + vs.first;
            }
          } else {
            msg += "\nNo variable declared.";
          }
          return ctx.registerErrorMessage(msg);
        }
        // depending on the type of p->second, we are branching
        // on one of the following procedure:
        if (std::holds_alternative<std::monostate>(p->second)) {
          return ctx.registerErrorMessage(
              "field holder associated with variable '" + d.name +
              "' is not initialized");
        } else if (std::holds_alternative<std::span<const real>>(p->second)) {
          const auto& variable_values =
              std::get<std::span<const real>>(p->second);
          //  here we test if the variable_values is uniform: its size is then
          //  the size of targeted external state variable
          if (variable_values.size() == vsize) {
            std::copy(variable_values.begin(), variable_values.end(),
                      v.begin() + i);
          } else {
            //  here we test if the variable_values is not uniform: its size is
            //  then the size of targeted external state variable multiplied by
            //  the number of integration point.
            if (variable_values.size() != vsize * (this->n)) {
              return ctx.registerErrorMessage(
                  "invalid number of variable values given for variable '" +
                  d.name + "'");
            }
            // if we have a span, we store in evs this span for future use
            evs.push_back(std::make_tuple(i, vsize, variable_values.data()));
          }
        } else {
          using MutableFieldHolder =
              mgis::behaviour::MaterialStateManager::MutableFieldHolder;
          if (!std::holds_alternative<MutableFieldHolder>(p->second)) {
            return ctx.registerErrorMessage(
                "unsupport field holder for variable '" + d.name + "'");
          }
          const auto& field_value =
              std::get<MutableFieldHolder>(p->second).value;
          if (std::holds_alternative<real>(field_value)) {
            if (vsize != 1) {
              return ctx.registerErrorMessage(
                  "invalid number of values given for variable '" + d.name +
                  "'");
            }
            // if uniform field, copy field_value into v[i]
            // `evs` will be untouched.
            v[i] = std::get<real>(field_value);
          } else if (std::holds_alternative<std::span<real>>(field_value)) {
            const auto& variable_values =
                std::get<std::span<real>>(field_value);
            if (variable_values.size() == vsize) {
              std::copy(variable_values.begin(), variable_values.end(),
                        v.begin() + i);
            } else {
              if (variable_values.size() != vsize * (this->n)) {
                return ctx.registerErrorMessage(
                    "invalid number of variable values given for variable '" +
                    d.name + "'");
              }
              // if we have a span, we store in evs this span for future use
              evs.push_back(std::make_tuple(i, vsize, variable_values.data()));
            }
          } else {
            const auto& variable_values =
                std::get<std::vector<real>>(field_value);
            if (variable_values.size() == vsize) {
              std::copy(variable_values.begin(), variable_values.end(),
                        v.begin() + i);
            } else {
              if (variable_values.size() != vsize * (this->n)) {
                return ctx.registerErrorMessage(
                    "invalid number of variable values given for variable '" +
                    d.name + "'");
              }
              // if we have a vector, we store in evs this vector for future
              // use
              evs.push_back(std::make_tuple(i, vsize, variable_values.data()));
            }
          }
        }
        i += vsize;
      }
      return evs;
    };  // end of dispatch
    // The `b` field refers to MGIS MaterialDataManager class, which
    // is an ancestor class of the current BehaviourIntegratorBase class.
    auto oevaluators =
        dispatch(this->wks.mps, this->s1.material_properties, this->b.mps);
    if (isInvalid(oevaluators)) {
      return false;
    }
    this->wks.mps_evaluators = *oevaluators;
    // external state variables at the beginning of the time step
    oevaluators = dispatch(this->wks.esvs0, this->s0.external_state_variables,
                           this->b.esvs);
    if (isInvalid(oevaluators)) {
      return false;
    }
    this->wks.esvs0_evaluators = *oevaluators;
    // external state variables at the end of the time step
    oevaluators = dispatch(this->wks.esvs1, this->s1.external_state_variables,
                           this->b.esvs);
    if (isInvalid(oevaluators)) {
      return false;
    }
    this->wks.esvs1_evaluators = *oevaluators;
    return true;
  }  // end of setup

  bool BehaviourIntegratorBase::performsLocalBehaviourIntegration(
      const size_type ip, const IntegrationType it) {
    this->wks.has_rotated_thermodynamic_forces = false;
    char error_msg[512];
    error_msg[0] = '\0';
    const auto g_offset = this->s0.gradients_stride * ip;
    const auto t_offset = this->s0.thermodynamic_forces_stride * ip;
    const auto isvs_offset = this->s0.internal_state_variables_stride * ip;
    // Fill vector `v` with material properties or external state variables.
    // Evaluators are used to deal with both uniform or spatially
    // variable quantities.
    auto eval =
        [](std::vector<real>& v,
           const std::vector<std::tuple<size_type, size_type, const real*>>&
               evs,
           const size_type i) {
          for (const auto& [offset, size, values] : evs) {
            if (size == 1) {
              v[offset] = values[i];
            } else {
              const auto b_ptr = values + i * size;
              const auto e_ptr = values + (i + 1) * size;
              std::copy(b_ptr, e_ptr, v.begin() + offset);
            }
          }
        };  // end of eval
    eval(this->wks.mps, this->wks.mps_evaluators, ip);
    eval(this->wks.esvs0, this->wks.esvs0_evaluators, ip);
    eval(this->wks.esvs1, this->wks.esvs1_evaluators, ip);
    //
    rdt = real{1};
    mgis::behaviour::BehaviourDataView v;
    v.error_message = error_msg;
    v.rdt = &rdt;
    v.dt = this->time_increment;
    v.K = this->K.data() + this->K_stride * ip;
    v.speed_of_sound = nullptr;
    v.s0.gradients = this->s0.gradients.data() + g_offset;
    v.s1.gradients = this->s1.gradients.data() + g_offset;
    v.s0.thermodynamic_forces = this->s0.thermodynamic_forces.data() + t_offset;
    v.s1.thermodynamic_forces = this->s1.thermodynamic_forces.data() + t_offset;
    v.s0.material_properties = this->wks.mps.data();
    v.s1.material_properties = this->wks.mps.data();
    v.s0.internal_state_variables =
        this->s0.internal_state_variables.data() + isvs_offset;
    v.s1.internal_state_variables =
        this->s1.internal_state_variables.data() + isvs_offset;
    if (this->b.computesStoredEnergy) {
      v.s0.stored_energy = this->s0.stored_energies.data() + ip;
      v.s1.stored_energy = this->s1.stored_energies.data() + ip;
    } else {
      v.s0.stored_energy = nullptr;
      v.s1.stored_energy = nullptr;
    }
    if (this->b.computesDissipatedEnergy) {
      v.s0.dissipated_energy = this->s0.dissipated_energies.data() + ip;
      v.s1.dissipated_energy = this->s1.dissipated_energies.data() + ip;
    } else {
      v.s0.dissipated_energy = nullptr;
      v.s1.dissipated_energy = nullptr;
    }
    v.s0.mass_density = nullptr;
    v.s1.mass_density = nullptr;
    v.s0.external_state_variables = this->wks.esvs0.data();
    v.s1.external_state_variables = this->wks.esvs1.data();
    v.K[0] = static_cast<int>(it);
    const auto r = mgis::behaviour::integrate(v, this->b);
    if (!((r == 0) || (r == 1))) {
      std::cerr << "behaviour integration failed: " << error_msg
                << "(exit value: " << r << ")" << std::endl;
    }
    return (r == 0) || (r == 1);
  }  // end of BehaviourIntegratorBase::performsLocalBehaviourIntegration

  /*!
   * \return true if the computation of the inner forces needs the
   * thermodynamic forces in the global frame of all the points (see
   * `BehaviourIntegratorBase::getRotatedThermodynamicForces`)
   * \param[in] b: behaviour
   */
  static bool requiresRotatedThermodynamicForces(const Behaviour& b) noexcept {
    return (b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) &&
           (b.rotate_thermodynamic_forces_ptr == nullptr);
  }  // end of requiresRotatedThermodynamicForces

  /*!
   * \return a view on the values of a material property
   * \param[in] h: material property
   */
  static mgis::behaviour::MaterialStateManager::FieldHolder getView(
      const mgis::behaviour::MaterialStateManager::FieldHolder& h) {
    using FieldHolder = mgis::behaviour::MaterialStateManager::FieldHolder;
    using MutableFieldHolder =
        mgis::behaviour::MaterialStateManager::MutableFieldHolder;
    if (const auto* const mh = std::get_if<MutableFieldHolder>(&h)) {
      if (const auto* const v = std::get_if<real>(&(mh->value))) {
        auto r = FieldHolder{};
        r = *v;
        return r;
      }
      if (const auto* const v = std::get_if<std::span<real>>(&(mh->value))) {
        return std::span<const real>(*v);
      }
      return std::span<const real>(std::get<std::vector<real>>(mh->value));
    }
    return h;
  }  // end of getView

  namespace {

    /*!
     * \brief gives the behaviour, for the integration by batches, the
     * material properties at the end of the time step as the ones at the
     * beginning of the time step (views, without copy), and restores the
     * material properties at the beginning of the time step at its
     * destruction
     */
    struct MaterialPropertiesAtTheBeginningOfTheTimeStepGuard {
      /*!
       * \param[in, out] s0: state at the beginning of the time step
       * \param[in] s1: state at the end of the time step
       */
      MaterialPropertiesAtTheBeginningOfTheTimeStepGuard(
          mgis::behaviour::MaterialStateManager& s0,
          const mgis::behaviour::MaterialStateManager& s1)
          : state(s0), saved(std::move(s0.material_properties)) {
        s0.material_properties.clear();
        for (const auto& [name, h] : this->saved) {
          s0.material_properties[name] = getView(h);
        }
        for (const auto& [name, h] : s1.material_properties) {
          s0.material_properties[name] = getView(h);
        }
      }  // end of MaterialPropertiesAtTheBeginningOfTheTimeStepGuard
      //! \brief destructor
      ~MaterialPropertiesAtTheBeginningOfTheTimeStepGuard() {
        this->state.material_properties = std::move(this->saved);
      }

     private:
      //! \brief state at the beginning of the time step
      mgis::behaviour::MaterialStateManager& state;
      //! \brief saved material properties
      decltype(mgis::behaviour::MaterialStateManager::material_properties)
          saved;
    };

  }  // end of namespace

  bool BehaviourIntegratorBase::integratesBatchesOnDevice() const noexcept {
    return mfem::Device::Allows(mfem::Backend::DEVICE_MASK) &&
           mgis::behaviour::isGenericParallelBehaviour(this->b) &&
           mgis::behaviour::isMemoryLocationAccessible(
               this->b, getDeviceMemoryLocation());
  }  // end of integratesBatchesOnDevice

  /*!
   * \return a view of a copy, readable on the device, of the values of a
   * material property, of an external state variable or of the mass
   * density, given by a view of their values on the host
   * \param[out] c: copy, stored point by point (one point for uniform
   * values)
   * \param[in] v: view of the values on the host
   */
  static mgis::FieldView<const real> copyToDevice(
      mfem::Vector& c, const mgis::FieldView<const real>& v) {
    const auto np = v.isUniform() ? size_type{1} : v.getNumberOfPoints();
    const auto nc = v.getNumberOfComponents();
    const auto size = static_cast<std::size_t>(np) * nc;
    if (size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      raise("BehaviourIntegratorBase: array too large for MFEM");
    }
    c.UseDevice(true);
    c.SetSize(static_cast<int>(size));
    auto* const values = c.HostWrite();
    const auto* const h = v.data();
    const auto ps = v.getPointStride();
    const auto cs = v.getComponentStride();
    for (size_type i = 0; i != np; ++i) {
      for (size_type j = 0; j != nc; ++j) {
        values[i * nc + j] = h[i * ps + j * cs];
      }
    }
    const auto d = std::span<const real>(c.Read(), size);
    if (v.isUniform()) {
      return mgis::makeUniformFieldView(d, v.getNumberOfPoints());
    }
    return mgis::makePointByPointFieldView(d, np, nc, nc);
  }  // end of copyToDevice

  bool BehaviourIntegratorBase::integrateBatchOnDevice(
      const IntegrationType it) {
    using StateArray = mfem_mgis::MaterialArrays::StateArray;
    const auto orthotropic =
        this->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC;
    const auto prediction =
        (it == IntegrationType::PREDICTION_TANGENT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_SECANT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_ELASTIC_OPERATOR);
    const auto tangent = it != IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR;
    const auto l = getDeviceMemoryLocation();
    this->wks.has_rotated_thermodynamic_forces = false;
    const auto opts = mgis::behaviour::BehaviourIntegrationOptions{
        .integration_type =
            static_cast<mgis::behaviour::IntegrationType>(static_cast<int>(it)),
        .compute_speed_of_sound = false};
    auto r = mgis::behaviour::BehaviourIntegrationResult{};
    try {
      auto& arrays = this->getArrays();
      // the material properties at the end of the time step are given as the
      // ones at the beginning of the time step (see `integrateBatch`)
      const auto guard = MaterialPropertiesAtTheBeginningOfTheTimeStepGuard(
          this->s0, this->s1);
      // view of the arrays on the host, whose views are replaced by views of
      // the copies on the device
      auto v = mgis::behaviour::make_unchecked_view(*this);
      v.memory_location = l;
      // arrays read by the behaviour (`Read`), which are not modified by it,
      // or read and written (`ReadWrite`)
      const auto relocateArray =
          [&arrays](mgis::FieldView<real>& view, const StateArray a,
                    const TimeStepStage s, const bool output) {
            const auto& array = arrays.getArray(a, s);
            auto* const p = output ? arrays.readWrite(a, s)
                                   : const_cast<real*>(array.Read());
            view = relocate(view, p, array.Size());
          };
      for (const auto& [view, a] :
           {std::pair{&v.s0.gradients, StateArray::GRADIENTS},
            std::pair{&v.s0.thermodynamic_forces,
                      StateArray::THERMODYNAMIC_FORCES},
            std::pair{&v.s0.internal_state_variables,
                      StateArray::INTERNAL_STATE_VARIABLES},
            std::pair{&v.s0.stored_energies, StateArray::STORED_ENERGIES},
            std::pair{&v.s0.dissipated_energies,
                      StateArray::DISSIPATED_ENERGIES}}) {
        relocateArray(*view, a, bts, false);
      }
      // the gradients at the end of the time step of an orthotropic behaviour
      // are rotated in place, as on the host
      relocateArray(v.s1.gradients, StateArray::GRADIENTS, ets,
                    orthotropic && !prediction);
      for (const auto& [view, a] :
           {std::pair{&v.s1.thermodynamic_forces,
                      StateArray::THERMODYNAMIC_FORCES},
            std::pair{&v.s1.internal_state_variables,
                      StateArray::INTERNAL_STATE_VARIABLES},
            std::pair{&v.s1.stored_energies, StateArray::STORED_ENERGIES},
            std::pair{&v.s1.dissipated_energies,
                      StateArray::DISSIPATED_ENERGIES}}) {
        relocateArray(*view, a, ets, true);
      }
      if (tangent) {
        v.K = relocate(v.K, arrays.writeTangentOperatorBlocks(),
                       arrays.getTangentOperatorBlocks().Size());
      } else {
        v.K = {};
      }
      v.speed_of_sound = {};
      // material properties, external state variables and mass density
      // copied to the device
      auto inputs = std::vector<mgis::FieldView<const real>*>{};
      for (auto* const st : {&v.s0, &v.s1}) {
        for (auto& mp : st->material_properties) {
          inputs.push_back(&mp);
        }
        for (auto& esv : st->external_state_variables) {
          inputs.push_back(&esv);
        }
        inputs.push_back(&(st->mass_density));
      }
      auto& copies = this->wks.device_inputs;
      copies.resize(inputs.size());
      for (std::size_t i = 0; i != inputs.size(); ++i) {
        if (!inputs[i]->empty()) {
          *(inputs[i]) = copyToDevice(copies[i], *(inputs[i]));
        }
      }
      // rotation matrices, computed on the host and copied to the device
      auto matrices = mgis::FieldView<const real>{};
      if (orthotropic) {
        // `computeRotationMatrices` resizes the array before its copy
        auto& rm = this->wks.rotation_matrices;
        const auto hmatrices = this->computeRotationMatrices(rm);
        matrices = relocate(hmatrices, rm.Read(), rm.Size());
      }
      if (orthotropic && !prediction) {
        mgis::behaviour::rotateGradients(v.s1.gradients, this->b, matrices, l);
      }
      r = mgis::behaviour::integrate(v.batch_execution_policy, v, opts,
                                     this->time_increment);
      if (!((r.exit_status == 0) || (r.exit_status == 1))) {
        std::cerr << "behaviour integration failed: " << r.error_message
                  << " (exit value: " << r.exit_status
                  << ", first point in failure: " << r.n << ")" << std::endl;
        return false;
      }
      if (orthotropic && tangent) {
        mgis::behaviour::rotateTangentOperatorBlocks(v.K, this->b, matrices, l);
      }
      if (requiresRotatedThermodynamicForces(this->b)) {
        const auto np = static_cast<size_type>(this->n);
        const auto nc = this->s1.thermodynamic_forces_stride;
        const auto size = static_cast<std::size_t>(np) * nc;
        if (size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
          raise("BehaviourIntegratorBase: array too large for MFEM");
        }
        auto& rf = this->wks.rotated_thermodynamic_forces;
        rf.UseDevice(true);
        rf.SetSize(static_cast<int>(size));
        const auto d = mgis::makePointByPointFieldView(
            std::span<real>(rf.Write(), size), np, nc, nc);
        mgis::behaviour::rotateThermodynamicForces(
            d, this->b, v.s1.thermodynamic_forces, matrices, l);
        this->wks.has_rotated_thermodynamic_forces = true;
      }
    } catch (std::exception& e) {
      std::cerr << "behaviour integration failed: " << e.what() << std::endl;
      return false;
    }
    return true;
  }  // end of integrateBatchOnDevice

  bool BehaviourIntegratorBase::integrateBatch(const IntegrationType it) {
    if (!this->supportsBatchIntegration()) {
      return false;
    }
    if (this->integratesBatchesOnDevice()) {
      return this->integrateBatchOnDevice(it);
    }
    const auto orthotropic =
        this->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC;
    const auto prediction =
        (it == IntegrationType::PREDICTION_TANGENT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_SECANT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_ELASTIC_OPERATOR);
    const auto rotateTangentOperatorBlocks =
        orthotropic && (it != IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR);
    this->wks.has_rotated_thermodynamic_forces = false;
    const auto opts = mgis::behaviour::BehaviourIntegrationOptions{
        .integration_type =
            static_cast<mgis::behaviour::IntegrationType>(static_cast<int>(it)),
        .compute_speed_of_sound = false};
    auto r = mgis::behaviour::BehaviourIntegrationResult{};
    try {
      auto matrices = mgis::FieldView<const real>{};
      if (orthotropic) {
        matrices = this->computeRotationMatrices(this->wks.rotation_matrices);
      }
      if (orthotropic && !prediction) {
        const auto g = mgis::behaviour::getGradientsView(this->s1);
        mgis::behaviour::rotateGradients(g, this->b, matrices);
      }
      {
        const auto guard = MaterialPropertiesAtTheBeginningOfTheTimeStepGuard(
            this->s0, this->s1);
        r = mgis::behaviour::integrate(*this, opts, this->time_increment);
      }
      if (!((r.exit_status == 0) || (r.exit_status == 1))) {
        std::cerr << "behaviour integration failed: " << r.error_message
                  << " (exit value: " << r.exit_status
                  << ", first point in failure: " << r.n << ")" << std::endl;
        return false;
      }
      if (rotateTangentOperatorBlocks) {
        const auto blocks =
            mgis::behaviour::getTangentOperatorBlocksView(*this);
        mgis::behaviour::rotateTangentOperatorBlocks(blocks, this->b, matrices);
      }
      if (requiresRotatedThermodynamicForces(this->b)) {
        this->computeRotatedThermodynamicForces(matrices);
      }
    } catch (std::exception& e) {
      std::cerr << "behaviour integration failed: " << e.what() << std::endl;
      return false;
    }
    return true;
  }  // end of integrateBatch

  void BehaviourIntegratorBase::computeRotatedThermodynamicForces(
      const mgis::FieldView<const real>& matrices) {
    const auto np = static_cast<size_type>(this->n);
    const auto nc = this->s1.thermodynamic_forces_stride;
    const auto size = static_cast<std::size_t>(np) * nc;
    if (size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      raise("BehaviourIntegratorBase: array too large for MFEM");
    }
    auto& v = this->wks.rotated_thermodynamic_forces;
    v.SetSize(static_cast<int>(size));
    if (size != 0) {
      const auto d = mgis::makePointByPointFieldView(
          std::span<real>(v.HostWrite(), size), np, nc, nc);
      mgis::behaviour::rotateThermodynamicForces(
          d, this->b, mgis::behaviour::getThermodynamicForcesView(this->s1),
          matrices);
    }
    this->wks.has_rotated_thermodynamic_forces = true;
  }  // end of computeRotatedThermodynamicForces

  bool BehaviourIntegratorBase::updateRotatedThermodynamicForces() {
    if (!requiresRotatedThermodynamicForces(this->b)) {
      return false;
    }
    if (!this->wks.has_rotated_thermodynamic_forces) {
      // on the host, or on the device if the computations of the behaviour
      // can't access the memory of the host
      this->computeThermodynamicForcesInGlobalFrame(
          this->wks.rotated_thermodynamic_forces);
      this->wks.has_rotated_thermodynamic_forces = true;
    }
    return true;
  }  // end of updateRotatedThermodynamicForces

  std::span<const real>
  BehaviourIntegratorBase::getRotatedThermodynamicForces() {
    if (!this->updateRotatedThermodynamicForces()) {
      return {};
    }
    const auto& v = this->wks.rotated_thermodynamic_forces;
    if (v.Size() == 0) {
      return {};
    }
    return {v.HostRead(), static_cast<std::size_t>(v.Size())};
  }  // end of getRotatedThermodynamicForces

  bool BehaviourIntegratorBase::cleanup(Context& ctx) noexcept {
    for (const auto& [name, e] : this->material_properties_evaluators_bts) {
      std::ignore = e;
      if (!mgis::behaviour::unsetMaterialProperty(ctx, this->s0, name)) {
        return false;
      }
    }
    for (const auto& [name, e] : this->material_properties_evaluators_ets) {
      std::ignore = e;
      if (!mgis::behaviour::unsetMaterialProperty(ctx, this->s1, name)) {
        return false;
      }
    }
    for (const auto& [name, e] :
         this->external_state_variables_evaluators_bts) {
      std::ignore = e;
      if (!mgis::behaviour::unsetExternalStateVariable(ctx, this->s0, name)) {
        return false;
      }
    }
    for (const auto& [name, e] :
         this->external_state_variables_evaluators_ets) {
      std::ignore = e;
      if (!mgis::behaviour::unsetExternalStateVariable(ctx, this->s1, name)) {
        return false;
      }
    }
    this->wks.pqfcts_mps_bts.clear();
    this->wks.pqfcts_mps_ets.clear();
    this->wks.pqfcts_esvs_bts.clear();
    this->wks.pqfcts_esvs_ets.clear();
    return true;
  }  // end of cleanup

  void BehaviourIntegratorBase::copyStatesOnDevice(const bool u) {
    using StateArray = mfem_mgis::MaterialArrays::StateArray;
    auto& arrays = this->getArrays();
    // view of the arrays on the host, whose views of the states and of the
    // tangent operator blocks are replaced by views of the copies on the
    // device; the material properties and the external state variables are
    // not used
    auto v = mgis::behaviour::make_unchecked_view(*this);
    v.memory_location = getDeviceMemoryLocation();
    for (auto& [st, stage] : {std::pair{&v.s0, bts}, std::pair{&v.s1, ets}}) {
      // the state at the end of the time step is copied in the one at the
      // beginning of the time step by an update, and conversely
      const auto output = (stage == bts) == u;
      for (const auto& [view, a] :
           {std::pair{&(st->gradients), StateArray::GRADIENTS},
            std::pair{&(st->thermodynamic_forces),
                      StateArray::THERMODYNAMIC_FORCES},
            std::pair{&(st->internal_state_variables),
                      StateArray::INTERNAL_STATE_VARIABLES},
            std::pair{&(st->stored_energies), StateArray::STORED_ENERGIES},
            std::pair{&(st->dissipated_energies),
                      StateArray::DISSIPATED_ENERGIES}}) {
        const auto& array = arrays.getArray(a, stage);
        auto* const p =
            output ? arrays.write(a, stage) : const_cast<real*>(array.Read());
        *view = relocate(*view, p, array.Size());
      }
      std::fill(st->material_properties.begin(), st->material_properties.end(),
                mgis::FieldView<const real>{});
      std::fill(st->external_state_variables.begin(),
                st->external_state_variables.end(),
                mgis::FieldView<const real>{});
      st->mass_density = {};
    }
    // the tangent operator blocks are set to zero
    v.K = relocate(v.K, arrays.writeTangentOperatorBlocks(),
                   arrays.getTangentOperatorBlocks().Size());
    v.speed_of_sound = {};
    auto ops = batch::MFEMFieldOperations(true);
    if (u) {
      mgis::behaviour::update(v, ops);
    } else {
      mgis::behaviour::revert(v, ops);
    }
  }  // end of copyStatesOnDevice

  void BehaviourIntegratorBase::revert() {
    this->wks.has_rotated_thermodynamic_forces = false;
    if (this->integratesBatchesOnDevice()) {
      this->copyStatesOnDevice(false);
    } else {
      this->revertMaterialStates();
    }
  }  // end of revert

  void BehaviourIntegratorBase::update() {
    if (this->integratesBatchesOnDevice()) {
      this->copyStatesOnDevice(true);
    } else {
      this->updateMaterialStates();
    }
  }  // end of update

  void BehaviourIntegratorBase::setMacroscopicGradients(
      std::span<const real> g) {
    Material::setMacroscopicGradients(g);
  }  // end of setMacroscopicGradients

  bool BehaviourIntegratorBase::requiresCurrentSolutionForResidualAssembly()
      const noexcept {
    return false;
  }  // end of requiresCurrentSolutionForResidualAssembly

  bool BehaviourIntegratorBase::requiresCurrentSolutionForJacobianAssembly()
      const noexcept {
    return false;
  }  // end of requiresCurrentSolutionForJacobianAssembly

  BehaviourIntegratorBase::~BehaviourIntegratorBase() = default;

}  // end of namespace mfem_mgis
