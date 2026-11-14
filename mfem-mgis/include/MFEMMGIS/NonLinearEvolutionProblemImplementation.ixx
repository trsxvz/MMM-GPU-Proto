/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemImplementation.ixx
 * \brief  This file implements the `computeResultantForceOnBoundary` and
 * `computeMeanThermodynamicForcesValues` functions
 * \author Thomas Helfer
 * \date   28/03/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX
#define LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX

#include <map>
#include <span>
#include <exception>
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis::internals {

  /*!
   * \return if two materials have the same thermodynamic forces
   * \param[in] m1: first material
   * \param[in] m2: second material
   */
  inline bool haveSameThermodynamicForces(const Material& m1,
                                          const Material& m2) noexcept {
    const auto& tfs1 = m1.b.thermodynamic_forces;
    const auto& tfs2 = m2.b.thermodynamic_forces;
    if (tfs1.size() != tfs2.size()) {
      return false;
    }
    for (std::size_t i = 0; i != tfs1.size(); ++i) {
      if ((tfs1[i].name != tfs2[i].name) || (tfs1[i].type != tfs2[i].type)) {
        return false;
      }
    }
    return m1.s1.thermodynamic_forces_stride ==
           m2.s1.thermodynamic_forces_stride;
  }  // end of haveSameThermodynamicForces

  /*!
   * \brief apply a function to the materials of the selected behaviour
   * integrators of all the materials, so that all the processes report the
   * same errors
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear evolution problem
   * \param[in] selection: selection of the behaviour integrators
   * \param[in] f: function called with the execution context and a
   * material, returning true on success
   */
  template <bool parallel, typename Function>
  bool applyToSelectedMaterials(
      Context& ctx,
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const BehaviourIntegratorsSelection& selection,
      const Function& f) noexcept {
    for (const auto& m : p.getAssignedMaterialsIdentifiers()) {
      const auto obis = getSelectedBehaviourIntegrators(ctx, p, m, selection);
      if (isInvalid(obis)) {
        return false;
      }
      for (auto b = obis->first; b != obis->second; ++b) {
        const auto obi = p.getBehaviourIntegrator(ctx, m, b);
        if (isInvalid(obi)) {
          return ctx.registerErrorMessage("invalid behaviour integrator");
        }
        if (!obi->hasMaterial()) {
          continue;
        }
        const auto om = obi->getMaterial(ctx);
        if (isInvalid(om)) {
          return false;
        }
        if (!f(ctx, *om)) {
          return false;
        }
      }
    }
    return true;
  }  // end of applyToSelectedMaterials

  /*!
   * \brief synchronize on the host, for reading, the thermodynamic forces at
   * the end of the time step of the materials of the selected behaviour
   * integrators (see `Material::synchronizeOnHostForReading`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear evolution problem
   * \param[in] selection: selection of the behaviour integrators
   */
  template <bool parallel>
  bool synchronizeSelectedMaterialsOnHostForReading(
      Context& ctx,
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const BehaviourIntegratorsSelection& selection) noexcept {
    return applyToSelectedMaterials(
        ctx, p, selection, [](Context& c, const Material& m) noexcept -> bool {
          try {
            m.synchronizeOnHostForReading(
                MaterialArrays::StateArray::THERMODYNAMIC_FORCES,
                TimeStepStage::END_OF_TIME_STEP);
          } catch (...) {
            return registerExceptionInErrorBacktrace(c);
          }
          return true;
        });
  }  // end of synchronizeSelectedMaterialsOnHostForReading

}  // end of namespace mfem_mgis::internals

namespace mfem_mgis {

  template <bool parallel>
  bool computeResultantForceOnBoundary(
      Context& ctx,
      mfem::Vector& F,
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const std::vector<
          std::pair<size_type, std::vector<std::vector<size_type>>>>& elements,
      const BehaviourIntegratorsSelection& selection) noexcept {
    auto fed = p.getFiniteElementDiscretization();
    const auto& fes = fed.template getFiniteElementSpace<parallel>();
    const auto nc = fes.GetVDim();
    mfem::Vector elt_forces;
    F.SetSize(nc);
    F = real{0};
    // the inner forces are computed element by element, from the
    // thermodynamic forces read on the host
    if (!internals::synchronizeSelectedMaterialsOnHostForReading(ctx, p,
                                                                 selection)) {
      return false;
    }
    for (const auto& e : elements) {
      const auto& fe = *(fes.GetFE(e.first));
      auto& tr = *(fes.GetElementTransformation(e.first));
      const auto nnodes = fe.GetDof();
      // inner forces of the selected behaviour integrators
      const auto obis =
          getSelectedBehaviourIntegrators(ctx, p, tr.Attribute, selection);
      if (isInvalid(obis)) {
        return false;
      }
      for (auto b = obis->first; b != obis->second; ++b) {
        const auto obi = p.getBehaviourIntegrator(ctx, tr.Attribute, b);
        if (isInvalid(obi)) {
          return ctx.registerErrorMessage("invalid behaviour integrator");
        }
        try {
          obi->computeInnerForces(elt_forces, fe, tr);
        } catch (std::exception& ex) {
          return ctx.registerErrorMessage(ex.what());
        } catch (...) {
          return ctx.registerErrorMessage("unknown exception");
        }
        for (size_type c = 0; c != nc; ++c) {
          const auto* const Fe = elt_forces.GetData() + c * nnodes;
          for (const auto& i : e.second[c]) {
            F[c] += Fe[i];
          }
        }
      }
    }
    return true;
  }  // end of computeResultantForceOnBoundary

  template <bool parallel>
  std::optional<std::pair<std::vector<std::vector<real>>, std::vector<real>>>
  computeMeanThermodynamicForcesValues(
      Context& ctx,
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const BehaviourIntegratorsSelection& selection) noexcept {
    // the thermodynamic forces are read on the host
    if (!internals::synchronizeSelectedMaterialsOnHostForReading(ctx, p,
                                                                 selection)) {
      return {};
    }
    auto nmax = p.getFiniteElementSpace().GetMesh()->attributes.Max() + 1;
    std::vector<std::vector<mfem_mgis::real>> stress_integrals(nmax);
    const auto& mis = p.getAssignedMaterialsIdentifiers();
    for (const auto& mi : mis) {
      const auto obis = getSelectedBehaviourIntegrators(ctx, p, mi, selection);
      if (isInvalid(obis)) {
        return {};
      }
      const auto [first, last] = *obis;
      const auto obi = p.getBehaviourIntegrator(ctx, mi, first);
      if (isInvalid(obi)) {
        return {};
      }
      const auto om = obi->getMaterial(ctx);
      if (isInvalid(om)) {
        return {};
      }
      // the thermodynamic forces of the selected behaviour integrators are
      // summed: they must be the same
      for (auto b = first + 1; b != last; ++b) {
        const auto obi2 = p.getBehaviourIntegrator(ctx, mi, b);
        if (isInvalid(obi2)) {
          return {};
        }
        const auto om2 = obi2->getMaterial(ctx);
        if (isInvalid(om2)) {
          return {};
        }
        if (!internals::haveSameThermodynamicForces(*om, *om2)) {
          return ctx.registerErrorMessage(
              "computeMeanThermodynamicForcesValues: the selected behaviour "
              "integrators of material '" +
              std::to_string(mi) +
              "' do not have the same thermodynamic forces");
        }
      }
      const auto& s1 = om->s1;
      const auto thsize = s1.thermodynamic_forces_stride;
      stress_integrals[mi].resize(thsize, mfem_mgis::real(0));
    }
    std::vector<mfem_mgis::real> volumes(stress_integrals.size(),
                                         mfem_mgis::real(0));
    // The points are summed material by material, with the weights of the
    // points computed once (`getIntegrationPointsWeights`), in the order of
    // the points, which is the order of the elements: the sums are the ones
    // of a loop on the elements of each material.
    // Thermodynamic forces in the global frame of the orthotropic behaviours,
    // rotated at all the points at once
    auto rotated_thermodynamic_forces = mfem::Vector{};
    for (const auto& mi : mis) {
      const auto obis = getSelectedBehaviourIntegrators(ctx, p, mi, selection);
      if (isInvalid(obis)) {
        return {};
      }
      auto& s = stress_integrals[mi];
      auto& v = volumes[mi];
      for (auto b = obis->first; b != obis->second; ++b) {
        const auto obi = p.getBehaviourIntegrator(ctx, mi, b);
        if (isInvalid(obi)) {
          return {};
        }
        const auto& om = obi->getMaterial(ctx);
        if (isInvalid(om)) {
          return {};
        }
        const auto ow = p.getIntegrationPointsWeights(ctx, mi, b);
        if (isInvalid(ow)) {
          return {};
        }
        const auto& w = *ow;
        const auto& s1 = om->s1;
        const auto& qspace = om->getPartialQuadratureSpace();
        const auto thsize =
            static_cast<mfem_mgis::size_type>(s1.thermodynamic_forces_stride);
        const auto npoints = qspace.getNumberOfIntegrationPoints();
        if (static_cast<mfem_mgis::size_type>(w.size()) != npoints) {
          return ctx.registerErrorMessage(
              "computeMeanThermodynamicForcesValues: inconsistent number of "
              "weights of the integration points of material '" +
              std::to_string(mi) + "'");
        }
        // thermodynamic forces, read in the layout of the arrays of the
        // material, or point by point in the global frame for an orthotropic
        // behaviour
        auto thforces = mgis::behaviour::getThermodynamicForcesView(s1);
        if (om->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {
          try {
            om->computeThermodynamicForcesInGlobalFrame(
                rotated_thermodynamic_forces);
          } catch (std::exception& ex) {
            return ctx.registerErrorMessage(ex.what());
          } catch (...) {
            return ctx.registerErrorMessage("unknown exception");
          }
          thforces = mgis::makePointByPointFieldView(
              std::span<const real>(rotated_thermodynamic_forces.HostRead(),
                                    static_cast<std::size_t>(
                                        rotated_thermodynamic_forces.Size())),
              npoints, thsize, thsize);
        }
        for (mfem_mgis::size_type o = 0; o != npoints; ++o) {
          const auto wo = w[static_cast<std::size_t>(o)];
          for (mfem_mgis::size_type k = 0; k != thsize; ++k) {
            s[k] += wo * thforces(o, k);
          }
          // the volume is computed with the first selected behaviour integrator
          if (b == obis->first) {
            v += wo;
          }
        }
      }
    }
    return std::pair<std::vector<std::vector<real>>, std::vector<real>>{
        std::move(stress_integrals), std::move(volumes)};
  }  // end of computeMeanThermodynamicForcesValues

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX */
