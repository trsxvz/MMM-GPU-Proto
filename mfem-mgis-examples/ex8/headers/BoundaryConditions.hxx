#pragma once
#include <functional>
#include <vector>
#include <memory>
#include <utility>
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/ModelBase.hxx"
#include "MFEMMGIS/TimeStep.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/UniformHeatSourceBoundaryCondition.hxx"
#include "MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx"

#include "Setup.hxx"
#include "RobinBC.hxx"

void apply_boundary_conditions(
    mfem_mgis::attributes::MayAbort,
    mfem_mgis::Context&,
    mfem_mgis::NonLinearEvolutionProblem& heat_transfer,
    mfem_mgis::NonLinearEvolutionProblem& mechanics,
    const TestParameters& p,
    const std::function<double(double)>& power_history,
    mfem::GridFunction* u_mech);

// Model used to manually update the source term before passing it to the
// constitutive model responsible for swelling, the displacement field used
// by the Robin boundary condition, and the copies on the host of the
// gradients of the mechanical materials, read by the heat transfer
// behaviours (see Setup.hxx).
class FieldUpdaterModel : public mfem_mgis::ModelBase {
 private:
  std::shared_ptr<std::vector<mfem_mgis::real>> Pow_s0;
  std::shared_ptr<std::vector<mfem_mgis::real>> Pow_s1;
  std::function<double(double)> power_history;
  mfem::GridFunction* u_mech;
  const mfem::Vector* u_mech_true;
  // mechanical materials whose gradients are read on the host by the heat
  // transfer behaviours as external state variables
  std::vector<const mfem_mgis::Material*> mechanical_materials;
  std::string name;

 public:
  FieldUpdaterModel(mfem_mgis::Context& ctx,
                    const mfem_mgis::MeshDiscretization& mesh,
                    std::shared_ptr<std::vector<mfem_mgis::real>> pow0,
                    std::shared_ptr<std::vector<mfem_mgis::real>> pow1,
                    std::function<double(double)> history_func,
                    mfem::GridFunction* u,
                    const mfem::Vector* u_true,
                    std::vector<const mfem_mgis::Material*> materials)
      : mfem_mgis::ModelBase(ctx, mesh),
        Pow_s0(pow0),
        Pow_s1(pow1),
        power_history(history_func),
        u_mech(u),
        u_mech_true(u_true),
        mechanical_materials(std::move(materials)),
        name("FieldUpdater") {}

  [[nodiscard]] std::string getName() const noexcept override { return name; }

  std::pair<mfem_mgis::ExitStatus,
            std::optional<mfem_mgis::ComputeNextStateOutput>>
  computeNextState(mfem_mgis::Context& ctx,
                   const mfem_mgis::TimeStep& ts) noexcept override {
    double p0 = power_history(ts.begin);
    double p1 = power_history(ts.end);

    for (auto& val : *Pow_s0) val = p0;
    for (auto& val : *Pow_s1) val = p1;
    u_mech->SetFromTrueDofs(*u_mech_true);
    // the mechanical behaviours may be integrated on a device, without
    // synchronization of their arrays on the host (see
    // mfem_mgis::Material::synchronizeOnHost)
    for (const auto* const m : mechanical_materials) {
      for (const auto s : {mfem_mgis::bts, mfem_mgis::ets}) {
        m->synchronizeOnHostForReading(
            mfem_mgis::MaterialArrays::StateArray::GRADIENTS, s);
      }
    }

    return {mfem_mgis::ExitStatus::success,
            mfem_mgis::ComputeNextStateOutput{}};
  }

  [[nodiscard]] bool update(mfem_mgis::Context&) noexcept override {
    return true;
  }
  [[nodiscard]] bool revert(mfem_mgis::Context&) noexcept override {
    return true;
  }
};
