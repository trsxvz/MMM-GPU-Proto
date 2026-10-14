/*!
 * \file   UniaxialMacroscopicStressPeriodicSimulation.cxx
 * \brief
 * \author Thomas Helfer
 * \date   16/09/2025
 */

#include <cmath>
#include <limits>
#include <algorithm>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MM_OPERA_HPC/Utilities.hxx"
#include "MM_OPERA_HPC/UniaxialMacroscopicStressPeriodicSimulation.hxx"

namespace mm_opera_hpc {

  [[nodiscard]] static const std::array<mfem_mgis::real, 3u>
  computeMacroscopicCauchyStress(
      mfem_mgis::attributes::MayAbort,
      mfem_mgis::Context &ctx,
      mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
      const std::array<mfem_mgis::real, 3u> &F) {
    auto or_die = ctx.getFatalFailureHandler();
    // integrals of the diagonal components of the First Piola-Kirchhoff stress
    auto pk1_integral = std::array<mfem_mgis::real, 3u>{};
    auto volume = mfem_mgis::real{};
    // get pk1 integral and volume on all materials //
    const auto [pk1_integrals, volumes] =
        mfem_mgis::computeMeanThermodynamicForcesValues<true>(
            ctx, problem.template getImplementation<true>()) |
        or_die;
    // sum volume of all materials
    // sum pk1 integrals
    for (const auto &m : problem.getAssignedMaterialsIdentifiers()) {
      volume += volumes[m];
      pk1_integral[0] += pk1_integrals[m][0];
      pk1_integral[1] += pk1_integrals[m][1];
      pk1_integral[2] += pk1_integrals[m][2];
    }
    assert(volume > 0.0);
    // -- sum volumes of all processess
    MPI_Allreduce(MPI_IN_PLACE, &volume, 1,
                  mfem_mgis::mpi_type<mfem_mgis::real>, MPI_SUM,
                  MPI_COMM_WORLD);
    // -- sum integrals from all processes
    MPI_Allreduce(MPI_IN_PLACE, pk1_integral.data(), 3,
                  mfem_mgis::mpi_type<mfem_mgis::real>, MPI_SUM,
                  MPI_COMM_WORLD);
    // macroscopic values of the First Piola-Kirchhoff stress
    const auto pk1 = std::array<mfem_mgis::real, 3u>{pk1_integral[0] / volume,
                                                     pk1_integral[1] / volume,
                                                     pk1_integral[2] / volume};
    // macroscopic change of volume
    const auto det = F[0] * F[1] * F[2];
    // computation of the diagonal components of the Cauchy stress
    return {pk1[0] * F[0] / det, pk1[1] * F[1] / det, pk1[2] * F[2] / det};
  }  // end of computeMacroscopicCauchyStress

  /*!
   * \return the smallest residual on the lateral stresses that the fixed
   * point can reduce. When both lateral stresses are close to S, the
   * correction of the lateral components of the deformation gradient is
   * about (1 + nu) (1 - 2 nu) S / E: it is rounded to zero when it is smaller
   * than half the spacing of the floating point numbers around these
   * components, and the fixed point then stagnates or oscillates.
   * Negligible with double. With float, about 1.5e4 Pa for the UO2 based
   * cases in tension (lateral components just below 1), twice more in
   * compression.
   * \param[in] mp: macroscopic elastic material properties
   * \param[in] F: diagonal components of the deformation gradient
   */
  [[nodiscard]] static mfem_mgis::real computeLateralStressResolution(
      const MacroscropicElasticMaterialProperties &mp,
      const std::array<mfem_mgis::real, 3u> &F) noexcept {
    constexpr auto infinity = std::numeric_limits<mfem_mgis::real>::infinity();
    const auto spacing = [](const mfem_mgis::real x) {
      const auto ax = std::abs(x);
      return std::nextafter(ax, infinity) - ax;
    };
    const auto c = (1 + mp.nueff) * (1 - 2 * mp.nueff);
    if (!(c > 0)) {
      return mfem_mgis::real{0};
    }
    const auto ds = std::max(spacing(F[0]), spacing(F[1]));
    return std::sqrt(mfem_mgis::real{2}) * (ds / 2) * mp.Eeff / c;
  }  // end of computeLateralStressResolution

  [[nodiscard]] static bool simulateOverATimeStep(
      mgis::Context &ctx,
      std::ostream &out,
      mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
      UniaxialMacroscopicStressPeriodicSimulation::MacroscopicUnknowns
          &macroscopic_unknowns,  //
      const std::function<mfem_mgis::real(mfem_mgis::real)>
          &imposed_axial_deformation_gradient_value,
      const UniaxialMacroscopicStressPeriodicSimulation::NumericalParameters
          &np,  //
      const bool post_processing,
      const mfem_mgis::real bts,  //
      const mfem_mgis::real ets) {
    using namespace mfem_mgis::Profiler::Utils;
    auto or_die = ctx.getFatalFailureHandler();
    Message("Solving time step from ", bts, " to ", ets);
    // --- fixed-point param
    const auto tolFP = np.macroscopic_stress_absolute_tolerance;
    const int maxitFP = np.macroscopic_stress_maximum_number_of_iterations;
    // Traction
    const auto dt = ets - bts;

    auto &F = macroscopic_unknowns.F;
    auto &S = macroscopic_unknowns.S;
    // deformation gradient at the beginning of the time step
    auto S0 = S;
    auto F0 = F;

    if (macroscopic_unknowns.previous_time_increment.has_value()) {
      auto &dF = macroscopic_unknowns.dF;
      const auto pdt = *(macroscopic_unknowns.previous_time_increment);
      // -- update F from input data, here def = 5e-4. -- //
      F[0] += dF[0] * (dt / pdt);
      F[1] += dF[1] * (dt / pdt);
    }
    F[2] = imposed_axial_deformation_gradient_value(ets);
    // --- fixed point iteration number
    int itFP = 0;
    // --- fixed point loop --- //
    bool converged = false;
    auto initial_residual = mfem_mgis::real{};
    while ((!converged) && (itFP < maxitFP)) {
      const auto &mp = np.macroscopic_elastic_material_properties;
      // -- compute correction -- //
      // -- note that S[0] and S[1] are updated at the end of this loop section,
      // default is 0 -- //
      F[0] += ((mp.nueff * mp.nueff - 1) * S[0]) / mp.Eeff +
              (mp.nueff * (mp.nueff + 1) * S[1]) / mp.Eeff;
      F[1] += ((mp.nueff * mp.nueff - 1) * S[1]) / mp.Eeff +
              (mp.nueff * (mp.nueff + 1) * S[0]) / mp.Eeff;
      // -- note that MacroscopicGradientsEvolution takes into account of new
      // F[0] and F[1] -- //
      if (itFP == 0) {
        if (!problem.setSolverParameters(
                ctx,
                {{"VerbosityLevel", np.microscopic_equilibrium_verbosity_level},
                 {"MaximumNumberOfIterations",
                  np.microscopic_equilibrium_maximum_number_of_iterations},
                 {"RelativeTolerance",
                  np.microscopic_equilibrium_relative_tolerance},
                 {"AbsoluteTolerance",
                  np.microscopic_equilibrium_absolute_tolerance}})) {
          Message("error: setSolverParameters failed.");
          return false;
        }
      } else {
        if (np.switch_microscopic_equilibrium_convergence_criterion) {
          const auto atol = std::max(
              initial_residual * np.microscopic_equilibrium_relative_tolerance,
              np.microscopic_equilibrium_absolute_tolerance);
          if (!problem.setSolverParameters(
                  ctx,
                  {{"VerbosityLevel",
                    np.microscopic_equilibrium_verbosity_level},
                   {"MaximumNumberOfIterations",
                    np.microscopic_equilibrium_maximum_number_of_iterations},
                   {"AbsoluteTolerance", atol}})) {
            Message("error: setSolverParameters failed.");
            return false;
          }
        }
      }
      auto statistics = problem.solve(ctx, bts, dt);
      if (!statistics.status) {
        Message("error: the resolution failed.");
        S = S0;
        F = F0;
        problem.revert(ctx) | or_die;
        return false;
      }
      if (itFP == 0) {
        initial_residual = statistics.initial_residual_norm;
      }
      S = computeMacroscopicCauchyStress(mfem_mgis::may_abort, ctx, problem, F);
      const auto r = std::sqrt(S[0] * S[0] + S[1] * S[1]);
      //
      Message("Fixed Point iteration", itFP, ": |res| =", r);
      // the tolerance is at least twice the resolution of the fixed point
      const auto tolerance =
          std::max(tolFP, 2 * computeLateralStressResolution(mp, F));
      converged = std::abs(r) < tolerance;
      if (converged && (std::abs(r) >= tolFP)) {
        Message("Fixed Point converged within the resolution of the",
                "deformation gradient: tolerance raised to", tolerance);
      }
      itFP++;
    }
    if (!converged) {
      Message(
          "error: maximum number of iterations for the "
          "fixed-point algorithm reached");
      S = S0;
      F = F0;
      problem.revert(ctx) | or_die;
      return false;
    }
    // the values of a rejected time step are not written. The macroscopic
    // values are the same on all the processes
    if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) {
      out << ets << " "                                 //
          << F[0] << " " << F[1] << " " << F[2] << " "  //
          << S[0] << " " << S[1] << " " << S[2] << std::endl;
    }
    Message("Solution at time", ets, ":", F[2], S[0], S[1], S[2]);
    //
    if (post_processing) {
      // add postprocessing
      CatchTimeSection(ctx, "common::post_processing_step");
      problem.executePostProcessings(ctx, bts, dt);
    }
    // update state variable for the next time step
    problem.update(ctx) | or_die;
    // update information for time extrapolation
    auto &dF = macroscopic_unknowns.dF;
    macroscopic_unknowns.previous_time_increment = dt;
    for (std::size_t i = 0; i != 3; ++i) {
      dF[i] = F[i] - F0[i];
    }
    return true;
  }

  [[nodiscard]] static bool simulateOverATemporalSequence(
      mgis::Context &ctx,
      std::ostream &out,
      mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
      UniaxialMacroscopicStressPeriodicSimulation::MacroscopicUnknowns
          &macroscopic_unknowns,  //
      const std::function<mfem_mgis::real(mfem_mgis::real)>
          &imposed_axial_deformation_gradient_value,
      const UniaxialMacroscopicStressPeriodicSimulation::NumericalParameters
          &np,  //
      const bool post_processing,
      const mfem_mgis::real bts,  //
      const mfem_mgis::real ets) {
    using namespace mfem_mgis::Profiler::Utils;
    constexpr auto maximum_of_substeps = std::size_t{10};
    auto nsubsteps = std::size_t{};
    auto nstep = std::size_t{1};
    auto t = bts;
    auto dt = ets - bts;
    while (nstep != 0) {
      const auto do_post_processing = post_processing && (nstep == 1);
      const auto success =
          simulateOverATimeStep(ctx, out, problem, macroscopic_unknowns,
                                imposed_axial_deformation_gradient_value, np,
                                do_post_processing, t, t + dt);
      if (success) {
        t += dt;
        --nstep;
      } else {
        ++nsubsteps;
        if (nsubsteps == maximum_of_substeps) {
          Message(
              "Maximum number of substeps reached for temporal sequence from",
              bts, "to", ets);
          return false;
        } else {
          Message("dividing time step by 2");
        }
        dt *= mfem_mgis::real{0.5};
        nstep *= 2;
      }
    }
    //
    printMemoryFootprint(
        "[At end of temporal sequence: " + std::to_string(ets) + "]");
    return true;
  }  // end of simulateOverATemporalSequence

  UniaxialMacroscopicStressPeriodicSimulation::
      UniaxialMacroscopicStressPeriodicSimulation(
          mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
          const std::function<mfem_mgis::real(mfem_mgis::real)> Fzz,
          const NumericalParameters &np,
          const bool b) noexcept
      : problem(p),
        imposed_axial_deformation_gradient_value(Fzz),
        numerical_parameters(np),
        post_processing(b) {
    // warning Fxx, F[1], and Fzz are defined during the following loop
    problem.setMacroscopicGradientsEvolution([this](const mfem_mgis::real) {
      const auto &F = this->macroscopic_unknowns.F;
      auto ret = std::vector<mfem_mgis::real>(9, mfem_mgis::real{});
      std::copy(F.begin(), F.end(), ret.begin());
      return ret;
    });
  }  // end of UniaxialMacroscopicStressPeriodicSimulation

  bool UniaxialMacroscopicStressPeriodicSimulation::run(
      mgis::Context &ctx,
      std::ostream &out,
      const std::vector<mfem_mgis::real> &temporal_sequences) noexcept {
    for (std::size_t i = 0; i + 1 != temporal_sequences.size(); i++) {
      if (!simulateOverATemporalSequence(
              ctx, out, this->problem, this->macroscopic_unknowns,
              this->imposed_axial_deformation_gradient_value,
              this->numerical_parameters, this->post_processing,
              temporal_sequences[i], temporal_sequences[i + 1])) {
        return false;
      }
    }
    return true;
  }

}  // end of namespace mm_opera_hpc
