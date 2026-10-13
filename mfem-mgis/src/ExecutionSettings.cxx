/*!
 * \file   src/ExecutionSettings.cxx
 * \brief  This file implements the choice of the execution of the
 * computations on the CPU, on the GPU, or of each phase on one of them
 * \date   13/10/2026
 */

#include <set>
#include <array>
#include <string>
#include <cstring>
#include <optional>
#include <string_view>
#include "mfem/config/config.hpp"
#include "mfem/general/device.hpp"
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/BehaviourDescription.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/ExecutionSettings.hxx"
#include "ExecutionSettings.hxx"

namespace mfem_mgis {

  //! \return the execution settings, modifiable by `readExecutionSettings`
  static ExecutionSettings& getMutableExecutionSettings() noexcept {
    static ExecutionSettings s;
    return s;
  }  // end of getMutableExecutionSettings

  const ExecutionSettings& getExecutionSettings() noexcept {
    return getMutableExecutionSettings();
  }  // end of getExecutionSettings

  std::string_view getExecutionLocationName(
      const ExecutionLocation l) noexcept {
    return l == ExecutionLocation::GPU ? "gpu" : "cpu";
  }  // end of getExecutionLocationName

  std::string getBehavioursLibrary(const std::string_view d,
                                   const std::string_view n) {
#if defined(__APPLE__)
    constexpr auto suffix = std::string_view{".dylib"};
#elif defined(_WIN32)
    constexpr auto suffix = std::string_view{".dll"};
#else
    constexpr auto suffix = std::string_view{".so"};
#endif
    const auto& s = getExecutionSettings();
    const auto gpu = (s.mode != ExecutionMode::UNSPECIFIED) &&
                     (s.behaviours == ExecutionLocation::GPU);
    auto l = std::string{d};
    l += gpu ? "/gpu/lib" : "/cpu/lib";
    l += n;
    l += suffix;
    return l;
  }  // end of getBehavioursLibrary

}  // end of namespace mfem_mgis

namespace mfem_mgis::internals {

  //! \brief option placing a phase of the computations
  struct LocationOption {
    //! \brief long name
    const char* name;
    //! \brief short name
    const char* short_name;
    //! \brief location, modified if the option is given
    ExecutionLocation ExecutionSettings::* location;
  };

  //! \brief options placing the phases of the computations
  static constexpr auto location_options = std::array<LocationOption, 4>{
      LocationOption{finite_elements_option, finite_elements_short_option,
                     &ExecutionSettings::finite_elements},
      LocationOption{behaviours_option, behaviours_short_option,
                     &ExecutionSettings::behaviours},
      LocationOption{linear_solver_device_option,
                     linear_solver_device_short_option,
                     &ExecutionSettings::linear_solver},
      LocationOption{post_processings_option, post_processings_short_option,
                     &ExecutionSettings::post_processings}};

  /*!
   * \return the value of an option, if the current argument is this option
   * \param[in, out] a: current argument, moved to the value of the option
   * \param[in] e: end of the arguments
   * \param[in] name: long name of the option
   * \param[in] short_name: short name of the option
   * \param[in] previous: value already given, if any
   */
  static const char* readOptionValue(MainFunctionArguments& a,
                                     const MainFunctionArguments e,
                                     const char* const name,
                                     const char* const short_name,
                                     const char* const previous) {
    if ((std::strcmp(*a, name) != 0) && (std::strcmp(*a, short_name) != 0)) {
      return previous;
    }
    if (previous != nullptr) {
      mgis::raise("initialize: option '" + std::string{name} + "' given twice");
    }
    ++a;
    if ((a == e) || (**a == '-')) {
      mgis::raise("initialize: no value given to '" + std::string{name} + "'");
    }
    return *a;
  }  // end of readOptionValue

  /*!
   * \return the location given to an option
   * \param[in] name: long name of the option
   * \param[in] v: value
   */
  static ExecutionLocation getLocation(const char* const name,
                                       const std::string_view v) {
    if (v == "cpu") {
      return ExecutionLocation::CPU;
    }
    if (v == "gpu") {
      return ExecutionLocation::GPU;
    }
    mgis::raise("initialize: invalid value '" + std::string{v} +
                "' for option '" + std::string{name} +
                "' (expected 'cpu' or 'gpu')");
  }  // end of getLocation

  void readExecutionSettings(const int argc, const MainFunctionArguments argv) {
    const char* mode = nullptr;
    auto locations = std::array<const char*, location_options.size()>{};
    const auto e = argv + argc;
    for (auto a = argv; a != e; ++a) {
      const auto* const previous = mode;
      mode =
          readOptionValue(a, e, execution_option, execution_short_option, mode);
      if (mode != previous) {
        continue;
      }
      for (std::size_t i = 0; i != location_options.size(); ++i) {
        const auto& o = location_options[i];
        const auto* const p = locations[i];
        locations[i] = readOptionValue(a, e, o.name, o.short_name, p);
        if (locations[i] != p) {
          break;
        }
      }
    }
    auto s = ExecutionSettings{};
    if (mode != nullptr) {
      const auto m = std::string_view{mode};
      if (m == "cpu") {
        s.mode = ExecutionMode::CPU;
      } else if (m == "gpu") {
        s.mode = ExecutionMode::GPU;
      } else if (m == "mixed") {
        s.mode = ExecutionMode::MIXED;
      } else {
        mgis::raise("initialize: invalid value '" + std::string{m} +
                    "' for option '" + std::string{execution_option} +
                    "' (expected 'cpu', 'gpu' or 'mixed')");
      }
    }
    // the locations of the phases are only given in the mixed mode: in the
    // other modes, all the computations are at the same place
    for (std::size_t i = 0; i != location_options.size(); ++i) {
      if ((locations[i] != nullptr) && (s.mode != ExecutionMode::MIXED)) {
        mgis::raise("initialize: option '" +
                    std::string{location_options[i].name} + "' requires '" +
                    std::string{execution_option} + " mixed'");
      }
    }
#if !defined(MFEM_USE_CUDA) && !defined(MFEM_USE_HIP)
    if ((s.mode == ExecutionMode::GPU) || (s.mode == ExecutionMode::MIXED)) {
      mgis::raise("initialize: MFEM is built without CUDA nor HIP, '" +
                  std::string{execution_option} + " " + std::string{mode} +
                  "' requires a GPU");
    }
#endif /* !defined(MFEM_USE_CUDA) && !defined(MFEM_USE_HIP) */
    const auto d = s.mode == ExecutionMode::CPU ? ExecutionLocation::CPU
                                                : ExecutionLocation::GPU;
    for (std::size_t i = 0; i != location_options.size(); ++i) {
      const auto& o = location_options[i];
      s.*(o.location) =
          locations[i] != nullptr ? getLocation(o.name, locations[i]) : d;
    }
#if defined(MFEM_USE_HIP) && !defined(MFEM_MGIS_HAVE_DEVICE_SOURCES)
    // the kernels of mfem-mgis are not compiled with HIP
    if ((s.finite_elements == ExecutionLocation::GPU) ||
        (s.post_processings == ExecutionLocation::GPU)) {
      mgis::raise(
          "initialize: the kernels of mfem-mgis are not compiled with HIP: "
          "the finite elements and the post-processings can't be computed on "
          "the GPU (use '--execution mixed --finite-elements cpu "
          "--post-processings cpu')");
    }
#endif /* defined(MFEM_USE_HIP) && !defined(MFEM_MGIS_HAVE_DEVICE_SOURCES) */
    getMutableExecutionSettings() = s;
  }  // end of readExecutionSettings

  bool isExecutionModeSpecified() noexcept {
    return getExecutionSettings().mode != ExecutionMode::UNSPECIFIED;
  }  // end of isExecutionModeSpecified

  const char* getPlacementEnvironmentVariable(const char* const n) noexcept {
    const auto* const e = std::getenv(n);
    if ((e == nullptr) || (!isExecutionModeSpecified())) {
      return e;
    }
    // warning printed once per variable, by the main process
    static auto warned = std::set<std::string>{};
    if ((warned.insert(n).second) &&
        (isMPIMainProcessForDefaultCommunicator())) {
      getErrorStream() << "warning: the environment variable " << n
                       << " is ignored, the computations being placed by "
                          "the option '"
                       << execution_option << "'\n";
    }
    return nullptr;
  }  // end of getPlacementEnvironmentVariable

  std::string_view getDefaultBehaviourIntegration() noexcept {
    // PETSc does not support the integration by batches
    return (isExecutionModeSpecified() && (!usePETSc())) ? "Batch"
                                                         : "PointByPoint";
  }  // end of getDefaultBehaviourIntegration

  std::string_view getDefaultFiniteElementsComputation() noexcept {
    const auto& s = getExecutionSettings();
    // PETSc does not support the computations by kernels
    if ((!isExecutionModeSpecified()) || (usePETSc())) {
      return "ElementByElement";
    }
    // with `--execution cpu`, the kernels compute on the host, faster than
    // the computations element by element (optional, see
    // `areFiniteElementsKernelsOptional`)
    if ((s.mode == ExecutionMode::CPU) ||
        (s.finite_elements == ExecutionLocation::GPU)) {
      return "Kernels";
    }
    return "ElementByElement";
  }  // end of getDefaultFiniteElementsComputation

  bool areFiniteElementsKernelsOptional() noexcept {
    return (getExecutionSettings().mode == ExecutionMode::CPU) && (!usePETSc());
  }  // end of areFiniteElementsKernelsOptional

  std::string_view getDefaultPostProcessingsComputation() noexcept {
    const auto& s = getExecutionSettings();
    if ((isExecutionModeSpecified()) &&
        (s.post_processings == ExecutionLocation::GPU)) {
      return "Kernels";
    }
    return "ElementByElement";
  }  // end of getDefaultPostProcessingsComputation

  /*!
   * \return the options to give to select a location of the behaviours or
   * of the linear solvers, for the error messages
   * \param[in] o: long name of the option placing the phase
   * \param[in] l: location
   */
  static std::string getLocationOptions(const char* const o,
                                        const ExecutionLocation l) {
    return "'" + std::string{execution_option} + " " +
           std::string{getExecutionLocationName(l)} + "' or '" +
           std::string{execution_option} + " mixed " + std::string{o} + " " +
           std::string{getExecutionLocationName(l)} + "'";
  }  // end of getLocationOptions

  bool checkBehaviourExecutionLocation(Context& ctx,
                                       const Behaviour& b) noexcept {
    if (!isExecutionModeSpecified()) {
      return true;
    }
    const auto requested = getExecutionSettings().behaviours;
    // the behaviours integrated on a device which is not a GPU (debug device
    // of MFEM) are integrated on the CPU; a library generated by the
    // generic-parallel interface for a device can't compute on the host
    const auto on_device = isIntegratedByBatchesOnDevice(b) &&
                           mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                                                mfem::Backend::HIP_MASK);
    const auto host_only = (!on_device) &&
                           mgis::behaviour::isGenericParallelBehaviour(b) &&
                           (!mgis::behaviour::isMemoryLocationAccessible(
                               b, mgis::MemoryLocation::HOST));
    const auto location =
        on_device ? ExecutionLocation::GPU : ExecutionLocation::CPU;
    if ((location == requested) && (!host_only)) {
      return true;
    }
    auto msg = "the behaviour '" + b.behaviour + "' of the library '" +
               b.library + "' ";
    if (host_only) {
      msg += "computes on a device, which is not enabled";
    } else {
      msg += "is integrated on the " +
             std::string{location == ExecutionLocation::GPU ? "GPU" : "CPU"};
    }
    msg += ", while the execution settings place the behaviours on the ";
    if (requested == ExecutionLocation::GPU) {
      msg +=
          "GPU: a library generated by the generic-parallel interface of "
          "MFront with a GPU backend (cuda or hip) is expected (or use " +
          getLocationOptions(behaviours_option, ExecutionLocation::CPU) + ")";
    } else {
      msg +=
          "CPU: a library computing on the host (generic interface, or "
          "generic-parallel interface with a host backend) is expected (or "
          "use " +
          getLocationOptions(behaviours_option, ExecutionLocation::GPU) + ")";
    }
    return ctx.registerErrorMessage(msg);
  }  // end of checkBehaviourExecutionLocation

  /*!
   * \brief where a linear solver or a preconditioner computes
   */
  enum struct LinearSolverPlacement {
    //! \brief on the host only (direct solvers of the CPU)
    HOST,
    //! \brief on the GPU only (AmgX, cuDSS)
    DEVICE,
    /*!
     * \brief where the device of MFEM computes (solvers and smoothers of
     * MFEM)
     */
    MFEM,
    //! \brief where the memory of hypre is (solvers of hypre)
    HYPRE
  };

  /*!
   * \return where a linear solver or a preconditioner computes
   * \param[in] n: name of the linear solver or of the preconditioner
   */
  static LinearSolverPlacement getLinearSolverPlacement(
      const std::string_view n) noexcept {
    if ((n == "MUMPSSolver") || (n == "UMFPackSolver")) {
      return LinearSolverPlacement::HOST;
    }
    if ((n == "AmgX") || (n == "AmgXPCG") || (n == "CuDSSSolver")) {
      return LinearSolverPlacement::DEVICE;
    }
    // preconditioners of hypre which only compute on the host
    if ((n == "HypreEuclid") || (n == "HypreParaSails")) {
      return LinearSolverPlacement::HOST;
    }
    if (n.starts_with("Hypre")) {
      return LinearSolverPlacement::HYPRE;
    }
    return LinearSolverPlacement::MFEM;
  }  // end of getLinearSolverPlacement

  bool checkLinearSolverExecutionLocation(Context& ctx,
                                          const std::string_view n) noexcept {
    const auto& s = getExecutionSettings();
    if (!isExecutionModeSpecified()) {
      return true;
    }
    const auto p = getLinearSolverPlacement(n);
    const auto requested = s.linear_solver;
    // computations of MFEM on the GPU in the gpu and mixed modes
    const auto mfem_location = s.mode == ExecutionMode::CPU
                                   ? ExecutionLocation::CPU
                                   : ExecutionLocation::GPU;
    auto location = requested;
    if (p == LinearSolverPlacement::HOST) {
      location = ExecutionLocation::CPU;
    } else if (p == LinearSolverPlacement::DEVICE) {
      location = ExecutionLocation::GPU;
    } else if (p == LinearSolverPlacement::MFEM) {
      location = mfem_location;
    } else if (p == LinearSolverPlacement::HYPRE) {
      // where the memory of hypre is: on the host if hypre is built without
      // GPU support
      location =
          isHypreUsingGPU() ? ExecutionLocation::GPU : ExecutionLocation::CPU;
    }
    if (location == requested) {
      return true;
    }
    auto msg = "the linear solver or preconditioner '" + std::string{n} +
               "' computes on the " +
               std::string{location == ExecutionLocation::GPU ? "GPU" : "CPU"};
    if (p == LinearSolverPlacement::MFEM) {
      msg += " (as the device of MFEM)";
    } else if (p == LinearSolverPlacement::HYPRE) {
      msg +=
          " (as the memory of hypre, on the host if hypre is built without "
          "GPU support)";
    }
    msg += ", while the execution settings place the linear solvers on the " +
           std::string{requested == ExecutionLocation::GPU ? "GPU" : "CPU"};
    if (p == LinearSolverPlacement::MFEM) {
      msg +=
          ": use a solver of hypre or a direct solver (if the non linear "
          "solver accepts it: the fixed point solvers require the solvers "
          "of MFEM), or place the linear solvers on the GPU";
    } else {
      msg += " (use " +
             getLocationOptions(linear_solver_device_option, location) + ")";
    }
    return ctx.registerErrorMessage(msg);
  }  // end of checkLinearSolverExecutionLocation

}  // end of namespace mfem_mgis::internals
