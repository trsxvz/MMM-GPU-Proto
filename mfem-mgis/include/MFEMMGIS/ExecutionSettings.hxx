/*!
 * \file   include/MFEMMGIS/ExecutionSettings.hxx
 * \brief  This file declares the choice of the execution of the
 * computations on the CPU, on the GPU, or of each phase of the computations
 * on one of them
 * \date   13/10/2026
 */

#ifndef LIB_MFEM_MGIS_EXECUTIONSETTINGS_HXX
#define LIB_MFEM_MGIS_EXECUTIONSETTINGS_HXX

#include <string>
#include <string_view>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  //! \brief processor computing a phase of the computations
  enum struct ExecutionLocation { CPU, GPU };

  //! \brief execution of the computations, given by the option `--execution`
  enum struct ExecutionMode {
    /*!
     * \brief no choice: the computations are placed by the other options
     * (`--device`, `--hypre-on-host`), the parameters of the problems and
     * the environment variables, as without this option
     */
    UNSPECIFIED,
    /*!
     * \brief no GPU: all the computations on the host (a device of MFEM
     * computing on the host, such as `debug`, may be given by `--device`)
     */
    CPU,
    /*!
     * \brief device `gpu` of MFEM by default (or any configuration using a
     * GPU given by `--device`): all the computations on the GPU
     */
    GPU,
    /*!
     * \brief device `gpu` of MFEM by default (or any configuration using a
     * GPU given by `--device`), each phase of the computations being placed
     * on the CPU or on the GPU
     */
    MIXED
  };

  /*!
   * \brief choice of the processor computing each phase of the computations
   *
   * The settings are given by the options declared by
   * `declareDefaultOptions` and read by `initialize`:
   *
   * - `--execution` (`-x`): `cpu`, `gpu` or `mixed`;
   * - in the `mixed` mode only, the location of each phase, `cpu` or `gpu`
   *   (`gpu` by default): `--finite-elements` (`-fe`), `--behaviours`
   *   (`-bh`), `--linear-solver-device` (`-lsd`) and `--post-processings`
   *   (`-ppd`).
   *
   * The execution mode gives the default integration of the behaviours
   * (`Batch`), and the locations the default computations of the finite
   * elements (`Kernels` on the GPU, and on the host with `--execution cpu`)
   * and of the post-processings (`Kernels` on the GPU), the configuration
   * of the device of MFEM and the memory of hypre. With PETSc, the defaults
   * of the integration of the behaviours and of the computations of the
   * finite elements are `PointByPoint` and `ElementByElement`, and the
   * linear solvers can't be placed on the GPU. The locations are only
   * meaningful if the execution mode is given (`mode` is not
   * `UNSPECIFIED`). A parameter given to a problem takes precedence; the
   * environment variables placing the computations (`MFEM_MGIS_DEVICE`,
   * `MFEM_MGIS_HYPRE_ON_HOST`, `MFEM_MGIS_BEHAVIOUR_INTEGRATION`,
   * `MFEM_MGIS_*_COMPUTATION`) are ignored, with a warning. The integration
   * by batches and, with `--execution cpu`, the computations of the finite
   * elements by kernels on the host are replaced by the integration point by
   * point and the computations element by element if the behaviour
   * integrators do not support them. The behaviours and the linear solvers
   * are checked against these locations.
   */
  struct ExecutionSettings {
    //! \brief execution mode
    ExecutionMode mode = ExecutionMode::UNSPECIFIED;
    /*!
     * \brief location of the computations of the finite elements: data of
     * MFEM, gradients, inner forces and tangent operator
     */
    ExecutionLocation finite_elements = ExecutionLocation::CPU;
    //! \brief location of the integration of the behaviours
    ExecutionLocation behaviours = ExecutionLocation::CPU;
    //! \brief location of the linear solvers (and of the memory of hypre)
    ExecutionLocation linear_solver = ExecutionLocation::CPU;
    //! \brief location of the post-processings
    ExecutionLocation post_processings = ExecutionLocation::CPU;
  };

  /*!
   * \return the execution settings, read by `initialize`. The mode is
   * `UNSPECIFIED` before the call to `initialize` or if the option
   * `--execution` is not given.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const ExecutionSettings&
  getExecutionSettings() noexcept;
  /*!
   * \return the name of a location, as given to the options: `cpu` or
   * `gpu`
   * \param[in] l: location
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::string_view getExecutionLocationName(
      const ExecutionLocation l) noexcept;
  /*!
   * \return the path of a behaviour library built by the CMake function
   * `mfem_mgis_add_behaviours_library` (module `MFEMMGISBehaviours.cmake`,
   * installed with mfem-mgis): `<d>/gpu/lib<n>.so` if the execution
   * settings place the behaviours on the GPU, `<d>/cpu/lib<n>.so` otherwise
   * (in particular if the execution mode is not given). The suffix is the
   * one of the shared libraries of the system.
   * \param[in] d: directory of the libraries (`<target>_BEHAVIOURS_DIRECTORY`
   * in CMake)
   * \param[in] n: name of the library (`Behaviour` by default, as in
   * `mfem_mgis_add_behaviours_library`)
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::string getBehavioursLibrary(
      std::string_view d, std::string_view n = "Behaviour");

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_EXECUTIONSETTINGS_HXX */
