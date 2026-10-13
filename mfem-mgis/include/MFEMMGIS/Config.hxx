/*!
 * \file   include/MFEMMGIS/Config.hxx
 * \brief  This file declares the basic types and functions of `mfem-mgis`
 * \author Thomas Helfer
 * \date   19/06/2018
 */

#ifndef LIB_MFEM_MGIS_CONFIG_HXX
#define LIB_MFEM_MGIS_CONFIG_HXX

#include <limits>
#include <cstdlib>
#include "MGIS/Config.hxx"
#include "MGIS/Raise.hxx"
#include "MGIS/Context.hxx"
#include "MGIS/LogStream.hxx"
#include "MGIS/InvalidResult.hxx"
#include "MGIS/Utilities/Construct.hxx"
#include "MGIS/Utilities/Invoke.hxx"
#include "MGIS/Utilities/OptionalReference.hxx"

#include "MFEMMGIS/MGISForward.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

#ifndef MGIS_BEHAVIOUR_API_VERSION
#error "Incompatible version of MGIS"
#endif /* MGIS_BEHAVIOUR_API_VERSION */

#if MGIS_BEHAVIOUR_API_VERSION != 1
#error "Incompatible version of MGIS"
#endif /* MGIS_BEHAVIOUR_API_VERSION */

//! \brief macro hiding a symbol from the shared library
#define MFEM_MGIS_VISIBILITY_LOCAL MGIS_VISIBILITY_LOCAL

//! \brief macro exporting a symbol from the shared library
#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#if defined MFEMMGIS_EXPORTS
#define MFEM_MGIS_EXPORT MGIS_VISIBILITY_EXPORT
#else /* defined MFEMMGIS_EXPORTS */
#ifndef MFEM_MGIS_STATIC_BUILD
#define MFEM_MGIS_EXPORT MGIS_VISIBILITY_IMPORT
#else /* MFEM_MGIS_STATIC_BUILD */
#define MFEM_MGIS_EXPORT
#endif /* MFEM_MGIS_STATIC_BUILD */
#endif /* defined MFEMMGIS_EXPORTS */
#else  /* defined _WIN32 || defined _WIN64 || defined __CYGWIN__ */
#define MFEM_MGIS_EXPORT MGIS_VISIBILITY_EXPORT
#endif /* */

namespace mfem_mgis {

  using mgis::abort;
  using mgis::AbstractErrorHandler;
  using mgis::areInvalid;
  using mgis::areValid;
  using mgis::construct;
  using mgis::Context;
  using mgis::getDefaultVerbosityLevel;
  using mgis::InvalidResult;
  using mgis::invoke;
  using mgis::isInvalid;
  using mgis::isValid;
  using mgis::make_shared;
  using mgis::make_shared_as;
  using mgis::make_unique;
  using mgis::make_unique_as;
  using mgis::OptionalReference;
  using mgis::registerExceptionInErrorBacktrace;
  using mgis::terminate;
  using mgis::verboseDebug;
  using mgis::verboseFull;
  using mgis::verboseLevel0;
  using mgis::verboseLevel1;
  using mgis::verboseLevel2;
  using mgis::verboseLevel3;
  using mgis::VerbosityLevel;

  using mgis::debug;
  using mgis::getDefaultLogStream;
  using mgis::setDefaultLogStream;
  using mgis::warning;

  namespace attributes {
    //! \brief a simple alias
    using Throwing = ::mgis::attributes::ThrowingAttribute<true>;
    //! \brief a simple alias
    using MayThrow = ::mgis::attributes::ThrowingAttribute<true>;
    //! \brief a simple alias
    using MayAbort = ::mgis::attributes::AbortingAttribute<true>;
    //! \brief a simple alias
    using Unsafe = ::mgis::attributes::UnsafeAttribute;
  }  // namespace attributes
  //! \brief tag marking a call as unsafe without precautions
  inline constexpr auto unsafe = ::mgis::attributes::UnsafeAttribute{};
  //! \brief deprecated, use may_throw
  [[deprecated]] inline constexpr auto throwing =
      ::mgis::attributes::ThrowingAttribute<true>{};
  //! \brief tag marking a call that may throw
  inline constexpr auto may_throw =
      ::mgis::attributes::ThrowingAttribute<true>{};
  //! \brief tag marking a call that may abort
  inline constexpr auto may_abort =
      ::mgis::attributes::AbortingAttribute<true>{};

  //! \brief a simple alias
  using size_type = int;
  /*!
   * \brief constant stating that a number of components is not known at
   * compile-time
   */
  inline constexpr size_type dynamic_extent =
      std::numeric_limits<size_type>::max();
  //! \brief alias to the numeric type used
  using real = mgis::real;
  /*!
   * \brief this function can be called to report that the parallel
   * computations are not supported.
   */
  MFEM_MGIS_EXPORT [[noreturn]] void reportUnsupportedParallelComputations();
  //! \brief a simple alias
  using MainFunctionArguments = char**;
  /*!
   * \brief function that must be called to initialize `mfem-mgis`.
   * \param[in, out] argc: number of arguments
   * \param[in, out] argv: arguments
   *
   * In parallel, this function calls the `MPI_Init` function.
   * It is safe to call this function multiple times.
   *
   * The option `--execution <cpu|gpu|mixed>` (or `-x`) and, in the `mixed`
   * mode only, the options `--finite-elements` (`-fe`), `--behaviours`
   * (`-bh`), `--linear-solver-device` (`-lsd`) and `--post-processings`
   * (`-ppd`) are read first, before `MPI_Init` (see `ExecutionSettings`).
   * In the `gpu` and `mixed` modes, the device `gpu` of MFEM is created
   * unless `--device` gives another configuration using a GPU or the
   * device is already configured, which must then use a GPU; with
   * `--execution cpu`, a device using a GPU is an error. When `--execution`
   * is given, the environment variables `MFEM_MGIS_DEVICE` and
   * `MFEM_MGIS_HYPRE_ON_HOST` are ignored with a warning; in the `gpu` and
   * `mixed` modes, the memory of hypre follows the location of the linear
   * solvers (an option `--hypre-on-host` or `--no-hypre-on-host`
   * contradicting it is an error); `--hypre-device-pool-size` is an error
   * with `--execution cpu`, and PETSc is an error if the linear solvers are
   * placed on the GPU.
   *
   * If the option `--device <configuration>` (or `-dev`) is given, the MFEM
   * device is created with this configuration (for instance `cuda`), after
   * the initialization of hypre and before the one of PETSc. The processes
   * of a node use the GPUs of the node in turn. The memory of the host is
   * the standard one, and the memory of the device the default one of MFEM
   * (Umpire if MFEM is compiled with it), unless the environment variable
   * `MFEM_MEMORY` of MFEM gives them: MFEM compiled with Umpire otherwise
   * allocates the memory of the host with Umpire too, which slows down the
   * computations made on the host; an application which creates its own
   * `mfem::Device` may choose the memory types with
   * `mfem::Device::SetMemoryTypes` before. With a GPU, the profiling
   * sections wait for the end of the kernels launched on the GPU (see
   * `mgis::setProfilingSynchronization`). An application which uses
   * this option must neither create its own `mfem::Device` (MFEM refuses a
   * second configuration of the device) nor declare its own `--device`,
   * `-dev` options. Without this option, the configuration may be given by
   * the environment variable `MFEM_MGIS_DEVICE`, ignored if the device is
   * already configured or if the option `--execution` is given (to run a
   * whole test suite on a GPU, for instance).
   * An invalid configuration makes MFEM abort (`Invalid backend name`).
   * The variable `MFEM_DEVICE` of MFEM configures the device as soon as
   * MFEM is loaded: it takes precedence over `MFEM_MGIS_DEVICE` and makes
   * the option `--device` fail, and the device is then configured before
   * hypre, without the choice of the GPUs by process nor the
   * synchronisation of the profiling sections.
   *
   * The option `--hypre-device-pool-size <size>` (or `-hdps`), or else the
   * environment variable `MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE`, gives the size
   * in GiB of the pool of memory allocated at once by hypre on the GPU
   * through Umpire (4 GiB by default in hypre 3.2), which grows if needed:
   * a smaller pool leaves more memory to MFEM and to mfem-mgis. Each process
   * has its own pool (n times the size for n processes sharing a GPU). The
   * option is an error if hypre does not use Umpire on the GPU; the
   * environment variable is then ignored with a warning.
   *
   * The option `--hypre-on-host` (or `-hoh`, `--no-hypre-on-host` and
   * `-no-hoh` cancelling it), or else the environment variable
   * `MFEM_MGIS_HYPRE_ON_HOST` (any non empty value other than `0`), makes
   * hypre allocate its memory and compute on the host, even if the device
   * of MFEM uses a GPU: the behaviours and the computations of MFEM and of
   * mfem-mgis stay on the device, and MFEM copies the matrices and the
   * vectors given to hypre on the host. This is faster for small linear
   * systems, whose resolutions on a GPU are dominated by the launches of
   * kernels. The request is read before the initialization of hypre, so
   * that neither this initialization nor the configuration of the device
   * (by this function, by the program or by the variable `MFEM_DEVICE`)
   * make hypre use the GPU. The size of the pool of memory of hypre on the
   * GPU can't be given by the option `--hypre-device-pool-size` (the
   * environment variable is accepted). The option is an error if hypre is
   * built with GPU support and older than 2.31.0; the environment
   * variable is then ignored with a warning. PETSc may make hypre use the
   * GPU again (conversion of its matrices located on a GPU for the
   * preconditioners of hypre).
   *
   * The option `--kokkos` (or `-kk`, `--no-kokkos` and `-no-kk` cancelling
   * it), or else the environment variable `MFEM_MGIS_KOKKOS` (any non empty
   * value other than `0`), initializes Kokkos after the device of MFEM, on
   * its GPU, for the behaviours generated by the kokkos backend of the
   * generic-parallel interface. The GPU is the one of the device of MFEM
   * (configured by `--device`, by `--execution gpu` or `mixed`, by
   * `MFEM_MGIS_DEVICE`, by `MFEM_DEVICE` or by the program before the
   * call); otherwise, Kokkos chooses it. It requires mfem-mgis built with
   * Kokkos (CMake option `enable-kokkos`): otherwise, the option is an
   * error and the environment variable is ignored with a warning.
   */
  MFEM_MGIS_EXPORT void initialize(int& argc, MainFunctionArguments& argv);
  /*!
   * \brief function that must be called to end `mfem-mgis`. It finalizes
   * Kokkos if `initialize` initialized it, then MPI.
   * This call is optional if the code exits normally, but recommended if
   * `initialize` initialized Kokkos, which is otherwise finalized during the
   * destruction of the static objects.
   */
  MFEM_MGIS_EXPORT void finalize();
  /*!
   * \return if the current process is the main MPI process
   *
   * \note in sequential, this function always returns true
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool
  isMPIMainProcessForDefaultCommunicator() noexcept;
  /*!
   * \return the MPI rank for the default communicator
   *
   * \note in sequential, this function always returns 0
   */
  MFEM_MGIS_EXPORT [[nodiscard]] int
  getMPIRankForDefaultCommunicator() noexcept;
  /*!
   * \return the total number of MPI processes for the default communicator.
   *
   * \note in sequential, this function always returns 1
   */
  MFEM_MGIS_EXPORT [[nodiscard]] int
  getMPISizeForDefaultCommunicator() noexcept;
#ifdef MFEM_USE_MPI
  /*!
   * \brief gather the error messages of all processes in the context
   *
   * After this call, the context holds the same error message on all
   * processes:
   *
   * - if no process has an error message, the context is left unchanged.
   * - if all the processes having an error message have the same one, this
   *   message is kept unchanged.
   * - otherwise, the error messages are merged and each one is preceded by
   *   the rank of the process that registered it. Processes without error
   *   message are not listed.
   *
   * \param[in, out] ctx: execution context
   * \param[in] c: MPI communicator
   *
   * \note this function is collective: it must be called by all the
   * processes of the communicator.
   */
  MFEM_MGIS_EXPORT void synchronize(Context& ctx, const MPI_Comm c) noexcept;
#endif /* MFEM_USE_MPI */
  /*!
   * \brief a small wrapper used to build the exception outside the
   * `throw` statement. As most exception classes' constructors may
   * throw, this avoids undefined behaviour as reported by the
   * `cert-err60-cpp` warning of `clang-tidy` (thrown exception type
   * is not nothrow copy constructible).
   * \tparam Exception: type of the exception to be thrown.
   */
  template <typename Exception = std::runtime_error>
  [[noreturn]] MFEM_MGIS_VISIBILITY_LOCAL void raise();

  /*!
   * \brief a small wrapper used to build the exception outside the
   * `throw` statement. As most exception classes' constructors may
   * throw, this avoids undefined behaviour as reported by the
   * `cert-err60-cpp` warning of `clang-tidy` (thrown exception type
   * is not nothrow copy constructible).
   * \tparam Exception: type of the exception to be thrown.
   * \tparam Args: type of the arguments passed to the exception's
   * constructor.
   * \param[in] a: arguments passed to the exception's constructor.
   */
  template <typename Exception = std::runtime_error, typename... Args>
  [[noreturn]] MFEM_MGIS_VISIBILITY_LOCAL void raise(Args&&... a);
  /*!
   * \brief raise an exception if the first argument is `true`.
   * \tparam Exception: type of the exception to be thrown.
   * \tparam Args: type of the arguments passed to the exception's
   * constructor.
   * \param[in] b: condition to be checked. If `true`, an exception is
   * thrown.
   * \param[in] a: arguments passed to the exception's constructor.
   */
  template <typename Exception = std::runtime_error, typename... Args>
  MFEM_MGIS_VISIBILITY_LOCAL inline void raise_if(const bool b, Args&&... a);

  /*!
   * \brief function that must be called if one MPI process detects a fatal
   * error.
   * \param[in] error: exit status
   */
  MFEM_MGIS_EXPORT [[noreturn]] void abort(const int error = EXIT_FAILURE);
  /*!
   * \brief function that must be called if one MPI process detects a fatal
   * error.
   * \param[in] msg: message displayed by the calling process
   * \param[in] error: exit status
   */
  MFEM_MGIS_EXPORT [[noreturn]] void abort(const char* const msg,
                                           const int error = EXIT_FAILURE);
  //! \return if the usage of PETSc has been requested by the user.
  MFEM_MGIS_EXPORT bool usePETSc();
  /*!
   * \brief activate PETSc with the configuration file given in parameter.
   * \param[in] petscrc_file: PETSc configuration file
   */
  MFEM_MGIS_EXPORT void setPETSc(const char* petscrc_file);
  /*!
   * \brief declare the default options, read by `initialize`: the execution
   * settings (`--execution`, `--finite-elements`, `--behaviours`,
   * `--linear-solver-device`, `--post-processings`), the configuration of
   * the MFEM device (`--device`), the size of the pool of memory of hypre on
   * the GPU (`--hypre-device-pool-size`), the computations of hypre on the
   * host (`--hypre-on-host`), the initialization of Kokkos (`--kokkos`), and
   * the PETSc options (`--use-petsc`, `--petsc-configuration-file`), if
   * PETSc is available.
   * \param[in, out] parser: options parser
   */
  MFEM_MGIS_EXPORT void declareDefaultOptions(mfem::OptionsParser& parser);
  /*!
   * \return true if hypre computes on a GPU (hypre built with the support of
   * GPUs, its memory being located on the device), false otherwise and
   * without MPI. Contrary to `mfem::HypreUsingGPU`, which is only declared
   * when MFEM is compiled with MPI, it can be called without guards.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isHypreUsingGPU() noexcept;
  /*!
   * \return true if Kokkos is initialized (for instance by
   * `mfem_mgis::initialize` with the option `--kokkos`, see the CMake option
   * `enable-kokkos`), false otherwise and if mfem-mgis is built without
   * Kokkos
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isKokkosInitialized() noexcept;

  //! \return the output stream
  MFEM_MGIS_EXPORT std::ostream& getOutputStream();
  //! \return the error stream
  MFEM_MGIS_EXPORT std::ostream& getErrorStream();

}  // namespace mfem_mgis

#include "MFEMMGIS/Config.ixx"

#endif /* LIB_MFEM_MGIS_CONFIG_HXX */
