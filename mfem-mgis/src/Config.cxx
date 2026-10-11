/*!
 * \file   src/Config.cxx
 * \brief  This file implements the functions declared in `MFEMMGIS/Config.hxx`
 * \author Thomas Helfer
 * \date   14/02/2021
 */

#include <cmath>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <string_view>
#include <iostream>
#ifdef MFEM_USE_MPI
#include "mpi.h"
#include "mfem/linalg/hypre.hpp"
#include "HYPRE_utilities.h"
#endif /* MFEM_USE_MPI */
#include "mfem/general/error.hpp"
#include "mfem/general/device.hpp"
#include "mfem/general/backends.hpp"
#include "mfem/general/optparser.hpp"
#include "mfem/general/communication.hpp"
#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /*MFEM_USE_PETSC */
#include "MGIS/Raise.hxx"
#include "MGIS/Profiling.hxx"
#include <MFEMMGIS/Profiler.hxx>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief structure in charge of freeing resources on exit.
   */
  struct MGIS_VISIBILITY_LOCAL Finalizer {
    //! \brief option used to select the PETSc configuration file
    static const char* const petsc_configuration_file_option;
    //! \brief option used to select the configuration of the MFEM device
    static const char* const device_option;
    //! \brief short name of the option `device_option`
    static const char* const device_short_option;
    /*!
     * \brief option giving the size, in GiB, of the pool of memory of
     * hypre on the GPU
     */
    static const char* const hypre_device_pool_size_option;
    //! \brief short name of the option `hypre_device_pool_size_option`
    static const char* const hypre_device_pool_size_short_option;
    //! \return the unique instance of this class
    static Finalizer& get();
    /*!
     * \brief initialize the execution of the mfem-mgis
     * \param[in] argc: number of arguments
     * \param[in] argv: arguments
     */
    void initialize(int&, MainFunctionArguments&);
    //! \return true if PETSc is used
    bool usePETSc() const;
    /*!
     * \brief activate PETSc and provide a configuration file
     * \param[in] petscrc_file: PETSc configuration file
     */
    void setPETSc(const char* petscrc_file);
    //! \brief finalize the execution of the mfem-mgis
    void finalize();
    /*!
     * \brief abort the process
     * \param[in] error: exit status
     */
    [[noreturn]] void abort(int error);

   private:
    //! \brief boolean stating if PETSc shall be used
    bool use_petsc = false;
    /*!
     * \brief MFEM device, created if the `--device` option is given. It is
     * kept until the end of the program, since MFEM objects using it may
     * be destroyed after `finalize`.
     */
    std::unique_ptr<mfem::Device> device;
    //! \brief boolean stating if MPI has been finalized or aborted
    bool pendingExit = false;
    /*!
     * \brief constructor
     */
    Finalizer();
    //! \brief destructor
    ~Finalizer();
  };  // end of Finalizer

  const char* const Finalizer::petsc_configuration_file_option =
      "--petsc-configuration-file";

  const char* const Finalizer::device_option = "--device";

  const char* const Finalizer::device_short_option = "-dev";

  const char* const Finalizer::hypre_device_pool_size_option =
      "--hypre-device-pool-size";

  const char* const Finalizer::hypre_device_pool_size_short_option = "-hdps";

  /*!
   * \brief sets the size of the pool of memory allocated at once by hypre
   * on the GPU through Umpire (4 GiB by default in hypre 3.2), before the
   * first allocation of hypre on the GPU, which creates the pool. Each
   * process has its own pool.
   * \param[in] size: size in GiB, as given by the user
   * \param[in] environment: if true, the size is given by the environment
   * variable `MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE`, ignored with a warning if
   * hypre does not use Umpire on the GPU (as for a whole test suite)
   */
  static void setHypreDevicePoolSize(const char* const size,
                                     const bool environment) {
    auto* end = static_cast<char*>(nullptr);
    const auto s = std::strtod(size, &end);
    // number of bytes, which must be representable
    const auto bytes = s * 1073741824.;
    if ((end == size) || (*end != '\0') || (!std::isfinite(s)) ||
        (!(bytes >= 1)) ||
        (!(bytes <
           static_cast<double>(std::numeric_limits<std::size_t>::max())))) {
      mgis::raise("initialize: invalid size of the pool of memory of hypre '" +
                  std::string{size} +
                  "' (strictly positive number of GiB expected)");
    }
#if defined(MFEM_USE_MPI) && defined(HYPRE_USING_UMPIRE_DEVICE)
    static_cast<void>(environment);
    // the function only stores the size: its result, the global error flag
    // of hypre, may report a previous error
    static_cast<void>(
        HYPRE_SetUmpireDevicePoolSize(static_cast<std::size_t>(bytes)));
#else
    if (!environment) {
      mgis::raise(
          "initialize: hypre does not use a pool of memory of Umpire on the "
          "GPU");
    }
    mfem_mgis::getErrorStream()
        << "initialize: hypre does not use a pool of memory of Umpire on the "
           "GPU, the variable MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE is ignored\n";
#endif
  }  // end of setHypreDevicePoolSize

  /*!
   * \return the identifier of the GPU used by this process: the processes
   * of a node use the GPUs of the node in turn (several processes may share
   * a GPU)
   * \param[in] configuration: configuration of the device
   */
  static int getDeviceIdentifier(const char* const configuration) {
#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
    // MFEM aborts if the number of GPUs is queried on a machine without GPU:
    // it is only queried for a configuration using a GPU (`cuda`, `hip`,
    // `gpu` and their variants `ceed-cuda`, `raja-gpu`, ...)
    const auto c = std::string_view{configuration};
    if ((c.find("cuda") == std::string_view::npos) &&
        (c.find("hip") == std::string_view::npos) &&
        (c.find("gpu") == std::string_view::npos)) {
      return 0;
    }
    const auto n = mfem::Device::GetDeviceCount();
    if (n <= 0) {
      return 0;
    }
#ifdef MFEM_USE_MPI
    auto local = MPI_Comm{};
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                        &local);
    auto r = int{};
    MPI_Comm_rank(local, &r);
    MPI_Comm_free(&local);
    return r % n;
#else  /* MFEM_USE_MPI */
    return 0;
#endif /* MFEM_USE_MPI */
#else  /* defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP) */
    // without GPU support, MFEM aborts if the number of GPUs is queried
    static_cast<void>(configuration);
    return 0;
#endif /* defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP) */
  }  // end of getDeviceIdentifier

  /*!
   * \brief waits for the end of the kernels launched on the GPU, so that
   * their duration is counted in the profiling section which launched them
   */
  static void synchronizeDevice() noexcept {
#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
    try {
      MFEM_DEVICE_SYNC;
    } catch (...) {
      // an error of the GPU is reported by the next call to the GPU
    }
#endif /* defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP) */
  }  // end of synchronizeDevice

  Finalizer::Finalizer() = default;

#ifdef MFEM_USE_PETSC
  void Finalizer::setPETSc(const char* petscrc_file) {
    this->use_petsc = true;
    mfem::MFEMInitializePetsc(nullptr, nullptr, petscrc_file, nullptr);
  }
#else  /* MFEM_USE_PETSC */
  void Finalizer::setPETSc(const char*) {
    mgis::raise("Finalizer::setPETSc: PETSc is not supported");
  }
#endif /*MFEM_USE_PETSC*/

  void Finalizer::initialize(int& argc, MainFunctionArguments& argv) {
#ifdef MFEM_USE_PETSC
    const char* petscrc_file = nullptr;
#endif /* MFEM_USE_PETSC */
    const char* device_configuration = nullptr;
    const char* hypre_device_pool_size = nullptr;
    auto hypre_device_pool_size_from_environment = false;
    for (const auto* a = argv; a != argv + argc; ++a) {
      if ((std::strcmp(*a, hypre_device_pool_size_option) == 0) ||
          (std::strcmp(*a, hypre_device_pool_size_short_option) == 0)) {
        if (hypre_device_pool_size != nullptr) {
          mgis::raise(
              "initialize: size of the pool of memory of hypre already "
              "specified");
        }
        ++a;
        if (a == argv + argc) {
          mgis::raise("initialize: no value given to " +
                      std::string{hypre_device_pool_size_option});
        }
        hypre_device_pool_size = *a;
        continue;
      }
      if ((std::strcmp(*a, device_option) == 0) ||
          (std::strcmp(*a, device_short_option) == 0)) {
        if (device_configuration != nullptr) {
          mgis::raise("initialize: MFEM device already specified");
        }
        ++a;
        if (a == argv + argc) {
          mgis::raise("initialize: no value given to --device");
        }
        device_configuration = *a;
        continue;
      }
#ifdef MFEM_USE_PETSC
      // long and short names, as declared by `declareDefaultOptions`
      if ((std::strcmp(*a, "--use-petsc") == 0) ||
          (std::strcmp(*a, "-up") == 0)) {
        this->use_petsc = true;
      }
      if ((std::strcmp(*a, "--donot-use-petsc") == 0) ||
          (std::strcmp(*a, "-no-up") == 0)) {
        this->use_petsc = false;
      }
      if ((std::strcmp(*a, petsc_configuration_file_option) == 0) ||
          (std::strcmp(*a, "petscrc_file") == 0)) {
        if (petscrc_file != nullptr) {
          mgis::raise("initialize: PETSc configuration file already specified");
        }
        ++a;
        if (a == argv + argc) {
          mgis::raise(
              "initialize: no value given to --petsc-configuration-file");
        }
        petscrc_file = *a;
      }
#endif /* MFEM_USE_PETSC */
    }
    // size of the pool of memory of hypre on the GPU, given by the option or
    // by the environment
    if (hypre_device_pool_size == nullptr) {
      const auto* const e = std::getenv("MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE");
      if ((e != nullptr) && (*e != '\0')) {
        hypre_device_pool_size = e;
        hypre_device_pool_size_from_environment = true;
      }
    }
    if (hypre_device_pool_size != nullptr) {
      setHypreDevicePoolSize(hypre_device_pool_size,
                             hypre_device_pool_size_from_environment);
    }
    // without the option, the configuration of the device may be given by
    // the environment, for instance to run a whole test suite on a GPU
    if (device_configuration == nullptr) {
      const auto* const e = std::getenv("MFEM_MGIS_DEVICE");
      if ((e != nullptr) && (*e != '\0') && (!mfem::Device::IsConfigured())) {
        device_configuration = e;
      }
    }
    // the device is configured after the initialization of hypre and
    // before the one of PETSc, which uses the GPU of the device, as in the
    // parallel examples of MFEM
    if (device_configuration != nullptr) {
      if (mfem::Device::IsConfigured()) {
        mgis::raise(
            "initialize: the MFEM device is already configured, the "
            "option '--device' can't be used");
      }
      this->device = std::make_unique<mfem::Device>(
          device_configuration, getDeviceIdentifier(device_configuration));
      // the kernels are launched asynchronously on a GPU: the profiling
      // sections wait for them
      if (mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                               mfem::Backend::HIP_MASK)) {
        mgis::setProfilingSynchronization(&synchronizeDevice);
      }
    }
#ifdef MFEM_USE_PETSC
    if (this->use_petsc) {
      if (petscrc_file == nullptr) {
        mgis::raise("initialize: no PETSc configuration file given");
      }
      mfem::MFEMInitializePetsc(nullptr, nullptr, petscrc_file, nullptr);
    }
#endif /* MFEM_USE_PETSC */
  }    // end of initialize

  Finalizer& Finalizer::get() {
    static Finalizer f;
    return f;
  }  // end of get

  bool Finalizer::usePETSc() const {
    return this->use_petsc;
  }  // end of usePETSc

  void Finalizer::finalize() {
    if (!this->pendingExit) {
#ifdef MFEM_USE_MPI
#ifdef MFEM_USE_PETSC
      if (this->use_petsc) {
        mfem::MFEMFinalizePetsc();
      }
#endif /* MFEM_USE_PETSC */
      mfem::Hypre::Finalize();
      MPI_Finalize();
      this->pendingExit = true;
#endif /* MFEM_USE_MPI */
    }
  }  // end of finalize

  [[noreturn]] void Finalizer::abort(int error) {
    // MPI_Abort does not flush the output streams
    mfem_mgis::getOutputStream().flush();
    mfem_mgis::getErrorStream().flush();
    if (!this->pendingExit) {
#ifdef MFEM_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, error);
      this->pendingExit = true;
#endif /* MFEM_USE_MPI */
    }
    std::exit(error);
  }  // end of abort

  Finalizer::~Finalizer() { this->finalize(); }

  [[noreturn]] void reportUnsupportedParallelComputations() {
    mfem::mfem_error(
        "reportUnsupportedParallelComputations: "
        "unsupported parallel computations");
  }  // end of reportUnsupportedParallelComputations

#ifdef MFEM_USE_MPI

  //   [[noreturn]] static void exit_on_failure() {
  //     try {
  //       throw;
  //     } catch (std::exception& e) {
  //       mfem_mgis::getErrorStream() << e.what() << '\n';
  //     } catch (...) {
  //       mfem_mgis::getErrorStream() << "unknown exception thrown";
  //     }
  //     abort();
  //   }  // end of exit_on_failure

  void initialize(int& argc, MainFunctionArguments& argv) {
    static bool first = true;
    if (first) {
      mgis::setDefaultLogStream(mfem_mgis::getOutputStream());
      MPI_Init(&argc, &argv);
      mfem::Hypre::Init();
      if (getMPIrank() != 0) {
        mfem::out.Disable();
        mfem::err.Disable();
      }
      Finalizer::get().initialize(argc, argv);
      first = false;
    }
  }  // end of initialize

#else /* MFEM_USE_MPI */

  void initialize(int& argc, MainFunctionArguments& argv) {
    static bool first = true;
    if (first) {
      Finalizer::get().initialize(argc, argv);
      first = false;
    }
  }  // end of initialize

#endif /* MFEM_USE_MPI */

  void finalize() { Finalizer::get().finalize(); }  // end of finalize

  void abort(const int error) {
    Finalizer::get().abort(error);
  }  // end of abort

  void abort(const char* const msg, const int error) {
    mfem_mgis::getErrorStream() << msg << '\n';
    Finalizer::get().abort(error);
  }  // end of abort

  bool isHypreUsingGPU() noexcept {
#ifdef MFEM_USE_MPI
    return mfem::HypreUsingGPU();
#else  /* MFEM_USE_MPI */
    return false;
#endif /* MFEM_USE_MPI */
  }  // end of isHypreUsingGPU

  void declareDefaultOptions(mfem::OptionsParser& parser) {
    // the device is created by `initialize`: the option is only declared
    // here, so that the parser accepts it
    static const char* device_configuration = "";
    parser.AddOption(&device_configuration, Finalizer::device_short_option,
                     Finalizer::device_option,
                     "Configuration of the MFEM device (for instance cpu, "
                     "cuda, hip, omp or debug), created by "
                     "mfem_mgis::initialize. No device is created by "
                     "default.");
    // read by `initialize`: only declared here, as a string, so that the
    // parser accepts it
    static const char* hypre_device_pool_size = "";
    parser.AddOption(&hypre_device_pool_size,
                     Finalizer::hypre_device_pool_size_short_option,
                     Finalizer::hypre_device_pool_size_option,
                     "Size in GiB of the pool of memory allocated at once by "
                     "hypre on the GPU (Umpire), set by "
                     "mfem_mgis::initialize. Default of hypre if not given.");
#ifdef MFEM_USE_PETSC
    static bool use_petsc = false;
    static const char* petscrc_file = "";
#endif /* MFEM_USE_PETSC */
#ifdef MFEM_USE_PETSC
    parser.AddOption(&use_petsc, "-up", "--use-petsc", "-no-up",
                     "--donot-use-petsc",
                     "Use or not PETSc to solve the nonlinear system.");
    parser.AddOption(&petscrc_file, "petscrc_file",
                     Finalizer::petsc_configuration_file_option,
                     "Path to the PETSc configuration file.");
#endif /* MFEM_USE_PETSC */
  }

  bool usePETSc() { return Finalizer::get().usePETSc(); }  // end of usePETSc

  void setPETSc(const char* petscrc_file) {
    Finalizer::get().setPETSc(petscrc_file);
  }  // end of setPETSc

  int getMPIrank() {
    int rank = 0;
#ifdef MFEM_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif /* MFEM_USE_MPI */
    return rank;
  }  // end of getMPIrank

  int getMPIsize() {
    int size = 1;
#ifdef MFEM_USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif /* MFEM_USE_MPI */
    return size;
  }  // end of getMPIsize

}  // namespace mfem_mgis
