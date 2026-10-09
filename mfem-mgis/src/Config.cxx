/*!
 * \file   src/Config.cxx
 * \brief  This file implements the functions declared in `MFEMMGIS/Config.hxx`
 * \author Thomas Helfer
 * \date   14/02/2021
 */

#include <cstring>
#include <memory>
#include <iostream>
#ifdef MFEM_USE_MPI
#include "mpi.h"
#include "mfem/linalg/hypre.hpp"
#endif /* MFEM_USE_MPI */
#include "mfem/general/error.hpp"
#include "mfem/general/device.hpp"
#include "mfem/general/optparser.hpp"
#include "mfem/general/communication.hpp"
#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /*MFEM_USE_PETSC */
#include "MGIS/Raise.hxx"
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
    for (const auto* a = argv; a != argv + argc; ++a) {
      if (std::strcmp(*a, device_option) == 0) {
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
      if (std::strcmp(*a, "--use-petsc") == 0) {
        this->use_petsc = true;
      }
      if (std::strcmp(*a, petsc_configuration_file_option) == 0) {
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
#ifdef MFEM_USE_PETSC
    if (this->use_petsc) {
      if (petscrc_file == nullptr) {
        mgis::raise("initialize: no PETSc configuration file given");
      }
      mfem::MFEMInitializePetsc(nullptr, nullptr, petscrc_file, nullptr);
    }
#endif /* MFEM_USE_PETSC */
    // the device is configured after the initialization of hypre, as in
    // the parallel examples of MFEM
    if (device_configuration != nullptr) {
      this->device = std::make_unique<mfem::Device>(device_configuration);
    }
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

  void declareDefaultOptions(mfem::OptionsParser& parser) {
    // the device is created by `initialize`: the option is only declared
    // here, so that the parser accepts it
    static const char* device_configuration = "cpu";
    parser.AddOption(&device_configuration, "-dev", Finalizer::device_option,
                     "Configuration of the MFEM device (for instance cpu, "
                     "cuda, hip, omp or debug), created by "
                     "mfem_mgis::initialize.");
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
