/*!
 * \file   mfront/include/MFront/GenericParallel/Behaviour/OpenACCBackend.hxx
 * \brief  This file declares the backend of the `generic-parallel` behaviour
 * interface based on OpenACC
 * \date   20/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_OPENACCBACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_OPENACCBACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/Behaviour/AbstractBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  /*!
   * \brief backend based on OpenACC.
   *
   * The points are processed by a `parallel loop` construct, on the device
   * selected by the compiler and the OpenACC runtime (a GPU with
   * `nvc++ -acc=gpu`, the cores of the host with `nvc++ -acc=multicore`).
   * The memory space of the library is deduced from the compilation options
   * (see `mfront::generic_parallel::getOpenACCMemorySpace`).
   *
   * Options:
   *
   * - `device_stack_size`: minimal stack size, in bytes, of the threads of the
   *   device when the computations are performed by an NVIDIA GPU (ignored
   *   otherwise).
   */
  struct MFRONT_VISIBILITY_EXPORT OpenACCBackend final : AbstractBackend {
    //! \brief default constructor
    OpenACCBackend();
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    OpenACCBackend(const tfel::utilities::DataMap&);
    //
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::string getSourceFileExtension() const override;
    void writeSpecificIncludes(std::ostream&) const override;
    [[nodiscard]] std::string getSharedDataCopy(
        std::string_view, std::string_view) const override;
    [[nodiscard]] std::string getApply(std::string_view) const override;
    [[nodiscard]] std::string getMemorySpace() const override;
    //! \brief destructor
    ~OpenACCBackend() noexcept override;

   private:
    //! \brief stack size of the threads of the device (0 if not specified)
    int device_stack_size = 0;
  };  // end of OpenACCBackend

}  // end of namespace mfront::generic_parallel::behaviour

#endif /* LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_OPENACCBACKEND_HXX */
