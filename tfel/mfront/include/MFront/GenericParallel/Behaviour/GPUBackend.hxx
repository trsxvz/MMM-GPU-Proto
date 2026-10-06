/*!
 * \file   mfront/include/MFront/GenericParallel/Behaviour/GPUBackend.hxx
 * \brief  This file declares the CUDA and HIP backends of the
 * `generic-parallel` behaviour interface
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_GPUBACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_GPUBACKEND_HXX

#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/Behaviour/AbstractBackend.hxx"

namespace mfront::generic_parallel::behaviour {

  /*!
   * \brief base class of the CUDA and HIP backends: one thread per
   * integration point.
   *
   * Options:
   *
   * - `number_of_threads_per_block`: number of threads per block, rounded up
   *   to a multiple of the warp size (128 by default).
   */
  struct MFRONT_VISIBILITY_EXPORT GPUBackend : AbstractBackend {
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    GPUBackend(const tfel::utilities::DataMap&);
    //
    void writeSpecificIncludes(std::ostream&) const override;
    [[nodiscard]] std::string getApply(std::string_view) const override;
    [[nodiscard]] std::string getSharedDataCopy(
        std::string_view, std::string_view) const override;
    //! \brief destructor
    ~GPUBackend() noexcept override;

   protected:
    //! \brief number of threads per block
    int number_of_threads_per_block = 128;
  };  // end of GPUBackend

  //! \brief CUDA backend
  struct MFRONT_VISIBILITY_EXPORT CUDABackend final : GPUBackend {
    using GPUBackend::GPUBackend;
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::string getSourceFileExtension() const override;
    //! \brief destructor
    ~CUDABackend() noexcept override;
  };  // end of CUDABackend

  //! \brief HIP backend
  struct MFRONT_VISIBILITY_EXPORT HIPBackend final : GPUBackend {
    using GPUBackend::GPUBackend;
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::string getSourceFileExtension() const override;
    void writeSpecificIncludes(std::ostream&) const override;
    //! \brief destructor
    ~HIPBackend() noexcept override;
  };  // end of HIPBackend

}  // end of namespace mfront::generic_parallel::behaviour

#endif /* LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_GPUBACKEND_HXX */
