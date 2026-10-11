/*!
 * \file   src/PostProcessingsComputation.cxx
 * \brief  This file implements the choice of the computation of the
 * post-processings
 * \date   11/10/2026
 */

#include <atomic>
#include <string>
#include <cstdlib>
#include <string_view>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/PostProcessingsComputation.hxx"
#include "PostProcessingsComputation.hxx"

namespace mfem_mgis::internals {

  //! \brief number of the computations made by kernels
  static std::atomic<size_type>
      number_of_post_processings_computations_by_kernels{0};

  void countPostProcessingComputationByKernels() noexcept {
    ++number_of_post_processings_computations_by_kernels;
  }  // end of countPostProcessingComputationByKernels

  bool usePostProcessingsKernels(attributes::MayThrow) {
    const auto* const v = "MFEM_MGIS_POST_PROCESSINGS_COMPUTATION";
    const auto* const e = std::getenv(v);
    if (e == nullptr) {
      return false;
    }
    const auto value = std::string_view{e};
    if (value == "ElementByElement") {
      return false;
    }
    if (value == "Kernels") {
      return true;
    }
    mgis::raise("invalid value '" + std::string{value} +
                "' for environment variable '" + std::string{v} +
                "' (expected 'ElementByElement' or 'Kernels')");
  }  // end of usePostProcessingsKernels

}  // end of namespace mfem_mgis::internals

namespace mfem_mgis {

  size_type getNumberOfPostProcessingsComputationsByKernels() noexcept {
    return internals::number_of_post_processings_computations_by_kernels;
  }  // end of getNumberOfPostProcessingsComputationsByKernels

}  // end of namespace mfem_mgis
