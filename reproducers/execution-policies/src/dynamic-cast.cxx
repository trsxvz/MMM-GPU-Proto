/*!
 * \file   src/dynamic-cast.cxx
 * \brief  For information only: the same check by `dynamic_cast`, which
 * the design does not use. Loaded at run time, as a generated library, but
 * linked to the runtime library: `dynamic_cast` needs the description of
 * the type of the base class, defined there.
 */

#include "MFront/GenericParallel/BatchData.h"
#include "MFront/GenericParallel/CUDABackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/KokkosBackendExecutionPolicy.hxx"

namespace mgp = mfront::generic_parallel;

extern "C" {

MFRONT_GENERICPARALLEL_GENERATED_VISIBILITY int reproducer_dynamic_cast_cuda(
    const mfront_gp_AbstractBackendExecutionPolicy* const p) {
  return dynamic_cast<const mgp::CUDABackendExecutionPolicy*>(p) != nullptr;
}

MFRONT_GENERICPARALLEL_GENERATED_VISIBILITY int
reproducer_dynamic_cast_kokkos_cuda(
    const mfront_gp_AbstractBackendExecutionPolicy* const p) {
  using Policy = mgp::KokkosBackendExecutionPolicy<Kokkos::Cuda>;
  return dynamic_cast<const Policy*>(p) != nullptr;
}

}  // end of extern "C"
