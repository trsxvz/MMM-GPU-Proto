/*!
 * \file   stand-ins/StandIns.hxx
 * \brief  Stand-ins for the types of CUDA, SYCL and Kokkos held by the
 * policies, so that the checks run on machines without these libraries:
 *
 * - `cudaStream_t` is a pointer to an opaque structure, as in
 *   `cuda_runtime_api.h`;
 * - `sycl::queue` holds a shared pointer to its implementation, as the
 *   queues of `icpx` and `AdaptiveCpp`;
 * - the execution spaces of Kokkos are declared without attribute of
 *   visibility, as in the headers of Kokkos, and give their name by a
 *   static function `name()`.
 */

#ifndef LIB_STANDINS_HXX
#define LIB_STANDINS_HXX

#include <memory>

// cuda_runtime_api.h
typedef struct CUstream_st* cudaStream_t;

namespace sycl {

  class queue {
   public:
    explicit queue(const int d) : impl(std::make_shared<int>(d)) {}
    int get_device() const { return *(this->impl); }
    long use_count() const { return this->impl.use_count(); }

   private:
    std::shared_ptr<int> impl;
  };

}  // end of namespace sycl

/*
 * name of the implementation of SYCL, read by the real header from
 * `SYCL_IMPLEMENTATION_ONEAPI` (icpx) or `SYCL_IMPLEMENTATION_ACPP`
 * (AdaptiveCpp)
 */
#define MFRONT_GENERICPARALLEL_SYCL_IMPLEMENTATION "stand-in"

namespace Kokkos {

  class Cuda {
   public:
    static const char* name() { return "Cuda"; }
    int instance = 0;
  };

  class Serial {
   public:
    static const char* name() { return "Serial"; }
    int instance = 0;
  };

}  // end of namespace Kokkos

#endif /* LIB_STANDINS_HXX */
