/*!
 * \file   src/Batch/HostThreads.hxx
 * \brief  This file implements the loops of the host shared between the
 * threads of OpenMP (internal header, host code)
 * \date   16/11/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_HOSTTHREADS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_HOSTTHREADS_HXX

#ifdef MFEM_MGIS_HAVE_OPENMP
#include <omp.h>
#include <cstdint>
#include <exception>
#endif /* MFEM_MGIS_HAVE_OPENMP */

namespace mfem_mgis::batch::internals {

  /*!
   * \brief calls `f(begin, end)` on consecutive ranges covering `[0, n)`.
   *
   * When mfem-mgis is built with OpenMP (option `enable-openmp`), each
   * thread treats one range, the ranges having the same size up to one
   * index. The threads are the ones of OpenMP: their number and their
   * placement are given by its settings (`OMP_NUM_THREADS`, by default the
   * processors given to the process by its launcher, `OMP_PROC_BIND`,
   * `OMP_PLACES`). Otherwise, `f` is called once, with the whole range.
   *
   * The results must not depend on the ranges: a call only writes values
   * which belong to the indices of its range, and objects of its own (work
   * arrays, declared in `f`). `f` is not called on an empty range. If calls
   * throw, one of the exceptions is thrown again when all of them returned.
   * \param[in] n: number of indices
   * \param[in] f: function treating a range of indices
   */
  template <typename Index, typename Function>
  void forEachThreadRange(const Index n, const Function& f) {
#ifdef MFEM_MGIS_HAVE_OPENMP
    auto error = std::exception_ptr{};
#pragma omp parallel
    {
      const auto threads = static_cast<std::uint64_t>(omp_get_num_threads());
      const auto thread = static_cast<std::uint64_t>(omp_get_thread_num());
      const auto size = static_cast<std::uint64_t>(n);
      const auto begin = static_cast<Index>(size * thread / threads);
      const auto end = static_cast<Index>(size * (thread + 1) / threads);
      if (begin != end) {
        try {
          f(begin, end);
        } catch (...) {
#pragma omp critical(mfem_mgis_for_each_thread_range)
          if (error == nullptr) {
            error = std::current_exception();
          }
        }
      }
    }
    if (error != nullptr) {
      std::rethrow_exception(error);
    }
#else  /* MFEM_MGIS_HAVE_OPENMP */
    if (n != Index{}) {
      f(Index{}, n);
    }
#endif /* MFEM_MGIS_HAVE_OPENMP */
  }  // end of forEachThreadRange

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_HOSTTHREADS_HXX */
