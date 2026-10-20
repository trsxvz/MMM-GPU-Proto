/*!
 * \file   reproducers/mgis-gpu-headers/check.cxx
 * \brief  check of the headers of MGIS compiled by the calling code with a
 * GPU programming model: field operations (MGIS/FieldOperations) and memory
 * resources (MGIS/MemoryResources). The backend is chosen by one of the
 * macros MGIS_CHECK_CUDA, MGIS_CHECK_HIP, MGIS_CHECK_SYCL, MGIS_CHECK_KOKKOS,
 * MGIS_CHECK_OPENMP, MGIS_CHECK_STLPAR (parallel algorithms of the C++
 * standard library compiled by nvc++ -stdpar=gpu) and MGIS_CHECK_RAJA (loops
 * of RAJA compiled by nvcc for CUDA), the last two with the memory of CUDA
 * (see scripts/check-mgis-gpu-headers.sh). The memory resource is the one of
 * the backend (managed memory of CUDA with nvc++, of the resource of camp
 * with RAJA).
 *
 * The field operations copy and fill views located in the memory of the
 * device only: contiguous copy, transposition between the point by point and
 * the component by component layouts, padding, fill of one component. The
 * memory resource allocates memory read by the host and written by the
 * field operations, directly and through a pool of the standard library,
 * which requests blocks aligned on their size.
 */

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <memory_resource>
#include "MGIS/FieldView.hxx"

#if defined(MGIS_CHECK_CUDA)
#include "MGIS/FieldOperations/CUDAFieldOperations.hxx"
#include "MGIS/MemoryResources/CUDAManagedMemoryResource.hxx"
#elif defined(MGIS_CHECK_HIP)
#include "MGIS/FieldOperations/HIPFieldOperations.hxx"
#include "MGIS/MemoryResources/HIPManagedMemoryResource.hxx"
#elif defined(MGIS_CHECK_SYCL)
#include "MGIS/FieldOperations/SYCLFieldOperations.hxx"
#include "MGIS/MemoryResources/SYCLSharedMemoryResource.hxx"
#elif defined(MGIS_CHECK_KOKKOS)
#include "MGIS/FieldOperations/KokkosFieldOperations.hxx"
#include "MGIS/MemoryResources/KokkosSharedSpaceMemoryResource.hxx"
#elif defined(MGIS_CHECK_OPENMP)
#include "MGIS/FieldOperations/OpenMPFieldOperations.hxx"
#elif defined(MGIS_CHECK_RAJA)
#include <cuda_runtime.h>
#include "MGIS/FieldOperations/RAJAFieldOperations.hxx"
#include "MGIS/MemoryResources/RAJAManagedMemoryResource.hxx"
#elif defined(MGIS_CHECK_STLPAR)
#include <cuda_runtime.h>
#include "MGIS/FieldOperations/ParallelSTLFieldOperations.hxx"
#include "MGIS/MemoryResources/CUDAManagedMemoryResource.hxx"
#else
#error "no backend selected"
#endif

using mgis::real;
using mgis::size_type;
using mgis::FieldView;
using mgis::FieldLayout;

static bool success = true;

static void check(const bool b, const char* const msg) {
  if (!b) {
    success = false;
    std::cerr << msg << '\n';
  }
}

//! \brief memory of the device and synchronization, for each backend
struct Backend {
#if defined(MGIS_CHECK_SYCL)
  sycl::queue queue{sycl::property::queue::in_order{}};
#endif
  real* allocate(const size_type n) {
    void* p = nullptr;
#if defined(MGIS_CHECK_CUDA) || defined(MGIS_CHECK_STLPAR) || \
    defined(MGIS_CHECK_RAJA)
    if (cudaMalloc(&p, n * sizeof(real)) != cudaSuccess) {
      throw std::runtime_error("cudaMalloc failed");
    }
#elif defined(MGIS_CHECK_HIP)
    if (hipMalloc(&p, n * sizeof(real)) != hipSuccess) {
      throw std::runtime_error("hipMalloc failed");
    }
#elif defined(MGIS_CHECK_SYCL)
    p = sycl::malloc_device<real>(n, this->queue);
#elif defined(MGIS_CHECK_KOKKOS)
    p = Kokkos::kokkos_malloc<Kokkos::DefaultExecutionSpace::memory_space>(
        n * sizeof(real));
#elif defined(MGIS_CHECK_OPENMP)
    p = omp_target_alloc(n * sizeof(real), omp_get_default_device());
#endif
    if (p == nullptr) {
      throw std::runtime_error("allocation failed");
    }
    return static_cast<real*>(p);
  }
  void free(real* const p) {
#if defined(MGIS_CHECK_CUDA) || defined(MGIS_CHECK_STLPAR) || \
    defined(MGIS_CHECK_RAJA)
    static_cast<void>(cudaFree(p));
#elif defined(MGIS_CHECK_HIP)
    static_cast<void>(hipFree(p));
#elif defined(MGIS_CHECK_SYCL)
    sycl::free(p, this->queue);
#elif defined(MGIS_CHECK_KOKKOS)
    Kokkos::kokkos_free<Kokkos::DefaultExecutionSpace::memory_space>(p);
#elif defined(MGIS_CHECK_OPENMP)
    omp_target_free(p, omp_get_default_device());
#endif
  }
  void toDevice(real* const d, const std::vector<real>& h) {
    [[maybe_unused]] const auto s = h.size() * sizeof(real);
#if defined(MGIS_CHECK_CUDA) || defined(MGIS_CHECK_STLPAR) || \
    defined(MGIS_CHECK_RAJA)
    static_cast<void>(cudaMemcpy(d, h.data(), s, cudaMemcpyHostToDevice));
#elif defined(MGIS_CHECK_HIP)
    static_cast<void>(hipMemcpy(d, h.data(), s, hipMemcpyHostToDevice));
#elif defined(MGIS_CHECK_SYCL)
    this->queue.memcpy(d, h.data(), s).wait();
#elif defined(MGIS_CHECK_KOKKOS)
    Kokkos::deep_copy(
        Kokkos::View<real*, Kokkos::DefaultExecutionSpace::memory_space,
                     Kokkos::MemoryUnmanaged>(d, h.size()),
        Kokkos::View<const real*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
            h.data(), h.size()));
#elif defined(MGIS_CHECK_OPENMP)
    omp_target_memcpy(d, h.data(), s, 0, 0, omp_get_default_device(),
                      omp_get_initial_device());
#endif
  }
  std::vector<real> toHost(const real* const d, const size_type n) {
    auto h = std::vector<real>(n);
    [[maybe_unused]] const auto s = n * sizeof(real);
#if defined(MGIS_CHECK_CUDA) || defined(MGIS_CHECK_STLPAR) || \
    defined(MGIS_CHECK_RAJA)
    static_cast<void>(cudaMemcpy(h.data(), d, s, cudaMemcpyDeviceToHost));
#elif defined(MGIS_CHECK_HIP)
    static_cast<void>(hipMemcpy(h.data(), d, s, hipMemcpyDeviceToHost));
#elif defined(MGIS_CHECK_SYCL)
    this->queue.memcpy(h.data(), d, s).wait();
#elif defined(MGIS_CHECK_KOKKOS)
    Kokkos::deep_copy(
        Kokkos::View<real*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
            h.data(), n),
        Kokkos::View<const real*, Kokkos::DefaultExecutionSpace::memory_space,
                     Kokkos::MemoryUnmanaged>(d, n));
#elif defined(MGIS_CHECK_OPENMP)
    omp_target_memcpy(h.data(), d, s, 0, 0, omp_get_initial_device(),
                      omp_get_default_device());
#endif
    return h;
  }
  void synchronize() {
#if defined(MGIS_CHECK_CUDA) || defined(MGIS_CHECK_STLPAR) || \
    defined(MGIS_CHECK_RAJA)
    static_cast<void>(cudaDeviceSynchronize());
#elif defined(MGIS_CHECK_HIP)
    static_cast<void>(hipDeviceSynchronize());
#elif defined(MGIS_CHECK_SYCL)
    this->queue.wait_and_throw();
#elif defined(MGIS_CHECK_KOKKOS)
    Kokkos::fence();
#endif
  }
#if defined(MGIS_CHECK_CUDA)
  mgis::field_operations::CUDAFieldOperations operations() {
    return mgis::field_operations::CUDAFieldOperations{};
  }
#elif defined(MGIS_CHECK_HIP)
  mgis::field_operations::HIPFieldOperations operations() {
    return mgis::field_operations::HIPFieldOperations{};
  }
#elif defined(MGIS_CHECK_SYCL)
  mgis::field_operations::SYCLFieldOperations operations() {
    return mgis::field_operations::SYCLFieldOperations{this->queue};
  }
#elif defined(MGIS_CHECK_KOKKOS)
  mgis::field_operations::KokkosFieldOperations<> operations() {
    return mgis::field_operations::KokkosFieldOperations<>{};
  }
#elif defined(MGIS_CHECK_OPENMP)
  mgis::field_operations::OpenMPFieldOperations operations() {
    return mgis::field_operations::OpenMPFieldOperations{};
  }
#elif defined(MGIS_CHECK_STLPAR)
  mgis::field_operations::ParallelSTLFieldOperations<> operations() {
    return mgis::field_operations::ParallelSTLFieldOperations<>{};
  }
#elif defined(MGIS_CHECK_RAJA)
  // default execution policy: the one of CUDA (code compiled by nvcc)
  mgis::field_operations::RAJAFieldOperations<> operations() {
    static_assert(
        RAJA::type_traits::is_cuda_policy<
            mgis::field_operations::RAJADefaultExecutionPolicy>::value,
        "the default execution policy of RAJA is not the one of CUDA");
    return mgis::field_operations::RAJAFieldOperations<>{};
  }
#endif
};

static void test(Backend& backend) {
  constexpr auto n = size_type{37};
  constexpr auto nc = size_type{3};
  auto values = std::vector<real>(n * nc);
  for (size_type k = 0; k != values.size(); ++k) {
    values[k] = static_cast<real>(k);
  }
  auto* const a = backend.allocate(n * nc);         // point by point
  auto* const c = backend.allocate(nc * (n + 5));   // component by component
  auto* const d = backend.allocate(n * nc);         // point by point
  backend.toDevice(a, values);
  auto ops = backend.operations();
  const auto va = mgis::makeFieldView(std::span<real>(a, n * nc),
                                      FieldLayout::POINT_BY_POINT, n, nc);
  const auto sc = std::span<real>(c, nc * (n + 5));
  const auto vc = mgis::makeComponentByComponentFieldView(sc, n, nc, n + 5);
  const auto vd = mgis::makeFieldView(std::span<real>(d, n * nc),
                                      FieldLayout::POINT_BY_POINT, n, nc);
  const auto all = FieldView<real>(sc, 1, sc.size(), sc.size(), 1);
  ops.fill(all, -1);  // contiguous fill
  ops.copy(vc, va);   // transposition
  ops.copy(vd, vc);   // transposition
  backend.synchronize();
  check(backend.toHost(d, n * nc) == values, "invalid copies");
  const auto pc = backend.toHost(c, sc.size());
  check((pc[n] == -1) && (pc[2 * (n + 5) + n + 4] == -1) &&
            (pc[(n + 5) + 2] == values[2 * nc + 1]),
        "invalid component by component copy");
  ops.fill(vd, 0);   // contiguous fill with zero
  ops.copy(va, vd);  // contiguous copy
  ops.fill(mgis::getComponents(vc, 1, 1), 2.5);  // strided fill
  backend.synchronize();
  check(backend.toHost(a, n * nc) == std::vector<real>(n * nc, 0),
        "invalid fill");
  const auto pc2 = backend.toHost(c, sc.size());
  check((pc2[(n + 5) + 3] == 2.5) && (pc2[3] == values[3 * nc]),
        "invalid fill of a component");
#if !defined(MGIS_CHECK_OPENMP)
  // memory resource: memory written by the device, read by the host
  {
#if defined(MGIS_CHECK_CUDA) || defined(MGIS_CHECK_STLPAR)
    auto r = mgis::memory_resources::CUDAManagedMemoryResource{};
#elif defined(MGIS_CHECK_RAJA)
    // managed memory of the default resource of camp (CUDA)
    auto r = mgis::memory_resources::RAJAManagedMemoryResource<>{};
#elif defined(MGIS_CHECK_HIP)
    auto r = mgis::memory_resources::HIPManagedMemoryResource{};
#elif defined(MGIS_CHECK_SYCL)
    auto r = mgis::memory_resources::SYCLSharedMemoryResource{backend.queue};
#elif defined(MGIS_CHECK_KOKKOS)
    auto r = mgis::memory_resources::KokkosSharedSpaceMemoryResource{};
#endif
    auto v = std::pmr::vector<real>(100, real{1}, &r);
    ops.fill(mgis::makeFieldView(std::span<real>(v),
                                 FieldLayout::POINT_BY_POINT, 100, 1),
             3);
    backend.synchronize();
    check(v == std::pmr::vector<real>(100, real{3}),
          "invalid memory resource");
    // allocations aligned up to and beyond the alignment guaranteed by the
    // memory of the backend (stricter ones being requested by the pools of
    // the standard library)
    for (const auto al : {std::size_t{32}, std::size_t{64}, std::size_t{256},
                          std::size_t{512}, std::size_t{4096}}) {
      auto* const q = r.allocate(1000, al);
      check(reinterpret_cast<std::uintptr_t>(q) % al == 0,
            "invalid alignment of the memory resource");
      r.deallocate(q, 1000, al);
    }
    // a memory resource returns an address even for zero bytes
    auto* const z = r.allocate(0, alignof(real));
    check(z != nullptr, "allocation of zero bytes: null pointer returned");
    r.deallocate(z, 0, alignof(real));
    auto pool = std::pmr::synchronized_pool_resource(&r);
    auto w = std::pmr::vector<real>(600, real{1}, &pool);
    ops.fill(mgis::makeFieldView(std::span<real>(w),
                                 FieldLayout::POINT_BY_POINT, 600, 1),
             4);
    backend.synchronize();
    check(w == std::pmr::vector<real>(600, real{4}),
          "invalid pool of the memory resource");
  }
#endif
  for (auto* const p : {a, c, d}) {
    backend.free(p);
  }
}

int main() {
#if defined(MGIS_CHECK_KOKKOS)
  Kokkos::initialize();
#endif
  try {
    auto backend = Backend{};
#if defined(MGIS_CHECK_SYCL)
    std::cout << "device: "
              << backend.queue.get_device().get_info<sycl::info::device::name>()
              << '\n';
    check(!backend.queue.get_device().is_cpu(), "the SYCL device is the host");
#elif defined(MGIS_CHECK_OPENMP)
    check(omp_get_num_devices() > 0, "no OpenMP device");
#endif
    // with CUDA, HIP and the parallel algorithms of nvc++, the views are in
    // the memory of the device only (cudaMalloc, hipMalloc): operations made
    // by the host would fail
    test(backend);
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    success = false;
  }
#if defined(MGIS_CHECK_KOKKOS)
  Kokkos::finalize();
#endif
  std::cout << (success ? "success" : "failure") << '\n';
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
