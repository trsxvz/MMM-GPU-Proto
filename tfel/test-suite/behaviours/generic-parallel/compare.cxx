/*!
 * \file   test-suite/behaviours/generic-parallel/compare.cxx
 * \brief  This program compares the results of a behaviour integrated over a
 * batch of integration points with the `generic-parallel` interface to the
 * results of the `generic` interface (one call per point).
 *
 * Two time steps of a random loading are computed: the second one starts from
 * the results of the first one. All the combinations of the following cases
 * are tested:
 *
 * - type of computation: consistent tangent operator (or the type of
 *   tangent operator given by the `--tangent-operator` option), elastic
 *   prediction operator, no tangent operator;
 * - storage of the fields: point by point, component by component.
 *
 * For finite strain behaviours, the stress measure and the tangent operator
 * are given by the `--finite-strain` option (for instance
 * `--finite-strain=PK1,DPK1_DF`), and the deformation gradients are
 * perturbations of the identity.
 *
 * The stored and dissipated energies are compared: with the
 * `--energies=true` option, the behaviour must compute non null energies. With
 * the
 * `--speed-of-sound=<mass density>` option, the speed of sound is requested
 * and compared, the mass density being uniform.
 *
 * The `--perturbation-of-internal-state-variables=<amplitude>` option adds a
 * random perturbation to the initial values of the internal state variables
 * and the `--out-of-bounds-policy=<behaviour>,<None|Warning|Strict>` option
 * sets the out of bounds policy in both libraries: those options allow to
 * test batches where the integration fails at some points.
 *
 * The `--initial-internal-state-variables=<values>` option gives the initial
 * values of the internal state variables of a point (null by default), as a
 * comma separated list, a single value being used for all the internal state
 * variables.
 *
 * With the `--initialize-function=<number of inputs>` option (resp. the
 * `--post-processing=<number of outputs>` option), the tested function is an
 * initialize function (resp. a post-processing) executed on random states,
 * for both storages of the fields. With the `--rotation=<number of values>`
 * option, the tested function is a rotation function of an orthotropic
 * behaviour.
 *
 * This file can be compiled by a C++ compiler, by a CUDA compiler or by a
 * HIP compiler. In the two last cases, the data are allocated in managed
 * memory. The HIP runtime is used if the HIP language is enabled or if the
 * `MFRONT_GPB_COMPARE_USE_HIP` macro is defined (required on NVIDIA GPUs).
 * If the `MFRONT_GPB_COMPARE_USE_KOKKOS` macro is defined, Kokkos is
 * initialized and the data are allocated in its shared memory space (if any).
 * If the `MFRONT_GPB_COMPARE_USE_SYCL` macro is defined, the data are
 * allocated in shared USM memory (or in host memory if the
 * `MFRONT_GPB_COMPARE_USE_HOST_MEMORY` macro is defined) and the
 * computations are performed on the device selected by the default SYCL
 * selector.
 *
 * \date   05/10/2026
 */

#include <dlfcn.h>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>
#include <utility>
#include <type_traits>
#include <algorithm>
#if defined(MFRONT_GPB_COMPARE_USE_KOKKOS)
#include <Kokkos_Core.hpp>
#elif defined(MFRONT_GPB_COMPARE_USE_SYCL)
#include <sycl/sycl.hpp>
#elif defined(__HIP__) || defined(MFRONT_GPB_COMPARE_USE_HIP)
#include <hip/hip_runtime.h>
#define GPU_API(f) hip##f
#elif defined(__CUDACC__)
#include <cuda_runtime.h>
#define GPU_API(f) cuda##f
#endif
#include "MFront/GenericBehaviour/BehaviourData.h"
#include "MFront/GenericParallelBehaviour/BatchData.h"

using GenericFunction = int (*)(mfront_gb_BehaviourData*);
using BatchFunction = int (*)(mfront_gpb_OutputStatus*,
                              const mfront_gpb_BatchData*);

//! \brief sizes of the data of an integration point
struct Sizes {
  std::size_t g, th, mp, isv, esv, K;
};

/*!
 * \brief options of finite strain behaviours, encoded as the second and
 * third values of the `K` array of the `generic` interface
 */
struct FiniteStrainOptions {
  bool finite_strain = false;
  double stress_measure = 0;
  double tangent_operator = 0;
};

//! \brief states of all the points, stored point by point
struct States {
  std::vector<double> g0, g1, th0, th1, mp, isv0, isv1, esv0, esv1, K, rdt;
  //! \brief stored and dissipated energies
  std::vector<double> stored_energy0, stored_energy1, dissipated_energy0,
      dissipated_energy1;
  //! \brief speed of sound
  std::vector<double> speed_of_sound;
  //! \brief uniform mass density, not given to the behaviour if null
  double mass_density = 0;
  std::vector<int> status;
  States(const std::size_t n, const Sizes& s)
      : g0(n * s.g),
        g1(n * s.g),
        th0(n * s.th),
        th1(n * s.th),
        mp(n * s.mp),
        isv0(n * s.isv),
        isv1(n * s.isv),
        esv0(n * s.esv),
        esv1(n * s.esv),
        K(n * s.K),
        rdt(n, 1),
        stored_energy0(n, 0),
        stored_energy1(n, 0),
        dissipated_energy0(n, 0),
        dissipated_energy1(n, 0),
        speed_of_sound(n, 0),
        status(n, 0) {}
};

[[noreturn]] static void fail(const std::string& msg) {
  std::fprintf(stderr, "%s\n", msg.c_str());
  std::exit(EXIT_FAILURE);
}

#if defined(MFRONT_GPB_COMPARE_USE_KOKKOS)

//! \brief initialization of Kokkos
struct Runtime {
  Kokkos::ScopeGuard guard;
};

//! \brief memory space of the data
using MemorySpace = std::conditional_t<Kokkos::has_shared_space,
                                       Kokkos::SharedSpace,
                                       Kokkos::HostSpace>;

template <typename T>
static T* allocate(const std::size_t n) {
  return static_cast<T*>(Kokkos::kokkos_malloc<MemorySpace>(
      std::max<std::size_t>(n, 1) * sizeof(T)));
}

static void deallocate(void* const p) { Kokkos::kokkos_free<MemorySpace>(p); }

template <typename T>
static void prefetch(T* const, const std::size_t) {}

static void synchronize() { Kokkos::fence(); }

//! \return the queue passed to the behaviour (default execution space)
static void* queue() { return nullptr; }

#elif defined(MFRONT_GPB_COMPARE_USE_SYCL)

//! \brief nothing to initialize
struct Runtime {};

//! \return the SYCL queue used by all the computations
static ::sycl::queue& getQueue() {
  static auto q = ::sycl::queue{};
  return q;
}

#if defined(MFRONT_GPB_COMPARE_USE_HOST_MEMORY)

template <typename T>
static T* allocate(const std::size_t n) {
  auto* const p =
      static_cast<T*>(std::malloc(std::max<std::size_t>(n, 1) * sizeof(T)));
  if (p == nullptr) {
    fail("allocation failed");
  }
  return p;
}

static void deallocate(void* const p) { std::free(p); }

template <typename T>
static void prefetch(T* const, const std::size_t) {}

#else /* defined(MFRONT_GPB_COMPARE_USE_HOST_MEMORY) */

template <typename T>
static T* allocate(const std::size_t n) {
  auto* const p =
      ::sycl::malloc_shared<T>(std::max<std::size_t>(n, 1), getQueue());
  if (p == nullptr) {
    fail("allocation failed");
  }
  return p;
}

static void deallocate(void* const p) { ::sycl::free(p, getQueue()); }

template <typename T>
static void prefetch(T* const p, const std::size_t n) {
  if (n != 0) {
    getQueue().prefetch(p, n * sizeof(T));
  }
}

#endif /* defined(MFRONT_GPB_COMPARE_USE_HOST_MEMORY) */

static void synchronize() { getQueue().wait_and_throw(); }

//! \return the queue passed to the behaviour
static void* queue() { return &getQueue(); }

#elif defined(GPU_API)

//! \brief nothing to initialize
struct Runtime {};

static void check(const GPU_API(Error_t) e) {
  if (e != GPU_API(Success)) {
    fail(std::string{"GPU error: "} + GPU_API(GetErrorString)(e));
  }
}

template <typename T>
static T* allocate(const std::size_t n) {
  T* p = nullptr;
  check(GPU_API(MallocManaged)(&p, std::max<std::size_t>(n, 1) * sizeof(T)));
  return p;
}

static void deallocate(void* const p) { check(GPU_API(Free)(p)); }

//! \brief move the data to the device before the computation
template <typename T>
static void prefetch(T* const p, const std::size_t n) {
  if (n == 0) {
    return;
  }
  int device = 0;
  check(GPU_API(GetDevice)(&device));
#if defined(__HIPCC__) || defined(__HIP__)
  check(hipMemPrefetchAsync(p, n * sizeof(T), device, 0));
#else
  auto location = cudaMemLocation{};
  location.type = cudaMemLocationTypeDevice;
  location.id = device;
  check(cudaMemPrefetchAsync(p, n * sizeof(T), location, 0, 0));
#endif
}

static void synchronize() { check(GPU_API(DeviceSynchronize)()); }

//! \return the queue passed to the behaviour (default stream)
static void* queue() { return nullptr; }

#else

//! \brief nothing to initialize
struct Runtime {};

template <typename T>
static T* allocate(const std::size_t n) {
  auto* const p =
      static_cast<T*>(std::malloc(std::max<std::size_t>(n, 1) * sizeof(T)));
  if (p == nullptr) {
    fail("allocation failed");
  }
  return p;
}

static void deallocate(void* const p) { std::free(p); }

template <typename T>
static void prefetch(T* const, const std::size_t) {}

static void synchronize() {}

//! \return the queue passed to the behaviour (not used)
static void* queue() { return nullptr; }

#endif

static void* load(const std::string& lib, const std::string& f) {
  auto* const l = dlopen(lib.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (l == nullptr) {
    fail(dlerror());
  }
  auto* const s = dlsym(l, f.c_str());
  if (s == nullptr) {
    fail(dlerror());
  }
  return s;
}

/*!
 * \brief call the given function at each point with the behaviour data of
 * the `generic` interface
 * \param[in,out] st: states
 * \param[in] s: sizes of the data of a point
 * \param[in] fs: finite strain options
 * \param[in] n: number of points
 * \param[in] dt: time increment
 * \param[in] K0: type of computation
 * \param[in] f: function taking the behaviour data and the index of the
 * point, and returning the status
 */
template <typename Function>
static void forEachReferencePoint(States& st,
                                  const Sizes& s,
                                  const FiniteStrainOptions& fs,
                                  const std::size_t n,
                                  const double dt,
                                  const double K0,
                                  const Function& f) {
  for (std::size_t i = 0; i != n; ++i) {
    auto d = mfront_gb_BehaviourData{};
    d.error_message = nullptr;
    d.dt = dt;
    // the type of computation is stored in the tangent operator, which is
    // empty for models
    auto K = std::array<double, 3>{K0, fs.stress_measure, fs.tangent_operator};
    if (s.K != 0) {
      std::copy_n(K.begin(), fs.finite_strain ? 3 : 1, st.K.begin() + i * s.K);
    }
    d.K = s.K != 0 ? st.K.data() + i * s.K : K.data();
    st.rdt[i] = 1;
    d.rdt = st.rdt.data() + i;
    d.speed_of_sound = st.speed_of_sound.data() + i;
    d.s0.gradients = st.g0.data() + i * s.g;
    d.s1.gradients = st.g1.data() + i * s.g;
    d.s0.thermodynamic_forces = st.th0.data() + i * s.th;
    d.s1.thermodynamic_forces = st.th1.data() + i * s.th;
    d.s0.material_properties = st.mp.data() + i * s.mp;
    d.s1.material_properties = st.mp.data() + i * s.mp;
    d.s0.internal_state_variables = st.isv0.data() + i * s.isv;
    d.s1.internal_state_variables = st.isv1.data() + i * s.isv;
    d.s0.external_state_variables = st.esv0.data() + i * s.esv;
    d.s1.external_state_variables = st.esv1.data() + i * s.esv;
    d.s0.stored_energy = st.stored_energy0.data() + i;
    d.s1.stored_energy = st.stored_energy1.data() + i;
    d.s0.dissipated_energy = st.dissipated_energy0.data() + i;
    d.s1.dissipated_energy = st.dissipated_energy1.data() + i;
    d.s0.mass_density = d.s1.mass_density =
        st.mass_density > 0 ? &(st.mass_density) : nullptr;
    st.status[i] = f(d, i);
  }
}  // end of forEachReferencePoint

//! \brief reference: one call to the `generic` interface per point
static void integrateReference(States& st,
                               const Sizes& s,
                               const FiniteStrainOptions& fs,
                               const std::size_t n,
                               const double dt,
                               const double K0,
                               GenericFunction f) {
  forEachReferencePoint(
      st, s, fs, n, dt, K0,
      [f](mfront_gb_BehaviourData& d, const std::size_t) { return f(&d); });
}  // end of integrateReference

/*!
 * \brief a field of the batch, stored point by point or component by
 * component
 */
/*!
 * \brief storages of the fields of a batch:
 *
 * - `by_point`: the values of a point are contiguous;
 * - `by_component`: the values of a component are contiguous;
 * - `by_variable`: as `by_component`, but the material properties and the
 *   external state variables are given variable by variable, with uniform
 *   values (null point stride). This storage assumes that the material
 *   properties and the external state variables are scalars and uniform,
 *   as in all the tests.
 */
enum struct Storage { by_point, by_component, by_variable };

//! \return a description of a storage
static const char* toString(const Storage s) {
  if (s == Storage::by_point) {
    return "by point";
  }
  return s == Storage::by_component ? "by component" : "by variable";
}

//! \brief all the storages
static constexpr Storage storages[] = {Storage::by_point, Storage::by_component,
                                       Storage::by_variable};

struct BatchField {
  BatchField(const std::vector<double>& v,
             const std::size_t n_,
             const std::size_t c_,
             const bool by_component_)
      : values(allocate<double>(v.size())),
        n(n_),
        c(c_),
        by_component(by_component_) {
    for (std::size_t i = 0; i != this->n; ++i) {
      for (std::size_t j = 0; j != this->c; ++j) {
        this->values[this->index(i, j)] = v[i * this->c + j];
      }
    }
    prefetch(this->values, this->n * this->c);
  }
  BatchField(const BatchField&) = delete;
  BatchField& operator=(const BatchField&) = delete;
  ~BatchField() { deallocate(this->values); }
  //! \brief copy the values in a vector storing them point by point
  void copy(std::vector<double>& v) const {
    for (std::size_t i = 0; i != this->n; ++i) {
      for (std::size_t j = 0; j != this->c; ++j) {
        v[i * this->c + j] = this->values[this->index(i, j)];
      }
    }
  }
  mfront_gpb_Field field() const {
    return this->by_component ? mfront_gpb_Field{this->values, 1, this->n}
                              : mfront_gpb_Field{this->values, this->c, 1};
  }
  mfront_gpb_ConstField const_field() const {
    const auto f = this->field();
    return {f.values, f.point_stride, f.component_stride};
  }
  //! \return a field giving the values of the `j`-th component of the first
  //! point at all points
  mfront_gpb_ConstField uniform(const std::size_t j) const {
    return {this->values + this->index(0, j), 0, 1};
  }

 private:
  std::size_t index(const std::size_t i, const std::size_t j) const {
    return this->by_component ? j * this->n + i : i * this->c + j;
  }
  double* const values;
  const std::size_t n, c;
  const bool by_component;
};

/*!
 * \brief fields of a batch, initialized from the states of the points.
 * The batch data refer to those fields.
 */
struct BatchFields {
  BatchFields(const States& st,
              const Sizes& s,
              const FiniteStrainOptions& fs,
              const std::size_t n,
              const double dt,
              const double K0,
              const Storage storage)
      : by_component(storage != Storage::by_point),
        g0(st.g0, n, s.g, by_component),
        g1(st.g1, n, s.g, by_component),
        th0(st.th0, n, s.th, by_component),
        th1(st.th1, n, s.th, by_component),
        mp(st.mp, n, s.mp, by_component),
        isv0(st.isv0, n, s.isv, by_component),
        isv1(st.isv1, n, s.isv, by_component),
        esv0(st.esv0, n, s.esv, by_component),
        esv1(st.esv1, n, s.esv, by_component),
        K(st.K, n, s.K, by_component),
        rdt(st.rdt, n, 1, false),
        stored_energy0(st.stored_energy0, n, 1, by_component),
        stored_energy1(st.stored_energy1, n, 1, by_component),
        dissipated_energy0(st.dissipated_energy0, n, 1, by_component),
        dissipated_energy1(st.dissipated_energy1, n, 1, by_component),
        speed_of_sound(st.speed_of_sound, n, 1, by_component),
        mass_density(std::vector<double>{st.mass_density}, 1, 1, false),
        status(allocate<int>(n)),
        number_of_points(n) {
    prefetch(this->status, n);
    synchronize();
    auto& d = this->data;
    d.queue = queue();
    d.number_of_points = n;
    d.dt = dt;
    d.options[0] = K0;
    d.options[1] = fs.stress_measure;
    d.options[2] = fs.tangent_operator;
    d.rdt = 1;
    d.tangent_operators = this->K.field();
    d.rdts = this->rdt.field();
    d.statuses = this->status;
    d.s0.gradients = this->g0.const_field();
    d.s1.gradients = this->g1.const_field();
    d.s0.thermodynamic_forces = this->th0.const_field();
    d.s1.thermodynamic_forces = this->th1.field();
    d.s0.material_properties = this->mp.const_field();
    d.s1.material_properties = this->mp.const_field();
    d.s0.internal_state_variables = this->isv0.const_field();
    d.s1.internal_state_variables = this->isv1.field();
    d.s0.external_state_variables = this->esv0.const_field();
    d.s1.external_state_variables = this->esv1.const_field();
    d.s0.stored_energy = this->stored_energy0.const_field();
    d.s1.stored_energy = this->stored_energy1.field();
    d.s0.dissipated_energy = this->dissipated_energy0.const_field();
    d.s1.dissipated_energy = this->dissipated_energy1.field();
    d.speed_of_sound = this->speed_of_sound.field();
    if (storage == Storage::by_variable) {
      for (std::size_t j = 0; j != s.mp; ++j) {
        this->mps.push_back(this->mp.uniform(j));
      }
      for (std::size_t j = 0; j != s.esv; ++j) {
        this->esvs0.push_back(this->esv0.uniform(j));
        this->esvs1.push_back(this->esv1.uniform(j));
      }
      // the group fields are ignored
      d.s0.material_properties = d.s1.material_properties = {};
      d.s0.external_state_variables = d.s1.external_state_variables = {};
      d.s0.material_properties_by_variable = this->mps.data();
      d.s1.material_properties_by_variable = this->mps.data();
      d.s0.external_state_variables_by_variable = this->esvs0.data();
      d.s1.external_state_variables_by_variable = this->esvs1.data();
    }
    if (st.mass_density > 0) {
      // uniform mass density
      d.s0.mass_density = d.s1.mass_density = this->mass_density.const_field();
      d.s0.mass_density.point_stride = d.s1.mass_density.point_stride = 0;
    }
  }
  BatchFields(const BatchFields&) = delete;
  BatchFields& operator=(const BatchFields&) = delete;
  ~BatchFields() { deallocate(this->status); }
  //! \brief copy the outputs in the states
  void copyResults(States& st) const {
    this->th1.copy(st.th1);
    this->isv1.copy(st.isv1);
    this->K.copy(st.K);
    this->rdt.copy(st.rdt);
    this->stored_energy1.copy(st.stored_energy1);
    this->dissipated_energy1.copy(st.dissipated_energy1);
    this->speed_of_sound.copy(st.speed_of_sound);
    std::copy(this->status, this->status + this->number_of_points,
              st.status.begin());
  }
  const bool by_component;
  const BatchField g0, g1, th0, th1, mp, isv0, isv1, esv0, esv1, K, rdt,
      stored_energy0, stored_energy1, dissipated_energy0, dissipated_energy1,
      speed_of_sound, mass_density;
  //! \brief fields of the material properties, variable by variable
  std::vector<mfront_gpb_ConstField> mps;
  //! \brief fields of the external state variables, variable by variable
  std::vector<mfront_gpb_ConstField> esvs0, esvs1;
  int* const status;
  const std::size_t number_of_points;
  mfront_gpb_BatchData data = {};
};

/*!
 * \brief call the given function and measure its execution time
 * \return the output status
 * \param[out] ms: execution time in milliseconds
 * \param[in] f: function taking an output status
 */
template <typename Function>
static mfront_gpb_OutputStatus measure(double& ms, const Function& f) {
  auto o = mfront_gpb_OutputStatus{};
  const auto start = std::chrono::steady_clock::now();
  f(o);
  const auto end = std::chrono::steady_clock::now();
  ms = std::chrono::duration<double, std::milli>(end - start).count();
  return o;
}  // end of measure

//! \brief batch integration with the `generic-parallel` interface
static mfront_gpb_OutputStatus integrateBatch(States& st,
                                              const Sizes& s,
                                              const FiniteStrainOptions& fs,
                                              const std::size_t n,
                                              const double dt,
                                              const double K0,
                                              const Storage storage,
                                              BatchFunction f,
                                              double& ms) {
  const auto b = BatchFields(st, s, fs, n, dt, K0, storage);
  const auto o =
      measure(ms, [&b, f](mfront_gpb_OutputStatus& r) { f(&r, &b.data); });
  b.copyResults(st);
  return o;
}

/*!
 * \brief tolerance on the relative differences. Both libraries are not
 * compiled by the same compiler, which may contract the floating point
 * operations differently: the results agree up to round-off errors, amplified
 * by the local Newton solvers of the behaviours. A defect gives differences of
 * the order of one.
 */
static constexpr double tolerance = 1e-10;

//! \brief relative difference over the points where the integration succeeded
static double relativeDifference(const std::vector<double>& a,
                                 const std::vector<double>& b,
                                 const std::vector<int>& status,
                                 const std::size_t c) {
  auto m = 0., r = 0.;
  for (std::size_t i = 0; i != status.size(); ++i) {
    if (status[i] == -1) {
      continue;
    }
    for (std::size_t j = 0; j != c; ++j) {
      m = std::max(m, std::abs(a[i * c + j] - b[i * c + j]));
      r = std::max(r, std::abs(a[i * c + j]));
    }
  }
  return r == 0 ? m : m / r;
}

//! \return the values of a comma separated list, `-` denoting an empty list
static std::vector<double> parseList(const std::string& l) {
  auto r = std::vector<double>{};
  if (l == "-") {
    return r;
  }
  for (auto v = l; !v.empty();) {
    const auto c = v.find(',');
    r.push_back(std::stod(v.substr(0, c)));
    v = c == std::string::npos ? "" : v.substr(c + 1);
  }
  return r;
}

/*!
 * \return the options of finite strain behaviours
 * \param[in] o: stress measure and tangent operator, separated by a comma
 */
static FiniteStrainOptions parseFiniteStrainOptions(const std::string& o) {
  const auto p = o.find(',');
  if (p == std::string::npos) {
    fail("invalid finite strain options '" + o + "'");
  }
  const auto sm = o.substr(0, p);
  const auto to = o.substr(p + 1);
  auto r = FiniteStrainOptions{.finite_strain = true};
  if (sm == "CAUCHY") {
    r.stress_measure = 0;
  } else if (sm == "PK2") {
    r.stress_measure = 1;
  } else if (sm == "PK1") {
    r.stress_measure = 2;
  } else {
    fail("invalid stress measure '" + sm + "'");
  }
  if (to == "DSIG_DF") {
    r.tangent_operator = 0;
  } else if (to == "DS_DEGL") {
    r.tangent_operator = 1;
  } else if (to == "DPK1_DF") {
    r.tangent_operator = 2;
  } else if (to == "DTAU_DDF") {
    r.tangent_operator = 3;
  } else {
    fail("invalid tangent operator '" + to + "'");
  }
  return r;
}

/*!
 * \return true if the summary of the integration of a batch is consistent
 * with the statuses of the points
 * \param[in] o: output status
 * \param[in] status: statuses of the points
 */
static bool isConsistent(const mfront_gpb_OutputStatus& o,
                         const std::vector<int>& status) {
  const auto nfailures =
      static_cast<std::size_t>(std::count(status.begin(), status.end(), -1));
  const auto nunreliable =
      static_cast<std::size_t>(std::count(status.begin(), status.end(), 0));
  const auto expected = nfailures != 0 ? -1 : (nunreliable != 0 ? 0 : 1);
  return (o.status == expected) && (o.number_of_failures == nfailures) &&
         (o.number_of_unreliable_results == nunreliable);
}

//! \brief function tested
enum struct Mode {
  INTEGRATION,
  INITIALIZE_FUNCTION,
  POST_PROCESSING,
  ROTATION
};

//! \brief optional arguments
struct Options {
  //! \brief options of finite strain behaviours
  FiniteStrainOptions finite_strain;
  //! \brief function tested
  Mode mode = Mode::INTEGRATION;
  //! \brief number of inputs of an initialize function, number of outputs
  //! of a post-processing or number of rotated values of a point
  std::size_t number_of_values = 0;
  /*!
   * \brief initial values of the internal state variables of a point (a
   * single value is used for all the internal state variables)
   */
  std::vector<double> initial_internal_state_variables = {0};
  //! \brief mass density, if the speed of sound is requested
  double mass_density = 0;
  //! \brief if true, the behaviour must compute non null energies
  bool energies = false;
  //! \brief type of the tangent operator (consistent tangent operator by
  //! default)
  double tangent_operator = 4;
  //! \brief amplitude of the random perturbation of the internal state
  //! variables
  double internal_state_variables_perturbation = 0;
  //! \brief name of the behaviour, required to set the out of bounds policy
  std::string behaviour;
  //! \brief out of bounds policy (0: None, 1: Warning, 2: Strict), if
  //! positive
  int out_of_bounds_policy = -1;
};

//! \return the optional arguments
static Options parseOptions(const int argc,
                            const char* const* const argv,
                            const int first) {
  auto r = Options{};
  for (auto i = first; i != argc; ++i) {
    const auto a = std::string{argv[i]};
    const auto p = a.find('=');
    if (p == std::string::npos) {
      fail("invalid option '" + a + "'");
    }
    const auto o = a.substr(0, p);
    const auto v = a.substr(p + 1);
    if (o == "--finite-strain") {
      r.finite_strain = parseFiniteStrainOptions(v);
    } else if (o == "--initialize-function") {
      r.mode = Mode::INITIALIZE_FUNCTION;
      r.number_of_values = std::stoul(v);
    } else if (o == "--post-processing") {
      r.mode = Mode::POST_PROCESSING;
      r.number_of_values = std::stoul(v);
    } else if (o == "--perturbation-of-internal-state-variables") {
      r.internal_state_variables_perturbation = std::stod(v);
    } else if (o == "--out-of-bounds-policy") {
      // name of the behaviour and policy, separated by a comma
      const auto c = v.find(',');
      if (c == std::string::npos) {
        fail("invalid out of bounds policy '" + v + "'");
      }
      r.behaviour = v.substr(0, c);
      const auto policy = v.substr(c + 1);
      r.out_of_bounds_policy =
          policy == "None" ? 0 : (policy == "Warning" ? 1 : 2);
      if ((policy != "None") && (policy != "Warning") && (policy != "Strict")) {
        fail("invalid out of bounds policy '" + policy + "'");
      }
    } else if (o == "--tangent-operator") {
      r.tangent_operator = std::stod(v);
    } else if (o == "--energies") {
      r.energies = v == "true";
    } else if (o == "--speed-of-sound") {
      r.mass_density = std::stod(v);
    } else if (o == "--initial-internal-state-variables") {
      r.initial_internal_state_variables = parseList(v);
    } else if (o == "--rotation") {
      r.mode = Mode::ROTATION;
      r.number_of_values = std::stoul(v);
    } else {
      fail("invalid option '" + o + "'");
    }
  }
  return r;
}

/*!
 * \brief description of a test: libraries, function, sizes, loading and
 * options
 */
struct Test {
  std::string generic_library, batch_library, function;
  Sizes sizes;
  Options options;
  //! \brief initial states of the points
  States initial;
  //! \brief identity, or null values for small strain behaviours
  std::vector<double> identity;
  //! \brief number of points
  std::size_t n;
  //! \brief time increment
  double dt;
};

//! \brief print the results of a comparison
static bool print(const char* const what,
                  const Storage storage,
                  const bool pass,
                  const mfront_gpb_OutputStatus& o,
                  const std::vector<int>& status,
                  const double d,
                  const double ms) {
  std::printf(
      "%-12s %-12s status=%2d failures=%6ld | relative difference %.1e | "
      "%.2f ms | %s\n",
      what, toString(storage), o.status,
      static_cast<long>(std::count(status.begin(), status.end(), -1)), d, ms,
      pass ? "OK" : "FAILED");
  return pass;
}

/*!
 * \brief compare two time steps integrated by both interfaces, for all
 * types of computation and storages of the fields
 */
static bool testIntegration(const Test& t) {
  const auto fg =
      reinterpret_cast<GenericFunction>(load(t.generic_library, t.function));
  const auto fb =
      reinterpret_cast<BatchFunction>(load(t.batch_library, t.function));
  const auto& s = t.sizes;
  const auto& fs = t.options.finite_strain;
  const auto n = t.n;
  // the speed of sound is requested by an offset of the type of computation
  const auto speed_of_sound = t.initial.mass_density > 0;
  const auto o0 = speed_of_sound ? 100. : 0.;
  auto ok = true;
  const auto Kt = t.options.tangent_operator;
  for (const auto K0 : {Kt, -1., 0.}) {
    // models have no tangent operator, hence no prediction
    if ((s.K == 0) && (K0 < -0.25)) {
      continue;
    }
    for (const auto storage : storages) {
      // first step from the initial state, second step from its results
      auto ref = t.initial;
      auto batch = t.initial;
      auto ms1 = 0., ms2 = 0.;
      integrateReference(ref, s, fs, n, t.dt, Kt + o0, fg);
      const auto o1 =
          integrateBatch(batch, s, fs, n, t.dt, Kt + o0, storage, fb, ms1);
      // the test is meaningless if the integration fails at all points
      const auto pass1 = (ref.status == batch.status) &&
                         isConsistent(o1, batch.status) &&
                         (o1.number_of_failures != n);
      for (auto* const st : {&ref, &batch}) {
        st->g0 = st->g1;
        st->th0 = st->th1;
        st->isv0 = st->isv1;
        st->stored_energy0 = st->stored_energy1;
        st->dissipated_energy0 = st->dissipated_energy1;
        // the perturbation of the gradients is increased by 50%
        for (std::size_t i = 0; i != n; ++i) {
          for (std::size_t j = 0; j != s.g; ++j) {
            auto& g = st->g1[i * s.g + j];
            g = t.identity[j] + 1.5 * (g - t.identity[j]);
          }
        }
      }
      integrateReference(ref, s, fs, n, t.dt, K0 + o0, fg);
      const auto o2 =
          integrateBatch(batch, s, fs, n, t.dt, K0 + o0, storage, fb, ms2);
      const auto nfailures =
          std::count(ref.status.begin(), ref.status.end(), -1);
      const auto dth = relativeDifference(ref.th1, batch.th1, ref.status, s.th);
      const auto disv =
          relativeDifference(ref.isv1, batch.isv1, ref.status, s.isv);
      // the tangent operator is not computed if K0 is in [-0.25, 0.5]
      const auto tangent = (K0 < -0.25) || (K0 > 0.5);
      const auto dK =
          tangent ? relativeDifference(ref.K, batch.K, ref.status, s.K) : 0.;
      // energies (null if the behaviour does not compute them) and speed
      // of sound
      const auto de =
          std::max(relativeDifference(ref.stored_energy1, batch.stored_energy1,
                                      ref.status, 1),
                   relativeDifference(ref.dissipated_energy1,
                                      batch.dissipated_energy1, ref.status, 1));
      const auto dsos =
          speed_of_sound
              ? relativeDifference(ref.speed_of_sound, batch.speed_of_sound,
                                   ref.status, 1)
              : 0.;
      // the comparison of the energies and of the speed of sound is
      // meaningless if they are not computed
      auto computed = [&ref](const std::vector<double>& v) {
        for (std::size_t i = 0; i != v.size(); ++i) {
          if ((ref.status[i] != -1) && (v[i] != 0)) {
            return true;
          }
        }
        return false;
      };
      // the energies and the speed of sound are only computed after an
      // integration
      const auto integration = K0 > -0.25;
      const auto energies_computed =
          (!integration) || (!t.options.energies) ||
          (computed(ref.stored_energy1) && computed(ref.dissipated_energy1));
      const auto speed_of_sound_computed =
          (!integration) || (!speed_of_sound) || computed(ref.speed_of_sound);
      const auto same_status = ref.status == batch.status;
      const auto pass = pass1 && same_status && energies_computed &&
                        speed_of_sound_computed &&
                        isConsistent(o2, batch.status) && (dth < tolerance) &&
                        (disv < tolerance) && (dK < tolerance) &&
                        (de < tolerance) && (dsos < tolerance);
      ok = ok && pass;
      std::printf(
          "K0=%4.1f %-12s status=%2d/%2d failures=%6ld same statuses=%d | "
          "relative differences: forces %.1e, isvs %.1e, K %.1e, energies "
          "%.1e, speed of sound %.1e | %.2f ms, %.2f ms | %s\n",
          K0, toString(storage), o1.status, o2.status,
          static_cast<long>(nfailures), same_status, dth, disv, dK, de, dsos,
          ms1, ms2, pass ? "OK" : "FAILED");
    }
  }
  return ok;
}  // end of testIntegration

/*!
 * \brief the functions of the `generic` interface executing an initialize
 * function or a post-processing return 0 on success, whereas the functions
 * of the `generic-parallel` interface report the same statuses as the
 * integration (1 on success).
 */
static void convertReferenceStatuses(std::vector<int>& status) {
  std::replace(status.begin(), status.end(), 0, 1);
}

/*!
 * \brief compare the execution of an initialize function by both
 * interfaces, for both storages of the fields
 */
static bool testInitializeFunction(const Test& t) {
  using GenericInitializeFunction =
      int (*)(mfront_gb_BehaviourData*, const mfront_gb_real*);
  using BatchInitializeFunction =
      int (*)(mfront_gpb_OutputStatus*, const mfront_gpb_BatchData*,
              const mfront_gpb_ConstField*);
  const auto fg = reinterpret_cast<GenericInitializeFunction>(
      load(t.generic_library, t.function));
  const auto fb = reinterpret_cast<BatchInitializeFunction>(
      load(t.batch_library, t.function));
  const auto& s = t.sizes;
  const auto ni = t.options.number_of_values;
  auto rng = std::mt19937{1};
  auto u = std::uniform_real_distribution<double>{-1e8, 1e8};
  auto inputs = std::vector<double>(t.n * ni);
  std::generate(inputs.begin(), inputs.end(), [&] { return u(rng); });
  auto ok = true;
  for (const auto storage : storages) {
    auto ref = t.initial;
    auto batch = t.initial;
    forEachReferencePoint(
        ref, s, t.options.finite_strain, t.n, t.dt, 0,
        [fg, &inputs, ni](mfront_gb_BehaviourData& d, const std::size_t i) {
          return fg(&d, ni == 0 ? nullptr : inputs.data() + i * ni);
        });
    convertReferenceStatuses(ref.status);
    auto ms = 0.;
    const auto b = BatchFields(batch, s, t.options.finite_strain, t.n, t.dt, 0,
                               storage);
    const auto in = BatchField(inputs, t.n, ni, b.by_component);
    const auto in_field = in.const_field();
    const auto o = measure(ms, [&](mfront_gpb_OutputStatus& r) {
      // no field is given if the function has no input
      fb(&r, &b.data, ni == 0 ? nullptr : &in_field);
    });
    b.copyResults(batch);
    const auto d =
        std::max(relativeDifference(ref.th1, batch.th1, ref.status, s.th),
                 relativeDifference(ref.isv1, batch.isv1, ref.status, s.isv));
    // the test is meaningless if the initialize function changes nothing
    const auto changed =
        (ref.th1 != t.initial.th1) || (ref.isv1 != t.initial.isv1);
    const auto pass = changed && (ref.status == batch.status) &&
                      isConsistent(o, batch.status) && (d < tolerance);
    ok = print("initialize", storage, pass, o, batch.status, d, ms) && ok;
  }
  return ok;
}  // end of testInitializeFunction

/*!
 * \brief compare the execution of a post-processing by both interfaces, for
 * both storages of the fields
 */
static bool testPostProcessing(const Test& t) {
  using GenericPostProcessing =
      int (*)(mfront_gb_real*, mfront_gb_BehaviourData*);
  using BatchPostProcessing =
      int (*)(mfront_gpb_OutputStatus*, const mfront_gpb_Field*,
              const mfront_gpb_BatchData*);
  const auto fg = reinterpret_cast<GenericPostProcessing>(
      load(t.generic_library, t.function));
  const auto fb =
      reinterpret_cast<BatchPostProcessing>(load(t.batch_library, t.function));
  const auto& s = t.sizes;
  const auto no = t.options.number_of_values;
  auto ok = true;
  for (const auto storage : storages) {
    auto ref = t.initial;
    auto batch = t.initial;
    auto ref_outputs = std::vector<double>(t.n * no);
    forEachReferencePoint(ref, s, t.options.finite_strain, t.n, t.dt, 0,
                          [fg, &ref_outputs, no](mfront_gb_BehaviourData& d,
                                                 const std::size_t i) {
                            return fg(ref_outputs.data() + i * no, &d);
                          });
    convertReferenceStatuses(ref.status);
    auto ms = 0.;
    const auto b = BatchFields(batch, s, t.options.finite_strain, t.n, t.dt, 0,
                               storage);
    const auto out =
        BatchField(std::vector<double>(t.n * no), t.n, no, b.by_component);
    const auto out_field = out.field();
    const auto o = measure(
        ms, [&](mfront_gpb_OutputStatus& r) { fb(&r, &out_field, &b.data); });
    b.copyResults(batch);
    auto batch_outputs = std::vector<double>(t.n * no);
    out.copy(batch_outputs);
    const auto d =
        relativeDifference(ref_outputs, batch_outputs, ref.status, no);
    // the test is meaningless if all the outputs are null
    const auto computed = std::any_of(ref_outputs.begin(), ref_outputs.end(),
                                      [](const double v) { return v != 0; });
    const auto pass = computed && (ref.status == batch.status) &&
                      isConsistent(o, batch.status) && (d < tolerance);
    ok = print("post-processing", storage, pass, o, batch.status, d, ms) &&
         ok;
  }
  return ok;
}  // end of testPostProcessing

/*!
 * \brief compare the rotation of the values of the points of a batch by both
 * interfaces, for both storages of the fields, with a rotation matrix per
 * point or a uniform rotation matrix
 */
static bool testRotation(const Test& t) {
  using GenericRotation =
      void (*)(mfront_gb_real*, const mfront_gb_real*, const mfront_gb_real*);
  using BatchRotation =
      int (*)(mfront_gpb_OutputStatus*, const mfront_gpb_RotationData*);
  const auto fg =
      reinterpret_cast<GenericRotation>(load(t.generic_library, t.function));
  const auto fb =
      reinterpret_cast<BatchRotation>(load(t.batch_library, t.function));
  const auto n = t.n;
  const auto ns = t.options.number_of_values;
  auto rng = std::mt19937{3};
  auto u = std::uniform_real_distribution<double>{-1, 1};
  auto src = std::vector<double>(n * ns);
  std::generate(src.begin(), src.end(), [&] { return u(rng); });
  // rotation matrices, stored by rows: rotations around the z axis, then
  // around the x axis, then around the z axis
  auto rotations = std::vector<double>(9 * n);
  for (std::size_t i = 0; i != n; ++i) {
    const auto a = 3.14 * u(rng), b = 3.14 * u(rng), c = 3.14 * u(rng);
    const auto ca = std::cos(a), sa = std::sin(a), cb = std::cos(b),
               sb = std::sin(b), cc = std::cos(c), sc = std::sin(c);
    const double r[9] = {ca * cc - sa * cb * sc,
                         -ca * sc - sa * cb * cc,
                         sa * sb,
                         sa * cc + ca * cb * sc,
                         -sa * sc + ca * cb * cc,
                         -ca * sb,
                         sb * sc,
                         sb * cc,
                         cb};
    std::copy(r, r + 9, rotations.begin() + 9 * i);
  }
  auto ok = true;
  for (const auto uniform : {false, true}) {
    // reference: in place rotation at each point
    auto ref = src;
    for (std::size_t i = 0; i != n; ++i) {
      const auto* const rv = rotations.data() + (uniform ? 0 : 9 * i);
      fg(ref.data() + i * ns, ref.data() + i * ns, rv);
    }
    for (const auto by_component : {false, true}) {
      const auto source = BatchField(src, n, ns, by_component);
      const auto destination =
          BatchField(std::vector<double>(n * ns), n, ns, by_component);
      const auto matrices = BatchField(
          uniform
              ? std::vector<double>(rotations.begin(), rotations.begin() + 9)
              : rotations,
          uniform ? 1 : n, 9, false);
      auto d = mfront_gpb_RotationData{};
      d.number_of_points = n;
      d.destination = destination.field();
      d.source = source.const_field();
      d.rotation_matrices = matrices.const_field();
      if (uniform) {
        d.rotation_matrices.point_stride = 0;
      }
      d.queue = queue();
      auto ms = 0.;
      const auto o =
          measure(ms, [&](mfront_gpb_OutputStatus& r) { fb(&r, &d); });
      auto batch = std::vector<double>(n * ns);
      destination.copy(batch);
      const auto status = std::vector<int>(n, 1);
      const auto diff = relativeDifference(ref, batch, status, ns);
      // the test is meaningless if the rotation changes nothing
      const auto rotated = ref != src;
      const auto pass = rotated && (o.status == 1) && (diff < tolerance);
      ok = pass && ok;
      std::printf(
          "rotation %-8s %-12s status=%2d | relative difference %.1e | %.2f "
          "ms | %s\n",
          uniform ? "uniform" : "by point",
          by_component ? "by component" : "by point", o.status, diff, ms,
          pass ? "OK" : "FAILED");
    }
  }
  return ok;
}  // end of testRotation

int main(const int argc, const char* const* const argv) {
  if (argc < 13) {
    fail(std::string{"usage: "} + argv[0] +
         " generic_library batch_library function gradients_size "
         "thermodynamic_forces_size material_properties_size "
         "internal_state_variables_size external_state_variables_size "
         "tangent_operator_size strain_amplitude material_properties "
         "temperature [--finite-strain=stress_measure,tangent_operator] "
         "[--initialize-function=number_of_inputs] "
         "[--post-processing=number_of_outputs] "
         "[--initial-internal-state-variables=values] "
         "[--speed-of-sound=mass_density] [--energies=true] "
         "[--tangent-operator=type] "
         "[--perturbation-of-internal-state-variables=amplitude] "
         "[--out-of-bounds-policy=behaviour,policy] "
         "[--rotation=number_of_values]");
  }
  [[maybe_unused]] const auto runtime = Runtime{};
  const auto s =
      Sizes{std::stoul(argv[4]), std::stoul(argv[5]), std::stoul(argv[6]),
            std::stoul(argv[7]), std::stoul(argv[8]), std::stoul(argv[9])};
  const auto amplitude = std::stod(argv[10]);
  const auto mp_values = parseList(argv[11]);
  const auto T = std::stod(argv[12]);
  if (mp_values.size() != s.mp) {
    fail("invalid number of material properties");
  }
  // not a multiple of the number of threads per block. The number of points
  // can be changed by the `MFRONT_GPB_COMPARE_NUMBER_OF_POINTS` environment
  // variable (for debugging).
  const auto* const np = std::getenv("MFRONT_GPB_COMPARE_NUMBER_OF_POINTS");
  const std::size_t n = np == nullptr ? 100003 : std::stoul(np);
  auto t = Test{.generic_library = argv[1],
                .batch_library = argv[2],
                .function = argv[3],
                .sizes = s,
                .options = parseOptions(argc, argv, 13),
                .initial = States(n, s),
                .identity = std::vector<double>(s.g, 0),
                .n = n,
                .dt = 1};
  // the identity, for the deformation gradients of finite strain behaviours
  if (t.options.finite_strain.finite_strain) {
    // the first three components of a tensor are its diagonal components
    for (std::size_t j = 0; j != std::min<std::size_t>(3, s.g); ++j) {
      t.identity[j] = 1;
    }
  }
  auto rng = std::mt19937{42};
  auto u = std::uniform_real_distribution<double>{-1, 1};
  auto& initial = t.initial;
  for (std::size_t i = 0; i != n; ++i) {
    for (std::size_t j = 0; j != s.g; ++j) {
      initial.g0[i * s.g + j] = t.identity[j];
      initial.g1[i * s.g + j] = t.identity[j] + amplitude * u(rng);
    }
    for (std::size_t j = 0; j != s.mp; ++j) {
      initial.mp[i * s.mp + j] = mp_values[j];
    }
    for (std::size_t j = 0; j != s.esv; ++j) {
      initial.esv0[i * s.esv + j] = initial.esv1[i * s.esv + j] = T;
    }
  }
  initial.mass_density = t.options.mass_density;
  const auto& isvs = t.options.initial_internal_state_variables;
  if ((isvs.size() != 1) && (isvs.size() != s.isv)) {
    fail("invalid number of initial values of the internal state variables");
  }
  const auto isv_perturbation = t.options.internal_state_variables_perturbation;
  for (std::size_t i = 0; i != n; ++i) {
    for (std::size_t j = 0; j != s.isv; ++j) {
      initial.isv0[i * s.isv + j] = initial.isv1[i * s.isv + j] =
          (isvs.size() == 1 ? isvs[0] : isvs[j]) + isv_perturbation * u(rng);
    }
  }
  // out of bounds policy, set in both libraries
  if (t.options.out_of_bounds_policy >= 0) {
    using SetOutOfBoundsPolicy = void (*)(const int);
    for (const auto& l : {t.generic_library, t.batch_library}) {
      const auto f = reinterpret_cast<SetOutOfBoundsPolicy>(
          load(l, t.options.behaviour + "_setOutOfBoundsPolicy"));
      f(t.options.out_of_bounds_policy);
    }
  }
  if (t.options.mode != Mode::INTEGRATION) {
    // initialize functions and post-processings are executed on arbitrary
    // states: stress-like thermodynamic forces and small positive internal
    // state variables
    for (auto* const v : {&initial.th0, &initial.th1}) {
      std::generate(v->begin(), v->end(), [&] { return 1e8 * u(rng); });
    }
    for (auto* const v : {&initial.isv0, &initial.isv1}) {
      std::generate(v->begin(), v->end(),
                    [&] { return 1e-3 * (1 + u(rng)) / 2; });
    }
  }
  const auto ok = [&t] {
    switch (t.options.mode) {
      case Mode::INITIALIZE_FUNCTION:
        return testInitializeFunction(t);
      case Mode::POST_PROCESSING:
        return testPostProcessing(t);
      case Mode::ROTATION:
        return testRotation(t);
      default:
        return testIntegration(t);
    }
  }();
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
