/*!
* \file   SaintVenantKirchhoffElasticity-generic-parallel-cuda.cu
* \brief  This file implements the generic-parallel interface for the SaintVenantKirchhoffElasticity behaviour (backend cuda)
* \date   19 / 10 / 2013
*/

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif /* NOMINMAX */
#include <windows.h>
#ifdef small
#undef small
#endif /* small */
#endif /* _WIN32 */

#ifndef MFRONT_SHAREDOBJ
#define MFRONT_SHAREDOBJ TFEL_VISIBILITY_EXPORT
#endif /* MFRONT_SHAREDOBJ */

#ifndef MFRONT_EXPORT_SYMBOL
#define MFRONT_EXPORT_SYMBOL(TYPE, NAME, VALUE) \
  MFRONT_SHAREDOBJ extern TYPE NAME;            \
  MFRONT_SHAREDOBJ TYPE NAME = VALUE
#endif /* MFRONT_EXPORT_SYMBOL*/

#ifndef MFRONT_EXPORT_ARRAY_ARGUMENTS
#define MFRONT_EXPORT_ARRAY_ARGUMENTS(...) __VA_ARGS__
#endif /* MFRONT_EXPORT_ARRAY_ARGUMENTS */

#ifndef MFRONT_EXPORT_ARRAY_OF_SYMBOLS
#define MFRONT_EXPORT_ARRAY_OF_SYMBOLS(TYPE, NAME, SIZE, VALUE) \
  MFRONT_SHAREDOBJ extern TYPE NAME[SIZE];                      \
  MFRONT_SHAREDOBJ TYPE NAME[SIZE] = {VALUE}
#endif /* MFRONT_EXPORT_ARRAY_OF_SYMBOLS*/

#include"TFEL/Material/OutOfBoundsPolicy.hxx"
#include"TFEL/Material/SaintVenantKirchhoffElasticity.hxx"
#include"MFront/GenericBehaviour/GenericBehaviourTraits.hxx"
#include"MFront/GenericParallelBehaviour/Integrate.hxx"
#include"MFront/GenericParallelBehaviour/FiniteStrainIntegrators.hxx"
#include "MFront/GenericParallelBehaviour/GPU.hxx"
#include"MFront/GenericParallelBehaviour/SaintVenantKirchhoffElasticity-generic-parallel-cuda.hxx"

static tfel::material::OutOfBoundsPolicy&
SaintVenantKirchhoffElasticity_getOutOfBoundsPolicy(){
static auto policy = tfel::material::None;
return policy;
}

namespace mfront::gb{

template<>
struct GenericBehaviourTraits<tfel::material::SaintVenantKirchhoffElasticity<tfel::material::ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

template<>
struct GenericBehaviourTraits<tfel::material::SaintVenantKirchhoffElasticity<tfel::material::ModellingHypothesis::AXISYMMETRICAL, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::AXISYMMETRICAL;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

template<>
struct GenericBehaviourTraits<tfel::material::SaintVenantKirchhoffElasticity<tfel::material::ModellingHypothesis::PLANESTRAIN, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::PLANESTRAIN;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

template<>
struct GenericBehaviourTraits<tfel::material::SaintVenantKirchhoffElasticity<tfel::material::ModellingHypothesis::GENERALISEDPLANESTRAIN, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::GENERALISEDPLANESTRAIN;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

template<>
struct GenericBehaviourTraits<tfel::material::SaintVenantKirchhoffElasticity<tfel::material::ModellingHypothesis::TRIDIMENSIONAL, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::TRIDIMENSIONAL;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

} // end of namespace mfront::gb

namespace {

//! \brief function object calling the initialize function 'StressFromInitialPressure'
struct SaintVenantKirchhoffElasticityInitializeFunction_StressFromInitialPressure {
template <typename Behaviour>
TFEL_HOST_DEVICE void operator()(Behaviour& b, const mfront::gb::real* const values) const {
b.executeStressFromInitialPressureInitializeFunction(values);
}
};

//! \brief function object calling the post-processing 'EquivalentStress'
struct SaintVenantKirchhoffElasticityPostProcessing_EquivalentStress {
template <typename Behaviour>
TFEL_HOST_DEVICE void operator()(Behaviour& b, mfront::gb::real* const outputs, const typename Behaviour::BehaviourData& initial_state) const {
b.executeEquivalentStressPostProcessing(outputs, initial_state);
}
};

//! \brief function object calling the post-processing 'PrincipalGreenLagrangeStrain'
struct SaintVenantKirchhoffElasticityPostProcessing_PrincipalGreenLagrangeStrain {
template <typename Behaviour>
TFEL_HOST_DEVICE void operator()(Behaviour& b, mfront::gb::real* const outputs, const typename Behaviour::BehaviourData& initial_state) const {
b.executePrincipalGreenLagrangeStrainPostProcessing(outputs, initial_state);
}
};

//! \brief description of the data of an integration point
struct SaintVenantKirchhoffElasticityAxisymmetricalGeneralisedPlaneStrainPointTraits {
//! \brief maximal sizes of the data of an integration point
static constexpr auto value = mfront::gpb::PointDataSizes{.gradients = 3, .thermodynamic_forces = 3, .material_properties = 2, .internal_state_variables = 0, .external_state_variables = 1, .tangent_operator = 9};
//! \brief integration strategy
using Integrator = mfront::gpb::FiniteStrainIntegrator;
};

//! \brief values of the parameters
struct SaintVenantKirchhoffElasticityAxisymmetricalGeneralisedPlaneStrainParameters{
double minimal_time_step_scaling_factor;
double maximal_time_step_scaling_factor;
};

//! \brief description of the data of an integration point
struct SaintVenantKirchhoffElasticityAxisymmetricalPointTraits {
//! \brief maximal sizes of the data of an integration point
static constexpr auto value = mfront::gpb::PointDataSizes{.gradients = 5, .thermodynamic_forces = 5, .material_properties = 2, .internal_state_variables = 0, .external_state_variables = 1, .tangent_operator = 25};
//! \brief integration strategy
using Integrator = mfront::gpb::FiniteStrainIntegrator;
};

//! \brief values of the parameters
struct SaintVenantKirchhoffElasticityAxisymmetricalParameters{
double minimal_time_step_scaling_factor;
double maximal_time_step_scaling_factor;
};

//! \brief description of the data of an integration point
struct SaintVenantKirchhoffElasticityPlaneStrainPointTraits {
//! \brief maximal sizes of the data of an integration point
static constexpr auto value = mfront::gpb::PointDataSizes{.gradients = 5, .thermodynamic_forces = 5, .material_properties = 2, .internal_state_variables = 0, .external_state_variables = 1, .tangent_operator = 25};
//! \brief integration strategy
using Integrator = mfront::gpb::FiniteStrainIntegrator;
};

//! \brief values of the parameters
struct SaintVenantKirchhoffElasticityPlaneStrainParameters{
double minimal_time_step_scaling_factor;
double maximal_time_step_scaling_factor;
};

//! \brief description of the data of an integration point
struct SaintVenantKirchhoffElasticityGeneralisedPlaneStrainPointTraits {
//! \brief maximal sizes of the data of an integration point
static constexpr auto value = mfront::gpb::PointDataSizes{.gradients = 5, .thermodynamic_forces = 5, .material_properties = 2, .internal_state_variables = 0, .external_state_variables = 1, .tangent_operator = 25};
//! \brief integration strategy
using Integrator = mfront::gpb::FiniteStrainIntegrator;
};

//! \brief values of the parameters
struct SaintVenantKirchhoffElasticityGeneralisedPlaneStrainParameters{
double minimal_time_step_scaling_factor;
double maximal_time_step_scaling_factor;
};

//! \brief description of the data of an integration point
struct SaintVenantKirchhoffElasticityTridimensionalPointTraits {
//! \brief maximal sizes of the data of an integration point
static constexpr auto value = mfront::gpb::PointDataSizes{.gradients = 9, .thermodynamic_forces = 9, .material_properties = 2, .internal_state_variables = 0, .external_state_variables = 1, .tangent_operator = 81};
//! \brief integration strategy
using Integrator = mfront::gpb::FiniteStrainIntegrator;
};

//! \brief values of the parameters
struct SaintVenantKirchhoffElasticityTridimensionalParameters{
double minimal_time_step_scaling_factor;
double maximal_time_step_scaling_factor;
};

} // end of anonymous namespace

#include <cstdio>
using Behaviour = tfel::material::SaintVenantKirchhoffElasticity<
    tfel::material::ModellingHypothesis::TRIDIMENSIONAL, double, false>;
using Parameters = SaintVenantKirchhoffElasticityTridimensionalParameters;
using Functor =
    SaintVenantKirchhoffElasticityInitializeFunction_StressFromInitialPressure;
using PointTraits = SaintVenantKirchhoffElasticityTridimensionalPointTraits;

template <bool direct>
__global__ void kernel(const mfront_gpb_BatchData d,
                       const mfront_gpb_ConstField inputs,
                       int* const status) {
  using Executor = mfront::gpb::InitializeFunctionExecutor<
      Behaviour, PointTraits, Functor, 1, direct, Parameters>;
  const auto e = Executor{d, inputs, Parameters{0.1, 1e300},
                          tfel::material::None};
  const auto r = e(threadIdx.x);
  status[threadIdx.x] = static_cast<int>(r.number_of_failures);
}

template <typename Executor>
__global__ void kernel2(const Executor e, int* const status) {
#ifdef LOCAL_COPY
  const auto local = e;
  const auto r = local(threadIdx.x);
#else
  const auto r = e(threadIdx.x);
#endif
  status[threadIdx.x] = static_cast<int>(r.number_of_failures);
}

// fields of two points, stored component by component (or point by point)
template <bool direct>
void run() {
  constexpr int n = 2;
  double *g, *th0, *th1, *mp, *esv, *in;
  int* status;
  cudaMallocManaged(&g, 9 * n * sizeof(double));
  cudaMallocManaged(&th0, 6 * n * sizeof(double));
  cudaMallocManaged(&th1, 6 * n * sizeof(double));
  cudaMallocManaged(&mp, 2 * n * sizeof(double));
  cudaMallocManaged(&esv, n * sizeof(double));
  cudaMallocManaged(&in, n * sizeof(double));
  cudaMallocManaged(&status, n * sizeof(int));
  // strides
  const auto ps = [](int c) { return direct ? c : 1; };
  const auto cs = [](int) { return direct ? 1 : n; };
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != 9; ++j) g[i * ps(9) + j * cs(9)] = j < 3 ? 1 : 0;
    for (int j = 0; j != 6; ++j) th0[i * ps(6) + j * cs(6)] = -(j + 1);
    for (int j = 0; j != 6; ++j) th1[i * ps(6) + j * cs(6)] = j + 1;
    mp[i * ps(2)] = 150e9;
    mp[i * ps(2) + cs(2)] = 0.3;
    esv[i] = 293.15;
    in[i] = 42 + i;
  }
  auto d = mfront_gpb_BatchData{};
  d.number_of_points = n;
  d.dt = 1;
  d.rdt = 1;
  d.s0.gradients = {g, (size_t)ps(9), (size_t)cs(9)};
  d.s1.gradients = {g, (size_t)ps(9), (size_t)cs(9)};
  d.s0.thermodynamic_forces = {th0, (size_t)ps(6), (size_t)cs(6)};
  d.s1.thermodynamic_forces = {th1, (size_t)ps(6), (size_t)cs(6)};
  d.s0.material_properties = {mp, (size_t)ps(2), (size_t)cs(2)};
  d.s1.material_properties = {mp, (size_t)ps(2), (size_t)cs(2)};
  d.s0.external_state_variables = {esv, 1, 1};
  d.s1.external_state_variables = {esv, 1, 1};
#ifdef USE_BACKEND
  d.statuses = status;
  auto os = mfront_gpb_OutputStatus{};
  const auto r = mfront::gpb::gpu::executeInitializeFunction<
      Behaviour, PointTraits, Functor, 1>(os, d, mfront_gpb_ConstField{in, 1, 1},
                                          Parameters{0.1, 1e300},
                                          tfel::material::None, 64);
  std::printf("backend status %d\n", r);
#elif defined(USE_PARAMETER)
  using Executor = mfront::gpb::InitializeFunctionExecutor<
      Behaviour, PointTraits, Functor, 1, direct, Parameters>;
  kernel2<<<1, n>>>(Executor{d, mfront_gpb_ConstField{in, 1, 1},
                             Parameters{0.1, 1e300}, tfel::material::None},
                    status);
  cudaDeviceSynchronize();
#else
  kernel<direct><<<1, n>>>(d, mfront_gpb_ConstField{in, 1, 1}, status);
  cudaDeviceSynchronize();
#endif
  for (int i = 0; i != n; ++i) {
    std::printf("direct=%d point %d failures %d th1", direct, i, status[i]);
    for (int j = 0; j != 6; ++j) std::printf(" %g", th1[i * ps(6) + j * cs(6)]);
    std::printf(" (expected %d %d %d 0 0 0)\n", 42 + i, 42 + i, 42 + i);
  }
}

int main() {
  run<true>();
  run<false>();
}
