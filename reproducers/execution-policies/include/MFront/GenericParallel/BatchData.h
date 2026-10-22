/*!
 * \file   MFront/GenericParallel/BatchData.h
 * \brief  Data given to a generated function, as `mfront_gpmp_BatchData`:
 * the policy is an opaque pointer in C. The structure
 * `mfront_gp_ReproducerReport` exists in this reproducer only: the
 * generated functions write in it what they have read from the policy.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_BATCHDATA_H
#define LIB_MFRONT_GENERICPARALLEL_BATCHDATA_H

#ifdef __cplusplus

#include "MFront/GenericParallel/AbstractBackendExecutionPolicy.hxx"

extern "C" {

typedef mfront::generic_parallel::AbstractBackendExecutionPolicy
    mfront_gp_AbstractBackendExecutionPolicy;

#else /* __cplusplus */

/*! \brief opaque structure, created by C++ code only */
typedef struct mfront_gp_AbstractBackendExecutionPolicy
    mfront_gp_AbstractBackendExecutionPolicy;

#endif /* __cplusplus */

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#define MFRONT_GENERICPARALLEL_GENERATED_VISIBILITY __declspec(dllexport)
#else
#define MFRONT_GENERICPARALLEL_GENERATED_VISIBILITY \
  __attribute__((visibility("default")))
#endif

/*! \brief statistics of an asynchronous call */
typedef struct mfront_gp_Statistics {
  int number_of_computed_points;
} mfront_gp_Statistics;

/*! \brief what a generated function has read from the policy */
typedef struct mfront_gp_ReproducerReport {
  /*! \brief the default policy of the backend was used (null pointer) */
  int default_policy;
  int preserve_stack_size;
  int asynchronous;
  /*! \brief stream, device of the queue or instance of the execution space */
  long long queue;
  /*! \brief number of owners of the queue seen by the generated function */
  long queue_owners;
} mfront_gp_ReproducerReport;

typedef struct mfront_gp_BatchData {
  int number_of_points;
  int* statuses;
  mfront_gp_Statistics* statistics;
  /*! \brief policy, null for the default policy of the backend */
  const mfront_gp_AbstractBackendExecutionPolicy* execution_policy;
  mfront_gp_ReproducerReport* report;
} mfront_gp_BatchData;

/*! \brief signature of the generated functions */
typedef int (*mfront_gp_GeneratedFunction)(const mfront_gp_BatchData*);

#ifdef __cplusplus
} /* end of extern "C" */
#endif /* __cplusplus */

#endif /* LIB_MFRONT_GENERICPARALLEL_BATCHDATA_H */
