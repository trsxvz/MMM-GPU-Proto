/*!
 * \file   src/Caller.h
 * \brief  Stand-in for MGIS or a user code: a library linked to the runtime
 * library, which creates the policies and calls the generated functions.
 */

#ifndef LIB_CALLER_H
#define LIB_CALLER_H

#if defined _WIN32 || defined _WIN64
#if defined caller_EXPORTS
#define REPRODUCER_CALLER_VISIBILITY __declspec(dllexport)
#else
#define REPRODUCER_CALLER_VISIBILITY __declspec(dllimport)
#endif
#else
#define REPRODUCER_CALLER_VISIBILITY __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \return the number of failed checks
 * \param[in] generated: path to the generated library
 * \param[in] dynamic_cast_library: path to the library checking the policies
 * by `dynamic_cast`, for information, or a null pointer
 * \param[in] c_types_only: only check the policies holding `C` types
 * (`cuda`, `stlpar`), for a generated library built by another toolchain
 */
REPRODUCER_CALLER_VISIBILITY int reproducer_run(
    const char* generated, const char* dynamic_cast_library, int c_types_only);

/*! \brief signature of `reproducer_run` */
typedef int (*reproducer_run_function)(const char*, const char*, int);

#ifdef __cplusplus
}
#endif

#endif /* LIB_CALLER_H */
