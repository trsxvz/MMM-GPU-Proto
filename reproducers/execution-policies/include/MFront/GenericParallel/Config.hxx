/*!
 * \file   MFront/GenericParallel/Config.hxx
 * \brief  Visibility of the symbols of the runtime library, as the macros of
 * `TFEL/Config/TFELConfig.hxx`: `dllexport` when building the library and
 * `dllimport` otherwise on Windows, default visibility elsewhere (all the
 * other symbols are hidden, as with the flags given by `tfel-config`).
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_CONFIG_HXX
#define LIB_MFRONT_GENERICPARALLEL_CONFIG_HXX

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#if defined MFrontGenericParallelRuntime_EXPORTS
#define MFRONT_GENERICPARALLEL_RUNTIME_VISIBILITY __declspec(dllexport)
#else
#define MFRONT_GENERICPARALLEL_RUNTIME_VISIBILITY __declspec(dllimport)
#endif
#else
#define MFRONT_GENERICPARALLEL_RUNTIME_VISIBILITY \
  __attribute__((visibility("default")))
#endif

#endif /* LIB_MFRONT_GENERICPARALLEL_CONFIG_HXX */
