/*!
 * \file   src/Loader.h
 * \brief  Loading of a library at run time: `dlopen` with `RTLD_NOW`, as
 * MGIS, or `LoadLibrary` on Windows.
 */

#ifndef LIB_LOADER_H
#define LIB_LOADER_H

#include <stdio.h>

#if defined _WIN32 || defined _WIN64
#include <windows.h>
typedef HMODULE reproducer_library;
static inline reproducer_library reproducer_load(const char* const path) {
  return LoadLibraryA(path);
}
static inline void* reproducer_symbol(reproducer_library l,
                                      const char* const n) {
  return (void*)GetProcAddress(l, n);
}
static inline void reproducer_print_error(const char* const path) {
  fprintf(stderr, "can't load '%s' (error %lu)\n", path,
          (unsigned long)GetLastError());
}
#else
#include <dlfcn.h>
typedef void* reproducer_library;
static inline reproducer_library reproducer_load(const char* const path) {
  return dlopen(path, RTLD_NOW);
}
static inline void* reproducer_symbol(reproducer_library l,
                                      const char* const n) {
  return dlsym(l, n);
}
static inline void reproducer_print_error(const char* const path) {
  fprintf(stderr, "can't load '%s' (%s)\n", path, dlerror());
}
#endif

#endif /* LIB_LOADER_H */
