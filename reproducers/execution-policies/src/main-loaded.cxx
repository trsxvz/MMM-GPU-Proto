/*!
 * \file   src/main-loaded.cxx
 * \brief  Program loading the library of the caller at run time, as Python
 * loads a module.
 * usage: execution-policies-loaded caller generated [dynamic-cast]
 */

#include <cstdio>
#include "Loader.h"
#include "Caller.h"

int main(const int argc, const char* const* const argv) {
  if ((argc != 3) && (argc != 4)) {
    std::fprintf(stderr, "usage: %s caller generated [dynamic-cast]\n",
                 argv[0]);
    return 1;
  }
  const auto l = reproducer_load(argv[1]);
  if (l == nullptr) {
    reproducer_print_error(argv[1]);
    return 1;
  }
  const auto run = reinterpret_cast<reproducer_run_function>(
      reproducer_symbol(l, "reproducer_run"));
  if (run == nullptr) {
    std::fprintf(stderr, "no function 'reproducer_run'\n");
    return 1;
  }
  return run(argv[2], argc == 4 ? argv[3] : nullptr, 0) == 0 ? 0 : 1;
}
