/*!
 * \file   src/main-linked.cxx
 * \brief  Program linked to the library of the caller.
 * usage: execution-policies-linked generated [dynamic-cast|-] [c-types]
 */

#include <cstdio>
#include <cstring>
#include "Caller.h"

int main(const int argc, const char* const* const argv) {
  if ((argc < 2) || (argc > 4)) {
    std::fprintf(stderr, "usage: %s generated [dynamic-cast|-] [c-types]\n",
                 argv[0]);
    return 1;
  }
  const auto* const dynamic_cast_library =
      ((argc > 2) && (std::strcmp(argv[2], "-") != 0)) ? argv[2] : nullptr;
  const auto c_types_only =
      (argc > 3) && (std::strcmp(argv[3], "c-types") == 0);
  return reproducer_run(argv[1], dynamic_cast_library, c_types_only) == 0 ? 0
                                                                          : 1;
}
