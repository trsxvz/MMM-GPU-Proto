/*!
 * \file   src/main-c.c
 * \brief  C caller: the header of the data is valid C, and a null policy
 * selects the default policy of the backend (except for `sycl`).
 * usage: execution-policies-c generated
 */

#include <stdio.h>
#include "Loader.h"
#include "MFront/GenericParallel/BatchData.h"

static int failures = 0;

static void check(const int b, const char* const what) {
  printf("  %-66s %s\n", what, b ? "ok" : "FAILED");
  if (!b) {
    ++failures;
  }
}

static int call(reproducer_library l, const char* const n) {
  int statuses[4] = {0, 0, 0, 0};
  mfront_gp_BatchData d;
  mfront_gp_GeneratedFunction f =
      (mfront_gp_GeneratedFunction)reproducer_symbol(l, n);
  if (f == NULL) {
    fprintf(stderr, "no function '%s'\n", n);
    return -100;
  }
  d.number_of_points = 4;
  d.statuses = statuses;
  d.statistics = NULL;
  d.execution_policy = NULL;
  d.report = NULL;
  return f(&d);
}

int main(const int argc, const char* const* const argv) {
  reproducer_library l;
  if (argc != 2) {
    fprintf(stderr, "usage: %s generated\n", argv[0]);
    return 1;
  }
  l = reproducer_load(argv[1]);
  if (l == NULL) {
    reproducer_print_error(argv[1]);
    return 1;
  }
  printf("C caller, null policies\n");
  check(call(l, "reproducer_cuda") == 0, "cuda: default policy");
  check(call(l, "reproducer_stlpar") == 0, "stlpar: default policy");
  check(call(l, "reproducer_sycl") == -5,
        "sycl: refused, the queue is required (-5)");
  printf("%s\n",
         failures == 0 ? "all the checks succeeded" : "SOME CHECKS FAILED");
  return failures == 0 ? 0 : 1;
}
