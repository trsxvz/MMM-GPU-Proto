#!/bin/bash
# Targeted test of some laws of the generic-parallel behaviour tests with
# the development tree: builds only the given laws (build-libraries.sh) for
# a configuration, then runs the comparisons of compare-all.sh matching the
# pattern. Reuses scratch/dev-test.sh, which builds all the laws, through
# DEVTEST_BEHAVIOURS. The pattern must only select comparisons of the given
# laws: a comparison of a law which is not built fails (for instance, the
# pattern ImplicitNorton also selects HenckyImplicitNorton). The other
# comparisons of the test (outside compare-all.sh) are skipped.
# usage: scratch/targeted-test.sh configuration pattern law...
c=$1
pattern=$2
shift 2
cd /home/t/Documents/GitHub/mfem-mgis-gpu-proto
DEVTEST_BEHAVIOURS="$*" DEVTEST_OPTIMIZATION_LEVEL=-O3 MAKEFLAGS=-j8 \
  scratch/dev-test.sh "$c" "$pattern"
