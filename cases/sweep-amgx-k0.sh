#!/bin/bash
# Sweep of configurations of AmgX on the elastic operator K0 of the beam r4
# (tool beam-k0, GPU): conjugate gradient of MFEM preconditioned by one
# cycle of AmgX described by each JSON file of the configurations
# directory, relative tolerance 1e-8, five resolutions.
# The configurations are copied with the raw outputs; summary
# (setup, iterations, fastest and slowest resolution after the first one,
# resolutions not converged) in summary.txt.
# To be run on a quiet machine, after build-mfem-mgis.sh cases.
#
# usage: cases/sweep-amgx-k0.sh configurations-directory [results directory]
#        cases/sweep-amgx-k0.sh --resume results-directory
# The second form computes summary.txt again from the raw outputs.
set -u
P=$(cd "$(dirname "$0")/.." && pwd)

# summary of the output of beam-k0 for one configuration: setup time,
# iterations and times of the resolutions after the first, resolutions not
# converged (decimal point of the C locale)
# usage: summarise name output-file
summarise() {
  LC_ALL=C awk -v name="$1" -F'|' '
    /^dofs/ { split($3, s, " "); setup = s[2] }
    /^resolution/ {
      split($2, t, " "); split($3, i, " ")
      if (n == 0 || t[2] < tmin) { tmin = t[2] }
      if (n == 0 || t[2] > tmax) { tmax = t[2] }
      if (n == 0 || i[2] < imin) { imin = i[2] }
      if (n == 0 || i[2] > imax) { imax = i[2] }
      if ($0 ~ /not converged/) { ++failures }
      ++n
    }
    END {
      if (n == 0) { printf "%s: failure, see %s.txt\n", name, name; exit }
      iterations = (imin == imax) ? imin : imin " to " imax
      printf "%s: setup %s s, %s iterations, resolution %.4f to %.4f s",
             name, setup, iterations, tmin, tmax
      if (failures > 0) {
        printf ", %d resolutions out of %d not converged", failures, n
      }
      printf "\n"
    }' "$2"
}

if [ "${1:-}" = --resume ]; then
  out=$(cd "${2:?results directory expected}" && pwd) || exit 1
  : > "$out/summary.txt"
  for c in "$out"/configurations/*.json; do
    name=$(basename "$c" .json)
    summarise "$name" "$out/$name.txt" | tee -a "$out/summary.txt"
  done
  exit 0
fi

configurations=$(cd "${1:?configurations directory expected}" && pwd) ||
  exit 1
if ! compgen -G "$configurations/*.json" > /dev/null; then
  echo "sweep-amgx-k0.sh: no configuration in $configurations" >&2
  exit 1
fi
source "$P/scripts/env.sh" > /dev/null
unset OMP_PROC_BIND OMP_PLACES
out=${2:-$P/results/sweep-amgx-k0-$(date +%Y-%m-%d-%H%M)}
mkdir -p "$out/configurations"
out=$(cd "$out" && pwd)
: > "$out/summary.txt"
B=$P/build/cases/beam
for c in "$configurations"/*.json; do
  name=$(basename "$c" .json)
  cp "$c" "$out/configurations/"
  (cd "$B" && env OMP_NUM_THREADS=1 timeout 600 ./beam-k0 -r 4 -n 5 \
    --device cuda -s cg-amgx -amgxc "$c") < /dev/null > "$out/$name.txt" 2>&1
  summarise "$name" "$out/$name.txt" | tee -a "$out/summary.txt"
done
echo "results in $out"
