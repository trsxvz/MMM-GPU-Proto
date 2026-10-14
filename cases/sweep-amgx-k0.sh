#!/bin/bash
# Balayage de configurations d'AmgX sur l'opérateur élastique K0 de la
# poutre r4 (outil poutre-k0, GPU) : gradient conjugué de MFEM
# préconditionné par un cycle d'AmgX décrit par chaque fichier JSON du
# dossier de configurations, tolérance relative 1e-8, cinq résolutions.
# Les configurations sont copiées avec les sorties brutes ; résumé
# (mise en place, itérations, résolution la plus rapide et la plus lente
# hors la première, résolutions non convergées) dans summary.txt.
# À lancer sur une machine calme, après build-mfem-mgis.sh cas.
#
# usage : cases/sweep-amgx-k0.sh dossier-de-configurations [dossier de résultats]
#         cases/sweep-amgx-k0.sh --resume dossier-de-résultats
# La seconde forme recalcule summary.txt à partir des sorties brutes.
set -u
P=$(cd "$(dirname "$0")/.." && pwd)

# summary of the output of poutre-k0 for one configuration: setup time,
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
      if (n == 0) { printf "%s : échec, voir %s.txt\n", name, name; exit }
      iterations = (imin == imax) ? imin : imin " à " imax
      printf "%s : mise en place %s s, %s itérations, résolution %.4f à %.4f s",
             name, setup, iterations, tmin, tmax
      if (failures > 0) {
        printf ", %d résolutions sur %d non convergées", failures, n
      }
      printf "\n"
    }' "$2"
}

if [ "${1:-}" = --resume ]; then
  out=$(cd "${2:?dossier de résultats attendu}" && pwd) || exit 1
  : > "$out/summary.txt"
  for c in "$out"/configurations/*.json; do
    name=$(basename "$c" .json)
    summarise "$name" "$out/$name.txt" | tee -a "$out/summary.txt"
  done
  exit 0
fi

configurations=$(cd "${1:?dossier de configurations attendu}" && pwd) ||
  exit 1
if ! compgen -G "$configurations/*.json" > /dev/null; then
  echo "sweep-amgx-k0.sh: aucune configuration dans $configurations" >&2
  exit 1
fi
source "$P/scripts/env.sh" > /dev/null
unset OMP_PROC_BIND OMP_PLACES
out=${2:-$P/results/balayage-amgx-k0-$(date +%Y-%m-%d-%H%M)}
mkdir -p "$out/configurations"
out=$(cd "$out" && pwd)
: > "$out/summary.txt"
B=$P/build/cas/poutre
for c in "$configurations"/*.json; do
  name=$(basename "$c" .json)
  cp "$c" "$out/configurations/"
  (cd "$B" && env OMP_NUM_THREADS=1 timeout 600 ./poutre-k0 -r 4 -n 5 \
    --device cuda -s cg-amgx -amgxc "$c") < /dev/null > "$out/$name.txt" 2>&1
  summarise "$name" "$out/$name.txt" | tee -a "$out/summary.txt"
done
echo "résultats dans $out"
