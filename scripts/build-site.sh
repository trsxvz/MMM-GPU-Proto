#!/bin/bash
# Builds the website of the prototype (docs/plan-ci-site.md, 9g):
# - the pages of site/ (Sphinx, MyST Markdown) at the root;
# - the web pages of TFEL (tfel/) and of MGIS (mgis/), built by their own
#   CMake target `website` (pandoc and pandoc-crossref);
# - the Sphinx documentation of mfem-mgis (mfem-mgis/) and its Doxygen
#   documentation (mfem-mgis/doxygen/).
# Requirements: the sources of MGIS with the patches of the prototype in
# mgis/ (only its series of patches is in this repository, see
# .github/workflows/pages.yml), and in the PATH: sphinx-build with the
# packages of site/requirements.txt, pandoc and pandoc-crossref (built with
# the same version of pandoc), doxygen, cmake, and C and C++ compilers (the
# configuration of TFEL and MGIS needs them, nothing is compiled).
# usage: scripts/build-site.sh [work directory]
# The site is written in <work directory>/html (build/site/html by default).
set -eu
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
W=$(mkdir -p "${1:-$P/build/site}" && cd "${1:-$P/build/site}" && pwd)
H=$W/html
# the work directory is either empty or a work directory of this script
if [ -n "$(ls -A "$W")" ] && [ ! -e "$W/.build-site" ]; then
  echo "build-site.sh: $W is not empty and not a work directory of" \
    "build-site.sh" >&2
  exit 1
fi
touch "$W/.build-site"
if [ ! -f "$P/mgis/CMakeLists.txt" ]; then
  echo "build-site.sh: sources of MGIS not found in $P/mgis" >&2
  exit 1
fi
for program in sphinx-build pandoc pandoc-crossref doxygen cmake; do
  command -v "$program" > /dev/null || {
    echo "build-site.sh: $program not found" >&2
    exit 1
  }
done
(cd "$W" && rm -rf html stage mfem-mgis-doxygen)
mkdir -p "$H" "$W/stage"

# runs a command, its outputs going to a log file of the work directory,
# whose end is printed if the command fails
# usage: run <log name> <command...>
run() {
  local log="$W/$1.log"
  shift
  "$@" > "$log" 2>&1 || {
    tail -n 30 "$log" >&2
    echo "build-site.sh: failure, see $log" >&2
    exit 1
  }
}

# web pages of TFEL and MGIS: the CMake target `website`, then the
# installation of its component (pages, style sheets, images)
website() {
  local name=$1
  shift
  run "$name-configure" cmake -S "$P/$name" -B "$W/$name" \
    -DCMAKE_BUILD_TYPE=Release -Denable-website=ON "$@"
  run "$name-build" cmake --build "$W/$name" --target website -j "$(nproc)"
  run "$name-install" cmake --install "$W/$name" --component website \
    --prefix "$W/stage/$name"
  cp -r "$W/stage/$name/share/doc/$name/web" "$H/$name"
}
website tfel -Denable-python=OFF -Denable-python-bindings=OFF \
  -Denable-fortran=OFF -Denable-testing=OFF
website mgis -Denable-testing=OFF

# mfem-mgis: Sphinx, then Doxygen (with the image of its style sheet,
# installed by the CMake of mfem-mgis)
run mfem-mgis-sphinx sphinx-build -b html -d "$W/doctrees/mfem-mgis" \
  "$P/mfem-mgis/docs/source" "$H/mfem-mgis"
run mfem-mgis-doxyfile cmake -DPROJECT_SOURCE_DIR="$P/mfem-mgis" \
  -DPROJECT_BINARY_DIR="$W/mfem-mgis-doxygen" \
  -P "$P/scripts/site-doxyfile.cmake"
(cd "$W/mfem-mgis-doxygen" && run mfem-mgis-doxygen doxygen Doxyfile)
cp -r "$W/mfem-mgis-doxygen/docs/doxygen/html" "$H/mfem-mgis/doxygen"
mkdir -p "$H/mfem-mgis/doxygen/img"
cp "$P/mfem-mgis/docs/doxygen/img/background.svg" "$H/mfem-mgis/doxygen/img/"

# pages of the prototype, at the root of the site (sphinx-build writes its
# own files only, the other directories are kept)
run site-sphinx sphinx-build -b html -d "$W/doctrees/site" "$P/site" "$H"
echo "website: $H/index.html"
