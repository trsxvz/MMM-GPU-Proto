#!/bin/bash
# Copy the last commit of the local TFEL repository (`tfel/`) into this
# repository, without the history of TFEL:
# - the series of patches since the base commit is regenerated in
#   `patches/tfel/`;
# - the sources of this commit replace the files of `tfel/` tracked by this
#   repository (the `tfel/` directory is both the local TFEL repository and
#   the copy tracked here).
# Only the index is modified: the changes are then committed as usual.
# usage: scripts/sync-tfel.sh
set -eu
PROTO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BASE=ff0959ece
cd "$PROTO"
if [ -n "$(git -C tfel status --porcelain --untracked-files=no)" ]; then
  echo "tfel/ has uncommitted changes" >&2
  exit 1
fi
# patches
rm -f patches/tfel/*.patch
git -C tfel format-patch -q --no-numbered "$BASE..HEAD" -o "$PROTO/patches/tfel/"
git add -A patches/tfel
# sources: only the objects which are not already known are copied
tree=$(git -C tfel rev-parse HEAD^{tree})
previous=$(git rev-parse -q --verify HEAD:tfel || true)
if [ -n "$previous" ] && git -C tfel cat-file -e "$previous" 2> /dev/null; then
  exclude=("--not" "$previous")
else
  exclude=()
fi
git -C tfel rev-list --objects "$tree" "${exclude[@]}" | cut -d' ' -f1 |
  git -C tfel pack-objects -q --stdout | git index-pack --stdin > /dev/null
git rm -r -q --cached --ignore-unmatch tfel
git read-tree --prefix=tfel/ "$tree"
echo "tfel/ synchronized with $(git -C tfel log -1 --format='%h %s')"
