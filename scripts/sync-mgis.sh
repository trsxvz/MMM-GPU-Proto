#!/bin/bash
# Regenerate in `patches/mgis/` the series of patches of the local MGIS
# repository (`mgis/`, branch `gpu-proto`) since the base commit (MGIS master
# used by the prototype). Only the patches are tracked by this repository.
# Only the index is modified: the changes are then committed as usual.
# usage: scripts/sync-mgis.sh
set -eu
PROTO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BASE=e8ef4b4
cd "$PROTO"
if [ -n "$(git -C mgis status --porcelain --untracked-files=no)" ]; then
  echo "mgis/ has uncommitted changes" >&2
  exit 1
fi
rm -f patches/mgis/*.patch
git -C mgis format-patch -q --no-numbered "$BASE..gpu-proto" -o "$PROTO/patches/mgis/"
git add -A patches/mgis
echo "patches/mgis synchronized with $(git -C mgis log -1 --format='%h %s' gpu-proto)"
