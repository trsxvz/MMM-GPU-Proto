#!/bin/bash
# Copy the last commit of a local repository (`mfem-mgis/`,
# `mfem-mgis-examples/` or `mm-opera-hpc/`) into this repository, without
# its history, as scripts/sync-tfel.sh does for TFEL:
# - the series of patches since the base commit is regenerated in
#   `patches/<repository>/`;
# - the sources of this commit replace the files of `<repository>/` tracked
#   by this repository (the directory is both the local repository and the
#   copy tracked here).
# Only the index is modified: the changes are then committed as usual.
# usage: scripts/sync-repository.sh repository
set -eu
PROTO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
repository=${1:?usage: scripts/sync-repository.sh repository}
# base commits: master of mfem-mgis (official repository) and of
# mfem-mgis-examples (latug0/mfem-mgis-examples, the reference repository of
# the examples, consistent with this commit of mfem-mgis), main of
# mm-opera-hpc (rprat-pro/mm-opera-hpc, its reference repository)
case "$repository" in
  mfem-mgis) BASE=5bbd5a6981 ;;
  mfem-mgis-examples) BASE=f256bd8 ;;
  mm-opera-hpc) BASE=f3fc124 ;;
  *)
    echo "unknown repository '$repository'" >&2
    exit 1
    ;;
esac
cd "$PROTO"
if [ -n "$(git -C "$repository" status --porcelain --untracked-files=no)" ]; then
  echo "$repository/ has uncommitted changes" >&2
  exit 1
fi
# patches
mkdir -p "patches/$repository"
rm -f "patches/$repository"/*.patch
git -C "$repository" format-patch -q --no-numbered "$BASE..HEAD" \
  -o "$PROTO/patches/$repository/"
git add -A "patches/$repository"
# sources: only the objects which are not already known are copied
tree=$(git -C "$repository" rev-parse "HEAD^{tree}")
previous=$(git rev-parse -q --verify "HEAD:$repository" || true)
if [ -n "$previous" ] &&
  git -C "$repository" cat-file -e "$previous" 2> /dev/null; then
  exclude=("--not" "$previous")
else
  exclude=()
fi
git -C "$repository" rev-list --objects "$tree" "${exclude[@]}" |
  cut -d' ' -f1 |
  git -C "$repository" pack-objects -q --stdout |
  git index-pack --stdin > /dev/null
# the index entries are replaced, even if a previous synchronization was
# not committed yet
git rm -r -q --cached --ignore-unmatch -f "$repository"
git read-tree --prefix="$repository/" "$tree"
echo "$repository/ synchronized with" \
  "$(git -C "$repository" log -1 --format='%h %s')"
