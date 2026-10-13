#!/bin/bash
# Checks that a series of patches of `patches/<directory>/` applies on its
# base commit and gives exactly the tree of the branch of the prototype:
# the patches are applied with `git am --keep-cr` (some files of TFEL have
# CRLF line endings) in a temporary worktree of the local repository, then
# the trees are compared.
#
# usage: scripts/verify-patches.sh repository patches-directory base branch
# example: scripts/verify-patches.sh mgis mgis e8ef4b4 gpu-proto
set -u
PROTO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
repository=${1:?usage: verify-patches.sh repository patches base branch}
patches=${2:?usage: verify-patches.sh repository patches base branch}
base=${3:?usage: verify-patches.sh repository patches base branch}
branch=${4:?usage: verify-patches.sh repository patches base branch}
cd "$PROTO" || exit 1
worktree=$PROTO/scratch/verify-patches-$patches
git -C "$repository" worktree add -q --detach "$worktree" "$base" || exit 1
status=0
if ! (cd "$worktree" &&
  git -c user.name=verify -c user.email=verify@localhost am -q --keep-cr \
    "$PROTO"/patches/"$patches"/*.patch); then
  echo "patches/$patches: the patches do not apply on $base" >&2
  status=1
else
  applied=$(git -C "$worktree" rev-parse "HEAD^{tree}")
  expected=$(git -C "$repository" rev-parse "$branch^{tree}")
  if [ "$applied" = "$expected" ]; then
    echo "patches/$patches: same tree as $branch ($applied)"
  else
    echo "patches/$patches: tree $applied, $branch has $expected" >&2
    status=1
  fi
fi
git -C "$repository" worktree remove --force "$worktree"
git -C "$repository" worktree prune
exit $status
