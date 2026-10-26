#!/bin/bash
# Side of the target machine, without access to the internet, of the
# installation of the prototype (docs/install.md). The deployment directory
# ROOT holds spack/, spack-packages/, bootstrap/ (mirror of the bootstrap
# of spack), the mirrors of the sources and the prototype (MMM-GPU-Proto/).
# ROOT is the directory of this script when it is the copy extracted from
# mmm-base.tar, otherwise the parent directory of the prototype; PROTO_DEPLOY
# overrides it. A deployment directory serves one machine profile (the
# first one given to the commands below).
#
# usage:
#   deploy.sh setup
#     after `tar -xf mmm-base.tar` in ROOT: extracts spack and its packages,
#     clones (or updates to the packed commit) the prototype, recreates
#     mgis/ (MGIS at its base and the patches of the prototype, branch
#     gpu-proto) and configures spack to work without network (builtin
#     repository of ROOT, bootstrap from ROOT/bootstrap).
#   deploy.sh concretize MACHINE
#     concretizes the environment deps/spack/mmm-light for the machine
#     profile MACHINE: the file spack.lock to give to `bundle.sh mirror`.
#   deploy.sh install MACHINE [MIRROR_TAR]
#     extracts the mirror of the sources made by `bundle.sh mirror` (if
#     given), then installs the environment. If PROTO_BUILDCACHE is set, the
#     installed packages and spack.lock are then copied in this binary cache
#     (a few files per package, for a filesystem with few inodes or whose
#     unused files are deleted).
#   deploy.sh cache MACHINE BUILDCACHE
#     copies the installed packages of the environment (dependencies
#     included) and its spack.lock in the binary cache BUILDCACHE (what
#     `install` does when PROTO_BUILDCACHE is set).
#   deploy.sh restore MACHINE BUILDCACHE
#     installs the environment from the binary cache of `install`, without
#     compiling, with the spack.lock kept in it.
#   deploy.sh build MACHINE [build-mfem-mgis.sh arguments]
#     builds and tests TFEL, MGIS, then mfem-mgis, the examples,
#     mm-opera-hpc and the cases (see the job build.slurm of the machine).
set -eu
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [ -f "$here/MANIFEST" ]; then
  ROOT=$here
else
  ROOT=$(cd "$here/../../.." && pwd)
fi
ROOT=${PROTO_DEPLOY:-$ROOT}
PROTO=$ROOT/MMM-GPU-Proto
export PROTO_DEPLOY=$ROOT
usage() {
  sed -n '11,36p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2
  exit 1
}
# spack of ROOT, isolated from the configuration of the user: its
# configuration is in the site scope ($ROOT/spack/etc/spack)
rspack() {
  SPACK_DISABLE_LOCAL_CONFIG=1 SPACK_USER_CACHE_PATH="$ROOT/user-cache" \
    "$ROOT/spack/bin/spack" "$@"
}
# environment of the machine profile. A deployment directory serves one
# machine profile, recorded in ROOT/MACHINE: the spack environment, its
# spack.lock and the builds of the prototype are those of this profile
machine() {
  local m=$1
  if [ ! -f "$PROTO/deps/machines/$m/machine.sh" ]; then
    echo "deploy.sh: unknown machine '$m'" >&2
    exit 1
  fi
  if [ -f "$ROOT/MACHINE" ] && [ "$(cat "$ROOT/MACHINE")" != "$m" ]; then
    echo "deploy.sh: $ROOT is the deployment of the machine" \
      "'$(cat "$ROOT/MACHINE")': deploy '$m' in another directory" >&2
    exit 1
  fi
  echo "$m" > "$ROOT/MACHINE"
  export PROTO_MACHINE=$m
  set +u
  source "$PROTO/scripts/env.sh"
  set -u
  E=$PROTO_SPACK_ENVIRONMENT
}
# environment spack of the machine profile, with the caches of ROOT
espack() {
  SPACK_USER_CACHE_PATH="$ROOT/user-cache" proto_spack -e "$E" "$@"
}
setup() {
  local t
  for t in spack spack-packages; do
    if [ ! -d "$ROOT/$t" ]; then
      echo "deploy.sh: extraction of $t"
      tar -xf "$ROOT/$t.tar" -C "$ROOT"
    fi
  done
  if [ ! -d "$PROTO/.git" ]; then
    git -c init.defaultBranch=main -c advice.detachedHead=false \
      clone --quiet "$ROOT/proto.bundle" "$PROTO"
  else
    if [ -n "$(git -C "$PROTO" status --porcelain --untracked-files=no)" ]; then
      echo "deploy.sh: $PROTO has changes which are not committed" >&2
      exit 1
    fi
    git -C "$PROTO" fetch --quiet "$ROOT/proto.bundle" HEAD
    git -C "$PROTO" -c advice.detachedHead=false checkout --quiet FETCH_HEAD
  fi
  echo "deploy.sh: prototype $(git -C "$PROTO" log -1 --format='%h %s')"
  # MGIS: base of the prototype and its patches
  if [ ! -d "$PROTO/mgis/.git" ]; then
    git -c init.defaultBranch=main init --quiet "$PROTO/mgis"
  fi
  git -C "$PROTO/mgis" fetch --quiet "$ROOT/mgis.bundle" base
  git -C "$PROTO/mgis" checkout --quiet --force -B gpu-proto FETCH_HEAD
  git -C "$PROTO/mgis" -c user.name=deploy -c user.email=deploy@localhost \
    am --quiet --keep-cr "$PROTO"/patches/mgis/*.patch
  echo "deploy.sh: mgis $(git -C "$PROTO/mgis" log -1 --format='%h %s')"
  # configuration of spack, without network
  rspack config --scope site add \
    "repos:builtin:$ROOT/spack-packages/repos/spack_repo/builtin"
  rspack config --scope site add "bootstrap:root:$ROOT/bootstrap-store"
  rspack config --scope site add "config:misc_cache:$ROOT/misc-cache"
  local s
  for s in binaries sources; do
    if ! rspack bootstrap list | grep -q "^Name: local-$s"; then
      rspack bootstrap add --scope site --trust "local-$s" \
        "$ROOT/bootstrap/metadata/$s"
    fi
  done
  for s in $(rspack bootstrap list | sed -n 's/^Name: \([^ ]*\) .*/\1/p'); do
    case "$s" in
      local-*) ;;
      *) rspack bootstrap disable --scope site "$s" ;;
    esac
  done
  rspack bootstrap now
  echo "deploy.sh: spack $(rspack --version) ready in $ROOT"
}
concretize() {
  machine "$1"
  espack concretize --force --fresh
  echo "deploy.sh: give $E/spack.lock to \`bundle.sh mirror\`"
}
install() {
  machine "$1"
  local mirror=$ROOT/mirror-$1
  if [ -n "${2:-}" ]; then
    tar -xf "$2" -C "$ROOT"
  fi
  if [ ! -d "$mirror" ]; then
    echo "deploy.sh: no mirror of the sources in $mirror" >&2
    exit 1
  fi
  if ! rspack mirror list | grep -q "^sources-$1 "; then
    rspack mirror add --scope site "sources-$1" "file://$mirror"
  fi
  # GNU patch first, in the PATH of the other builds: spack applies the
  # patches of the recipes with the command `patch`, missing on some compute
  # nodes
  espack install --fail-fast -j "$PROTO_JOBS" patch
  PATH=$(espack location -i patch)/bin:$PATH
  espack install --fail-fast -j "$PROTO_JOBS"
  if [ -n "${PROTO_BUILDCACHE:-}" ]; then
    cache "$1" "$PROTO_BUILDCACHE"
  fi
}
cache() {
  machine "$1"
  mkdir -p "$2"
  local c
  c=$(cd "$2" && pwd)
  # the packages and their dependencies, external packages excepted
  espack buildcache push --unsigned --update-index "file://$c"
  cp "$E/spack.lock" "$c/spack.lock-$1"
  echo "deploy.sh: binary cache $c ($(find "$c" -type f | wc -l) files)"
}
restore() {
  machine "$1"
  local cache
  cache=$(cd "$2" && pwd)
  cp "$cache/spack.lock-$1" "$E/spack.lock"
  if ! rspack mirror list | grep -q "^binaries-$1 "; then
    rspack mirror add --scope site --unsigned "binaries-$1" "file://$cache"
  fi
  espack install --use-buildcache only
}
build() {
  local m=$1
  shift
  machine "$m"
  "$PROTO/build-tfel.sh"
  "$PROTO/build-mgis.sh"
  "$PROTO/build-mfem-mgis.sh" "$@"
}
[ $# -ge 1 ] || usage
command=$1
shift
case "$command" in
  setup)
    [ $# -eq 0 ] || usage
    setup
    ;;
  concretize)
    [ $# -eq 1 ] || usage
    concretize "$1"
    ;;
  install)
    [ $# -ge 1 ] && [ $# -le 2 ] || usage
    install "$@"
    ;;
  cache)
    [ $# -eq 2 ] || usage
    cache "$@"
    ;;
  restore)
    [ $# -eq 2 ] || usage
    restore "$@"
    ;;
  build)
    [ $# -ge 1 ] || usage
    build "$@"
    ;;
  *)
    usage
    ;;
esac
