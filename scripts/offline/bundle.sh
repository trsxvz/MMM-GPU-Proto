#!/bin/bash
# Side of the machine with access to the internet of the installation of the
# prototype on a machine without access to the internet (docs/install.md).
# Nothing is compiled here: this script only downloads and packs.
#
# usage:
#   scripts/offline/bundle.sh base WORK
#     WORK/mmm-base.tar: spack, its packages (spack-packages), the
#     bootstrap mirror of spack (clingo and its dependencies, binaries and
#     sources), the prototype (git bundle of the current commit), MGIS at
#     the base of the patches of the prototype (git bundle), and a copy of
#     deploy.sh. The current commit of the prototype is packed: the changes
#     which are not committed are not.
#   scripts/offline/bundle.sh mirror WORK MACHINE SPACK_LOCK
#     WORK/mmm-mirror-MACHINE.tar: mirror of the sources of the packages of
#     the spack environment of MACHINE (deps/spack/mmm-light) concretized
#     on the target machine (SPACK_LOCK, given by `deploy.sh concretize`).
#     Requires the work directory WORK of `base` (instance of spack).
#
# WORK is a directory of this machine, reused from one call to the next.
set -eu
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
source "$P/scripts/offline/versions.sh"
usage() {
  sed -n '6,19p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2
  exit 1
}
# spack of the work directory, isolated from the configuration and the
# caches of the user
wspack() {
  SPACK_DISABLE_LOCAL_CONFIG=1 SPACK_USER_CACHE_PATH="$W/user-cache" \
    "$W/spack/bin/spack" "$@"
}
# clone (or update) of a repository at a given commit
checkout() {
  local url=$1 directory=$2 commit=$3
  if [ ! -d "$directory/.git" ]; then
    git clone --quiet "$url" "$directory"
  fi
  git -C "$directory" fetch --quiet --tags origin
  git -C "$directory" -c advice.detachedHead=false checkout --quiet "$commit"
}
base() {
  echo "bundle.sh: spack $PROTO_SPACK_VERSION"
  checkout https://github.com/spack/spack.git "$W/spack" "$PROTO_SPACK_VERSION"
  echo "bundle.sh: spack-packages $PROTO_SPACK_PACKAGES_COMMIT"
  checkout https://github.com/spack/spack-packages.git "$W/spack-packages" \
    "$PROTO_SPACK_PACKAGES_COMMIT"
  # builtin repository of this instance: the clone above, not the one
  # spack downloads by default
  wspack config --scope site add \
    "repos:builtin:$W/spack-packages/repos/spack_repo/builtin"
  echo "bundle.sh: bootstrap mirror of spack"
  wspack bootstrap mirror --binary-packages "$W/bootstrap"
  echo "bundle.sh: MGIS $PROTO_MGIS_BASE"
  checkout "$PROTO_MGIS_URL" "$W/mgis" "$PROTO_MGIS_BASE"
  git -C "$W/mgis" branch --force base "$PROTO_MGIS_BASE"
  git -C "$W/mgis" bundle create --quiet "$W/stage/mgis.bundle" base
  local commit
  commit=$(git -C "$P" rev-parse --short HEAD)
  echo "bundle.sh: prototype $commit"
  if [ -n "$(git -C "$P" status --porcelain --untracked-files=no)" ]; then
    echo "bundle.sh: warning: the changes of the prototype which are not" \
      "committed are not packed" >&2
  fi
  git -C "$P" bundle create --quiet "$W/stage/proto.bundle" HEAD
  cp "$P/scripts/offline/deploy.sh" "$W/stage/deploy.sh"
  cat > "$W/stage/MANIFEST" << EOF
prototype $commit
spack $PROTO_SPACK_VERSION
spack-packages $PROTO_SPACK_PACKAGES_COMMIT
mgis $PROTO_MGIS_BASE
created $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname)
EOF
  # spack and spack-packages without their history
  git -C "$W/spack" archive --prefix=spack/ HEAD > "$W/stage/spack.tar"
  git -C "$W/spack-packages" archive --prefix=spack-packages/ HEAD \
    > "$W/stage/spack-packages.tar"
  tar -cf "$W/mmm-base.tar" -C "$W/stage" MANIFEST deploy.sh proto.bundle \
    mgis.bundle spack.tar spack-packages.tar -C "$W" bootstrap
  echo "bundle.sh: $W/mmm-base.tar ($(du -h "$W/mmm-base.tar" | cut -f1))"
}
mirror() {
  local machine=$1 lock=$2
  if [ ! -x "$W/spack/bin/spack" ]; then
    echo "bundle.sh: no instance of spack in $W (run \`bundle.sh base\`)" >&2
    exit 1
  fi
  if [ ! -f "$P/deps/machines/$machine/spack-config.yaml" ]; then
    echo "bundle.sh: unknown machine '$machine'" >&2
    exit 1
  fi
  # copy of the environment with the lock file of the target machine; the
  # repository mmm_gpu is found at $env/../repos
  local e=$W/environments/$machine
  mkdir -p "$e"
  ln -sfn "$P/deps/spack/repos" "$W/environments/repos"
  cp "$P/deps/spack/mmm-light/spack.yaml" "$e/spack.yaml"
  cp "$lock" "$e/spack.lock"
  export PROTO_MACHINE_DIR=$P/deps/machines/$machine
  export PROTO_SPACK_PACKAGES=$W/spack-packages
  local m=$W/mirror-$machine
  wspack -e "$e" mirror create --directory "$m" --all
  tar -cf "$W/mmm-mirror-$machine.tar" -C "$W" "mirror-$machine"
  echo "bundle.sh: $W/mmm-mirror-$machine.tar" \
    "($(du -h "$W/mmm-mirror-$machine.tar" | cut -f1))"
}
command=${1:-}
[ $# -ge 2 ] || usage
W=$(mkdir -p "$2" && cd "$2" && pwd)
mkdir -p "$W/stage"
case "$command" in
  base)
    [ $# -eq 2 ] || usage
    base
    ;;
  mirror)
    [ $# -eq 4 ] || usage
    mirror "$3" "$(cd "$(dirname "$4")" && pwd)/$(basename "$4")"
    ;;
  *)
    usage
    ;;
esac
