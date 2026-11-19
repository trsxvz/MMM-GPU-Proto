#!/bin/bash
# Gathers the files of one time step of a ParaView export (directory
# CycleNNNNNN) in the archive CycleNNNNNN.tar next to it, then removes the
# directory, once the archive is checked (as many entries as files, from the
# listing written while archiving: the archive is not read again). Records of
# 1 MiB (-b 2048): the default records of 10 KiB make tar very slow on Lustre.
# usage: archive-cycle.sh <directory CycleNNNNNN>
set -u
cycle=${1%/}
[ -d "$cycle" ] || exit 0
name=$(basename "$cycle")
dir=$(dirname "$cycle")
listing=$(tar -b 2048 -cvf "$cycle.tar.part" -C "$dir" "$name") || {
  rm -f "$cycle.tar.part"
  echo "archive of $cycle failed, its files are kept" >&2
  exit 1
}
if [ "$(grep -vc '/$' <<< "$listing")" = "$(find "$cycle" -type f | wc -l)" ]; then
  mv "$cycle.tar.part" "$cycle.tar" && rm -rf "$cycle"
else
  rm -f "$cycle.tar.part"
  echo "archive of $cycle incomplete, its files are kept" >&2
  exit 1
fi
