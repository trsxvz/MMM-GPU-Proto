#!/bin/bash
# usage (on Adastra): deploy-overlay.sh <name> <archive of cases/ files>
# Builds /lus/scratch/.../proto/<name>: links to every entry of the deployed
# prototype (/lus/work, whose quota of files is reached), except cases/, a
# directory of links to its files, the files of the archive replacing them.
# The jobs use it with FM_PROTO=<that directory>: the scripts of a running
# job are never overwritten (bash reads a script while it runs it).
set -eu
src=/lus/work/RES1/genden15/tchenaille/mmm/MMM-GPU-Proto
dst=/lus/scratch/RES1/genden15/tchenaille/proto/$1
[ -e "$dst" ] && { echo "$dst exists" >&2; exit 1; }
mkdir -p "$dst/cases"
for e in "$src"/* "$src"/.git; do
  b=$(basename "$e")
  [ "$b" = cases ] || ln -s "$e" "$dst/$b"
done
for e in "$src"/cases/*; do ln -s "$e" "$dst/cases/$(basename "$e")"; done
tar xzf "$2" -C "$dst/cases" --unlink-first
echo "$dst"
