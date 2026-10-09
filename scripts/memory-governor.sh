#!/bin/bash
# Pauses a process tree (SIGSTOP) when the available memory falls below a
# low mark, and resumes it (SIGCONT) above a high mark, so that a long
# background build does not lead the kernel to kill other processes.
#
# A paused process keeps its memory. So the process of the tree using the
# most memory, and its ancestors, keep running: it ends and frees its
# memory, and the other ones are resumed. Pausing it too would wait forever
# (a ptxas process of MFEM takes about 22 GB on its own). If the memory is
# still low after `timeout` seconds, the tree is resumed anyway.
#
# The paused processes are the ones recorded at the pause: they are resumed
# even if they left the tree (parent killed), and when the root or this
# script ends. A paused process with a pending signal (for instance SIGTERM
# sent by spack) is resumed at once, so that it receives it.
# usage: scripts/memory-governor.sh root-pid [low MiB] [high MiB] [timeout s]
root=$1
low=${2:-3500}
high=${3:-7000}
timeout=${4:-1800}
paused=0
since=0
# processes paused by this script
stopped=""
resume() {
  if [ -n "$stopped" ]; then
    kill -CONT $stopped 2> /dev/null
  fi
  stopped=""
  paused=0
}
trap resume EXIT
trap 'exit 1' INT TERM HUP
tree() {
  echo "$1"
  for c in $(ps -o pid= --ppid "$1"); do
    tree "$c"
  done
}
ancestors() {
  local p=$1
  while [ -n "$p" ] && [ "$p" != 0 ] && [ "$p" != "$root" ]; do
    echo "$p"
    p=$(ps -o ppid= -p "$p" | tr -d ' ')
  done
  echo "$root"
}
while kill -0 "$root" 2> /dev/null; do
  available=$(awk '/MemAvailable/ {print int($2 / 1024)}' /proc/meminfo)
  if [ "$paused" = 0 ] && [ "$available" -lt "$low" ]; then
    processes=$(tree "$root")
    largest=$(ps -o pid=,rss= -p "$(echo $processes | tr ' ' ',')" 2> /dev/null |
      sort -k2 -n -r | awk 'NR == 1 {print $1}')
    keep=" $(ancestors "$largest" | tr '\n' ' ') "
    for p in $processes; do
      case "$keep" in
        *" $p "*) ;;
        *) kill -STOP "$p" 2> /dev/null && stopped="$stopped $p" ;;
      esac
    done
    paused=1
    since=$(date +%s)
    echo "$(date +%H:%M:%S) paused all but $largest (${available} MiB available)"
  elif [ "$paused" = 1 ]; then
    for p in $stopped; do
      if awk '/^(SigPnd|ShdPnd):/ && $2 !~ /^0+$/ {f = 1} END {exit !f}' \
        "/proc/$p/status" 2> /dev/null; then
        kill -CONT "$p" 2> /dev/null
      fi
    done
  fi
  if [ "$paused" = 1 ] && { [ "$available" -gt "$high" ] ||
    [ $(($(date +%s) - since)) -gt "$timeout" ]; }; then
    resume
    echo "$(date +%H:%M:%S) resumed (${available} MiB available)"
  fi
  sleep 2
done
