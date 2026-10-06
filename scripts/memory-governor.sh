#!/bin/bash
# Pauses a process tree (SIGSTOP) when the available memory falls below a
# low mark, and resumes it (SIGCONT) above a high mark, so that a long
# background build does not lead the kernel to kill other processes.
# A paused process keeps its memory: this only helps against transient
# peaks. A build whose processes need a lot of memory (for instance the
# CUDA kernels of MFEM, about 10 GB per ptxas process) must rather be run
# with fewer parallel jobs, or alone.
# usage: scripts/memory-governor.sh root-pid [low MiB] [high MiB]
root=$1
low=${2:-3500}
high=${3:-7000}
paused=0
tree() {
  echo "$1"
  for c in $(ps -o pid= --ppid "$1"); do
    tree "$c"
  done
}
while kill -0 "$root" 2> /dev/null; do
  available=$(awk '/MemAvailable/ {print int($2 / 1024)}' /proc/meminfo)
  if [ "$paused" = 0 ] && [ "$available" -lt "$low" ]; then
    kill -STOP $(tree "$root") 2> /dev/null
    paused=1
    echo "$(date +%H:%M:%S) paused (${available} MiB available)"
  elif [ "$paused" = 1 ] && [ "$available" -gt "$high" ]; then
    kill -CONT $(tree "$root") 2> /dev/null
    paused=0
    echo "$(date +%H:%M:%S) resumed (${available} MiB available)"
  fi
  sleep 2
done
