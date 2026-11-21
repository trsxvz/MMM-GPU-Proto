#!/bin/bash
# Cost per call of the handling of the statistics of a synchronous batch
# without caller statistics, for several designs (see the header of
# statistics-blocks.cu): trivial kernel alone, copies before and after the
# kernel, read back, flag in pinned memory, managed memory; alone, then
# after 200 µs of work enqueued on the same stream. A last check gives the
# number of failures of three successive calls (8, 0, 8 expected): the
# variant V2 gives wrong counts, its reset by a copy on the default stream
# being not ordered with the non-blocking stream of the calls (a block must
# be reset synchronously, on the stream of the call).
# usage: reproducers/statistics-blocks/run.sh
set -e
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
source "$P/scripts/env.sh" > /dev/null 2>&1
W=$P/scratch/reproducers/statistics-blocks
rm -rf "$W" && mkdir -p "$W" && cd "$W"
nvcc -O3 -std=c++20 -arch=native "$R/statistics-blocks.cu" -o statistics-blocks
./statistics-blocks
