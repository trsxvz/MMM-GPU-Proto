#!/bin/bash
# Build diagonal-rows.cxx and structures.cxx with the compiler and the flags
# of the MFEM of the prototype (config.mk of MFEM, nvcc) and run them on one,
# two and three processes:
# - diagonal-rows: tridiagonal rows, except the second half of the rows of
#   the last process, which only have a diagonal coefficient, for a coupled
#   or block diagonal matrix, with a symmetric pattern or not, with 4 and
#   400 rows per process;
# - structures: tridiagonal matrix coupling the processes, and matrices
#   whose rows are coupled with rows of the next process, distributed
#   evenly or not.
# usage: reproducers/cudss-distributed/run.sh
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
source "$P/scripts/env.sh" > /dev/null 2>&1
W=$P/scratch/reproducers/cudss-distributed
mkdir -p "$W" && cd "$W" || exit 1
{
  printf 'include %s/share/mfem/config.mk\n\n' "$PROTO_SPACK_VIEW"
  for p in diagonal-rows structures; do
    printf '%s: %s/%s.cxx\n\t$(MFEM_CXX) $(MFEM_FLAGS) -x cu $< -o $@ $(MFEM_LIBS)\n\n' \
      "$p" "$R" "$p"
  done
} > Makefile
make -s diagonal-rows structures > build.log 2>&1 || { cat build.log; exit 1; }
export OMP_NUM_THREADS=1
result() {
  timeout 60 mpirun -n "$@" 2>&1 | grep -aE "relative|code:|Segmentation" |
    head -1
}
for np in 1 2; do
  for coupled in 1 0; do
    for symmetric in 1 0; do
      for n in 4 400; do
        printf "diagonal rows, np %s, coupled %s, symmetric %s, n %s: %s\n" \
          "$np" "$coupled" "$symmetric" "$n" \
          "$(result "$np" ./diagonal-rows "$coupled" "$symmetric" "$n")"
      done
    done
  done
done
for rows in "12" "4 4" "4 4 4" "40 40 40"; do
  set -- $rows
  printf "tridiagonal, rows %s: %s\n" "$rows" \
    "$(result $# ./structures tridiagonal $rows)"
done
for rows in "4 4" "4 6" "4 4 4" "4 4 6" "4 3 5" "4 2 6" "5 2 5"; do
  set -- $rows
  printf "next, rows %s: %s\n" "$rows" "$(result $# ./structures next $rows)"
done
